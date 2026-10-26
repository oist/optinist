"""
Authentication helper functions for extracting user identity from tokens

This module provides shared utilities for extracting user_id (uid) from
authentication tokens, supporting both Firebase and JWT authentication methods.
"""

import asyncio
import hashlib
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Optional, Tuple

from fastapi import Request
from fastapi.security import HTTPAuthorizationCredentials
from firebase_admin import auth as firebase_auth

from studio.app.common.core.auth.auth_config import AUTH_CONFIG
from studio.app.common.core.auth.security import validate_access_token

# Thread pool for running synchronous Firebase calls asynchronously
_thread_pool = ThreadPoolExecutor(max_workers=10, thread_name_prefix="firebase_auth")

# Token cache with TTL (Time To Live)
# Structure: {token_hash: {"uid": str, "expires_at": float}}
#
# IMPORTANT: Multi-process considerations
# ========================================
# This is an in-memory, process-local cache. In multi-process deployments
# (e.g., Gunicorn/Uvicorn with multiple workers), each process maintains
# its own independent cache.
#
# Implications:
# - Security: No security issues. Each process correctly validates tokens.
# - Performance: Cache hit rate decreases proportionally to the number of processes.
#   For example, with 4 workers, the effective cache hit rate is ~25% of single-process.
# - Firebase API calls: Increases by a factor of the number of processes.
#
# Alternative solutions for multi-process environments:
# - Use a shared cache backend (Redis, Memcached) for cross-process cache sharing
# - Rely on Firebase Admin SDK's built-in caching (may already provide some caching)
# - Use sticky sessions at the load balancer to improve cache hit rates
#
# Current decision: Acceptable for current deployment as Firebase API has
#  sufficient quota and the performance impact is manageable.
_token_cache: Dict[str, Dict[str, any]] = {}

# Cache TTL in seconds (5 minutes)
_CACHE_TTL_SECONDS = 300


def _compute_token_hash(token: str) -> str:
    """
    Compute a hash of the token for use as a cache key

    Args:
        token: Authentication token string

    Returns:
        SHA256 hash of the token
    """
    return hashlib.sha256(token.encode()).hexdigest()


def _get_cached_uid(token: str) -> Optional[str]:
    """
    Get cached uid for a token if it exists and hasn't expired

    Args:
        token: Authentication token string

    Returns:
        Cached uid if found and valid, None otherwise
    """
    token_hash = _compute_token_hash(token)
    cached_data = _token_cache.get(token_hash)

    if cached_data is None:
        return None

    # Check if cache has expired
    if time.time() > cached_data["expires_at"]:
        # Remove expired cache entry
        del _token_cache[token_hash]
        return None

    return cached_data["uid"]


def _cache_uid(token: str, uid: Optional[str]) -> None:
    """
    Cache the uid for a token with TTL

    Args:
        token: Authentication token string
        uid: User ID to cache
    """
    if uid is None:
        return

    token_hash = _compute_token_hash(token)
    _token_cache[token_hash] = {
        "uid": uid,
        "expires_at": time.time() + _CACHE_TTL_SECONDS,
    }


def _verify_firebase_token_sync(token: str) -> Optional[str]:
    """
    Synchronous Firebase token verification (for use in thread pool)

    Args:
        token: Firebase ID token

    Returns:
        User ID if verification succeeds, None otherwise
    """
    try:
        user = firebase_auth.verify_id_token(token)
        return user.get("uid")
    except Exception:
        return None


async def _verify_firebase_token_async(token: str) -> Optional[str]:
    """
    Asynchronous Firebase token verification using thread pool

    Args:
        token: Firebase ID token

    Returns:
        User ID if verification succeeds, None otherwise
    """
    # Check cache first
    cached_uid = _get_cached_uid(token)
    if cached_uid is not None:
        return cached_uid

    # Run synchronous Firebase call in thread pool
    loop = asyncio.get_event_loop()
    uid = await loop.run_in_executor(_thread_pool, _verify_firebase_token_sync, token)

    # Cache the result
    _cache_uid(token, uid)

    return uid


def extract_uid_from_firebase_credential(
    credential: HTTPAuthorizationCredentials,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract user_id (uid) from Firebase authentication credential

    Args:
        credential: HTTPAuthorizationCredentials from FastAPI

    Returns:
        Tuple of (uid, error_message). If successful, error_message is None.
    """
    try:
        user = firebase_auth.verify_id_token(credential.credentials)
        uid = user.get("uid")
        return uid, None
    except Exception as e:
        return None, f"Firebase token validation failed: {e}"


def extract_uid_from_jwt_token(ex_token: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract user_id (uid) from JWT ExToken

    Args:
        ex_token: JWT token string

    Returns:
        Tuple of (uid, error_message). If successful, error_message is None.
    """
    try:
        payload, err = validate_access_token(ex_token)
        if err is not None:
            return None, err
        if payload:
            uid = payload.get("sub")
            return uid, None
        return None, "Empty payload"
    except Exception as e:
        return None, f"JWT token validation failed: {e}"


def _extract_token_from_request(request: Request) -> Optional[str]:
    """
    Extract authentication token from request headers based on configured auth method

    This function checks AUTH_CONFIG.USE_FIREBASE_TOKEN to determine which
    authentication method to use and extracts the appropriate token.

    Args:
        request: FastAPI Request object

    Returns:
        Token string if found, None otherwise
    """
    if AUTH_CONFIG.USE_FIREBASE_TOKEN:
        # Firebase authentication via Bearer token
        auth_header = request.headers.get("Authorization")
        if auth_header and auth_header.startswith("Bearer "):
            return auth_header.replace("Bearer ", "")
    else:
        # JWT authentication via ExToken header
        ex_token = request.headers.get("ExToken")
        if ex_token:
            return ex_token

    return None


def extract_uid_from_request(request: Request) -> Optional[str]:
    """
    Extract user_id (uid) from HTTP request based on configured authentication method

    This is a SYNCHRONOUS convenience function for middleware and other non-dependency
    contexts where you need to extract uid without raising exceptions.

    Note: This function uses caching to reduce Firebase API calls.
    For async contexts, prefer extract_uid_from_request_async() for better performance.

    Note: This function directly accesses request.headers because it's designed for
    middleware context where FastAPI's Dependency Injection is not available.
    For route handlers, use extract_uid_from_firebase_credential() or
    extract_uid_from_jwt_token() with FastAPI's Depends(APIKeyHeader/HTTPBearer).

    Args:
        request: FastAPI Request object

    Returns:
        user_id if found, None otherwise (errors are silently ignored)
    """
    try:
        token = _extract_token_from_request(request)
        if token is None:
            return None

        if AUTH_CONFIG.USE_FIREBASE_TOKEN:
            # Firebase authentication
            # Check cache first
            cached_uid = _get_cached_uid(token)
            if cached_uid is not None:
                return cached_uid

            # Verify token and cache result
            try:
                user = firebase_auth.verify_id_token(token)
                uid = user.get("uid")
                _cache_uid(token, uid)
                return uid
            except Exception:
                pass

        else:
            # JWT authentication
            try:
                payload, err = validate_access_token(token)
                if err is None and payload:
                    return payload.get("sub")
            except Exception:
                pass

    except Exception:
        pass

    return None


async def extract_uid_from_request_async(request: Request) -> Optional[str]:
    """
    ASYNC version: Extract user_id (uid) from HTTP request
      based on configured auth method

    This is an ASYNCHRONOUS convenience function that should be used in async contexts
    like middleware. It runs Firebase token verification in a thread pool to avoid
    blocking the event loop.

    Benefits over sync version:
    - Non-blocking Firebase API calls (runs in thread pool)
    - Automatic caching with TTL
    - Better performance for high-concurrency scenarios

    Recommended usage: Use this in all async contexts
      (middleware, async route handlers, etc.)

    Args:
        request: FastAPI Request object

    Returns:
        user_id if found, None otherwise (errors are silently ignored)
    """
    try:
        token = _extract_token_from_request(request)
        if token is None:
            return None

        if AUTH_CONFIG.USE_FIREBASE_TOKEN:
            # Firebase authentication
            # Use async verification with caching and thread pool
            uid = await _verify_firebase_token_async(token)
            return uid

        else:
            # JWT authentication
            # JWT validation is typically fast,
            #   run in thread pool anyway to avoid blocking
            try:
                loop = asyncio.get_event_loop()
                payload, err = await loop.run_in_executor(
                    _thread_pool, validate_access_token, token
                )
                if err is None and payload:
                    return payload.get("sub")
            except Exception:
                pass

    except Exception:
        pass

    return None
