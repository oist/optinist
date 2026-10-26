"""
Logging middleware for capturing user context in logs
"""

from typing import Optional

from fastapi import Request
from starlette.types import ASGIApp, Receive, Scope, Send

from studio.app.common.core.auth.auth_helper import extract_uid_from_request_async
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.mode import MODE


class ClientIdLoggingMiddleware:
    """
    Pure ASGI Middleware to extract client_id from request and set it in logging context

    This middleware:
    - Extracts uid from authentication tokens (Firebase or JWT (ExToken))
    - Generate client_id from uid (by hashing, etc.)
    - Sets the client_id in the logging context for all log messages during the request
    - In standalone mode, skips client_id extraction and logging (client_id is blank)

    Note: This uses pure ASGI middleware instead of BaseHTTPMiddleware to avoid
    known performance issues and connection handling problems with BaseHTTPMiddleware.
    See: https://github.com/encode/starlette/issues/1012
    """

    def __init__(self, app: ASGIApp) -> None:
        self.app = app

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        if scope["type"] != "http":
            # Skip non-HTTP requests (e.g., websocket, lifespan)
            # Note: scope["type"] == "http" for both HTTP and HTTPS requests
            await self.app(scope, receive, send)
            return

        # Extract client_id from request headers
        client_id: Optional[str] = None

        # Skip uid extraction for standalone mode
        if not MODE.IS_STANDALONE:
            request = Request(scope, receive)
            # Use async version with caching and thread pool for better performance
            uid = await extract_uid_from_request_async(request)
            client_id = AppLogger.generate_client_id(uid)

        # Set client_id in logging context
        AppLogger.set_client_id(client_id)

        # Process the request
        await self.app(scope, receive, send)
