import logging
from datetime import datetime, timedelta
from typing import Any, Dict, Optional, Tuple, Union

from jose import ExpiredSignatureError, JWTError, jwt
from pydantic import ValidationError

from studio.app.common.core.auth.auth_config import AUTH_CONFIG


def _create_token(
    subject: Union[str, Any],
    token_type: str,
    expires_delta: timedelta = None,
    user_claims: Dict[str, Any] = None,
) -> str:
    token_data = {"sub": subject, "token_type": token_type}

    token_data.update(
        {"exp": datetime.utcnow() + expires_delta} if expires_delta else {}
    )
    token_data.update(user_claims if user_claims else {})

    return jwt.encode(
        token_data,
        AUTH_CONFIG.SECRET_KEY,
        algorithm=AUTH_CONFIG.ALGORITHM,
    )


def create_access_token(
    subject: Union[str, Any],
    expires_delta: timedelta = None,
    user_claims: Dict[str, Any] = None,
) -> str:
    return _create_token(
        subject,
        token_type="access_token",
        expires_delta=expires_delta,
        user_claims=user_claims,
    )


def create_refresh_token(
    subject: Union[str, Any],
    expires_delta: timedelta = None,
    user_claims: Dict[str, Any] = None,
) -> str:
    exp = (
        expires_delta
        if expires_delta
        else timedelta(minutes=AUTH_CONFIG.REFRESH_TOKEN_EXPIRE_MINUTES)
    )
    return _create_token(
        subject,
        token_type="refresh_token",
        expires_delta=exp,
        user_claims=user_claims,
    )


def _validate_token(
    token: str,
    token_type: str,
) -> Tuple[Dict[str, Any], Optional[str]]:
    payload = err = None
    try:
        payload = jwt.decode(
            token, AUTH_CONFIG.SECRET_KEY, algorithms=[AUTH_CONFIG.ALGORITHM]
        )
        if payload.get("token_type") != token_type:
            err = "Invalid token type"
    except ExpiredSignatureError as e:
        logging.getLogger().error(e)
        err = "Credentials has expired"
    except (JWTError, ValidationError) as e:
        logging.getLogger().error(e)
        err = "Could not validate credentials"
    except Exception as e:
        logging.getLogger().error(e)
        err = "Bad token"

    return payload, err


def validate_access_token(token: str) -> Tuple[Dict[str, Any], Optional[str]]:
    return _validate_token(token, "access_token")


def validate_refresh_token(token: str) -> Tuple[Dict[str, Any], Optional[str]]:
    return _validate_token(token, "refresh_token")
