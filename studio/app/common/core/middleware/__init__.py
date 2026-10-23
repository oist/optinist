"""
Middleware modules for the Studio application
"""

from studio.app.common.core.middleware.logging_middleware import (
    ClientIdLoggingMiddleware,
)

__all__ = ["ClientIdLoggingMiddleware"]
