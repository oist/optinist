"""
Helper utilities for managing logging context across process boundaries

This module provides utilities to propagate client_id context across:
- ProcessPoolExecutor subprocess boundaries
- Snakemake script execution contexts
"""

from functools import wraps
from typing import Any, Callable, Optional, TypeVar

from studio.app.common.core.logger import LOGGING_CLIENT_ID_KEY, AppLogger

R = TypeVar("R")


def with_client_id_context(func: Callable[..., R]) -> Callable[..., R]:
    """
    Decorator to automatically set client_id
      in logging context for subprocess functions.

    This decorator extracts 'client_id' from function kwargs and sets it in the
    logging context before executing the function. Useful for ProcessPoolExecutor
    subprocess functions.

    Usage:
        @with_client_id_context
        def subprocess_function(arg1, arg2):
            # client_id is automatically set in logging context
            logger.info("This log will include client_id")

        # Call from parent process
        client_id = get_client_id_for_subprocess()
        executor.submit(subprocess_function, arg1, arg2, client_id=client_id)

    Args:
        func: The function to wrap

    Returns:
        Wrapped function with client_id context management
    """

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> R:
        # Get client_id from kwargs (keep it in kwargs for function to use)
        client_id = kwargs.get(LOGGING_CLIENT_ID_KEY, None)

        # Set in logging context
        if client_id is not None:
            AppLogger.set_client_id(client_id)

        # Execute original function (client_id remains in kwargs)
        return func(*args, **kwargs)

    return wrapper


def get_client_id_for_subprocess() -> Optional[str]:
    """
    Get client_id from current context for passing to subprocess.

    Use this function in the parent process to retrieve the current client_id
    before spawning a subprocess via ProcessPoolExecutor.

    Usage:
        client_id = get_client_id_for_subprocess()
        future = executor.submit(subprocess_func, args, client_id=client_id)

    Returns:
        Current client_id or None if not set
    """
    return AppLogger.get_client_id()


def init_client_id_from_snakemake_config(config: dict):
    """
    Initialize client_id from snakemake config.

    This function should be called at the beginning of the main() function
    in snakemake script files to set up the logging context with client_id
    passed through snakemake config.

    Args:
        config: Snakemake config dictionary (typically snakemake.config)

    Usage (in snakemake script main function):
        def main():
            init_client_id_from_snakemake_config(snakemake.config)
            # ... rest of processing
    """
    client_id = config.get(LOGGING_CLIENT_ID_KEY)
    AppLogger.set_client_id(client_id)
