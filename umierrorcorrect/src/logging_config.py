#!/usr/bin/env python3
"""Logging configuration using loguru for umierrorcorrect."""

import sys
from pathlib import Path
from typing import Literal

from loguru import logger

# Default log format
LOG_FORMAT = (
    "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
    "<level>{message}</level>"
)

# Simple format for console output
SIMPLE_FORMAT = "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{message}</level>"


def setup_logging(
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO",
    log_file: str | Path | None = None,
    simple_format: bool = True,
    rotation: str = "10 MB",
    retention: str = "7 days",
) -> None:
    """Configure loguru logging for the application.

    Args:
        level: Minimum log level to display.
        log_file: Optional path to log file. If None, logs only to console.
        simple_format: Use simple format for console (default True).
        rotation: When to rotate log file (e.g., "10 MB", "1 day").
        retention: How long to keep old log files.
    """
    # Remove default handler
    logger.remove()

    # Console handler with colored output
    log_format = SIMPLE_FORMAT if simple_format else LOG_FORMAT
    logger.add(
        sys.stderr,
        format=log_format,
        level=level,
        colorize=True,
    )

    # File handler if specified
    if log_file:
        log_path = Path(log_file)
        logger.add(
            str(log_path),
            format=LOG_FORMAT,
            level=level,
            rotation=rotation,
            retention=retention,
            compression="gz",
        )


def get_logger(name: str | None = None):
    """Get a logger instance.

    Args:
        name: Optional name for the logger context.

    Returns:
        Configured logger instance.
    """
    if name:
        return logger.bind(name=name)
    return logger


# Intercept standard logging to redirect to loguru
class InterceptHandler:
    """Handler to intercept standard logging and redirect to loguru."""

    def __init__(self, level: int = 0) -> None:
        self.level = level

    def emit(self, record) -> None:
        """Emit a log record."""
        # Get corresponding Loguru level if it exists
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # Find caller from where originated the logged message
        frame, depth = sys._getframe(6), 6
        while frame and frame.f_code.co_filename == __file__:
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())


def redirect_standard_logging() -> None:
    """Redirect standard logging module to use loguru."""
    import logging

    logging.basicConfig(handlers=[InterceptHandler()], level=0, force=True)
