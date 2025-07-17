"""
Warning configuration utility for the Codon-GO pipeline.
"""

import warnings
import logging

logger = logging.getLogger(__name__)


def suppress_common_warnings():
    """
    Suppress common warnings that clutter the output without being actionable.
    """
    # Suppress pkg_resources deprecation warning
    warnings.filterwarnings(
        "ignore", 
        message="pkg_resources is deprecated", 
        category=UserWarning
    )
    
    # Suppress setuptools deprecation warnings
    warnings.filterwarnings(
        "ignore",
        message=".*pkg_resources.*",
        category=DeprecationWarning
    )
    
    # Suppress specific matplotlib warnings that are common but not actionable
    warnings.filterwarnings(
        "ignore",
        message=".*tight_layout.*",
        category=UserWarning
    )
    
    # Suppress pandas future warnings that are not immediately actionable
    warnings.filterwarnings(
        "ignore",
        message=".*DataFrame.applymap.*",
        category=FutureWarning
    )
    
    logger.debug("Common warnings suppressed")


def configure_warnings(verbose: bool = False):
    """
    Configure warning behavior based on verbosity level.
    
    Args:
        verbose: If True, show more warnings; if False, suppress common ones
    """
    if not verbose:
        suppress_common_warnings()
    else:
        # In verbose mode, show warnings but format them nicely
        warnings.formatwarning = lambda message, category, filename, lineno, file=None, line=None: \
            f"Warning ({category.__name__}): {message}\n"
        
        logger.debug("Verbose warning mode enabled")


def with_warnings_suppressed(func):
    """
    Decorator to suppress warnings for a specific function.
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function with warnings suppressed
    """
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            suppress_common_warnings()
            return func(*args, **kwargs)
    
    return wrapper