#!/usr/bin/env python3
"""
Main entry point for the Codon-GO Analysis Pipeline.
"""

import sys
import os

# Add src directory to path to allow imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Import and configure warnings suppression
from codon_go.utils.warnings_config import suppress_common_warnings
suppress_common_warnings()

from codon_go.cli import cli

if __name__ == '__main__':
    cli()