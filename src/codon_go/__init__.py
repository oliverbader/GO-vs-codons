"""
Codon-GO Analysis Pipeline

A modular Python pipeline for analyzing codon usage and GO term enrichment
in eukaryotic genomes.
"""

from .utils.warnings_config import suppress_common_warnings

# Suppress common warnings on import
suppress_common_warnings()

__version__ = "1.0.0"
__author__ = "Codon-GO Pipeline"
__email__ = "codon-go@example.com"

from .parsers import genome_parser, go_parser
from .analysis import codon_usage, stats
from .viz import boxplots, heatmap, pca_scatter
from .utils import config_loader, file_utils