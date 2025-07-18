"""
Analysis module for codon usage and statistical analysis.
"""

from .codon_usage import compute_relative_usage_by_aa, get_codon_usage_stats
from .stats import adaptive_go_analysis_by_codon, perform_enrichment_test