"""
Parsers module for genome and GO data processing.
"""

from .genome_parser import load_genome_annotations
from .go_parser import load_go_ontology, parse_gaf_file