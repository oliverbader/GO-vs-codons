"""
GO parser module for loading GO ontology and gene annotations.
"""

import os
from typing import Dict, List, Set, Tuple
import pandas as pd
from goatools import obo_parser
from goatools.associations import read_gaf
import logging

logger = logging.getLogger(__name__)


def load_go_ontology(obo_file: str) -> obo_parser.GODag:
    """
    Load GO ontology from OBO file using GOATOOLS.
    
    Args:
        obo_file: Path to GO ontology OBO file
        
    Returns:
        GODag object containing the ontology
    """
    if not os.path.exists(obo_file):
        raise FileNotFoundError(f"GO ontology file not found: {obo_file}")
    
    logger.info(f"Loading GO ontology from {obo_file}")
    
    try:
        go_dag = obo_parser.GODag(obo_file)
        logger.info(f"Loaded {len(go_dag)} GO terms")
        return go_dag
    except Exception as e:
        logger.error(f"Error loading GO ontology: {e}")
        raise


def parse_gaf_file(gaf_file: str) -> pd.DataFrame:
    """
    Parse GAF (Gene Association File) to extract gene-GO mappings.
    
    Args:
        gaf_file: Path to GAF file
        
    Returns:
        DataFrame with columns: gene_id, go_id, evidence, aspect
    """
    if not os.path.exists(gaf_file):
        raise FileNotFoundError(f"GAF file not found: {gaf_file}")
    
    logger.info(f"Parsing GAF file: {gaf_file}")
    
    try:
        # Read GAF file using GOATOOLS
        associations = read_gaf(gaf_file)
        
        rows = []
        for gene_id, go_terms in associations.items():
            for go_term in go_terms:
                rows.append({
                    'gene_id': gene_id,
                    'go_id': go_term.id,
                    'evidence': getattr(go_term, 'evidence_code', 'Unknown'),
                    'aspect': getattr(go_term, 'aspect', 'Unknown')
                })
        
        df = pd.DataFrame(rows)
        logger.info(f"Parsed {len(df)} gene-GO associations for {df['gene_id'].nunique()} genes")
        return df
        
    except Exception as e:
        logger.error(f"Error parsing GAF file: {e}")
        # Fallback to manual parsing
        return _parse_gaf_manual(gaf_file)


def _parse_gaf_manual(gaf_file: str) -> pd.DataFrame:
    """
    Manual GAF parsing as fallback.
    
    Args:
        gaf_file: Path to GAF file
        
    Returns:
        DataFrame with gene-GO mappings
    """
    rows = []
    
    with open(gaf_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if line.startswith('!') or not line:
                continue
            
            fields = line.split('\t')
            
            # GAF format has at least 15 fields
            if len(fields) < 15:
                continue
            
            try:
                gene_id = fields[1]  # DB_Object_ID
                go_id = fields[4]    # GO_ID
                evidence = fields[6] # Evidence_Code
                aspect = fields[8]   # Aspect
                
                rows.append({
                    'gene_id': gene_id,
                    'go_id': go_id,
                    'evidence': evidence,
                    'aspect': aspect
                })
                
            except IndexError:
                continue
    
    df = pd.DataFrame(rows)
    logger.info(f"Manually parsed {len(df)} gene-GO associations")
    return df


def filter_go_associations(df: pd.DataFrame, 
                          evidence_codes: List[str] = None,
                          aspects: List[str] = None) -> pd.DataFrame:
    """
    Filter GO associations by evidence codes and aspects.
    
    Args:
        df: DataFrame with gene-GO associations
        evidence_codes: List of evidence codes to include (None = all)
        aspects: List of aspects to include (None = all)
        
    Returns:
        Filtered DataFrame
    """
    filtered_df = df.copy()
    
    if evidence_codes:
        filtered_df = filtered_df[filtered_df['evidence'].isin(evidence_codes)]
        logger.info(f"Filtered by evidence codes: {len(filtered_df)} associations remaining")
    
    if aspects:
        filtered_df = filtered_df[filtered_df['aspect'].isin(aspects)]
        logger.info(f"Filtered by aspects: {len(filtered_df)} associations remaining")
    
    return filtered_df


def get_go_term_info(go_dag: obo_parser.GODag, go_id: str) -> Dict[str, str]:
    """
    Get information about a GO term.
    
    Args:
        go_dag: GODag object
        go_id: GO term ID
        
    Returns:
        Dictionary with term information
    """
    if go_id not in go_dag:
        return {'id': go_id, 'name': 'Unknown', 'namespace': 'Unknown'}
    
    term = go_dag[go_id]
    return {
        'id': go_id,
        'name': term.name,
        'namespace': term.namespace,
        'definition': getattr(term, 'defn', ''),
        'is_obsolete': term.is_obsolete
    }


def get_go_ancestors(go_dag: obo_parser.GODag, go_id: str) -> Set[str]:
    """
    Get all ancestor GO terms for a given GO term.
    
    Args:
        go_dag: GODag object
        go_id: GO term ID
        
    Returns:
        Set of ancestor GO term IDs
    """
    if go_id not in go_dag:
        return set()
    
    term = go_dag[go_id]
    ancestors = set()
    
    # Get all parents recursively
    def _get_parents(term_id):
        if term_id in go_dag:
            term = go_dag[term_id]
            for parent in term.parents:
                ancestors.add(parent.id)
                _get_parents(parent.id)
    
    _get_parents(go_id)
    return ancestors


def propagate_go_annotations(df: pd.DataFrame, go_dag: obo_parser.GODag) -> pd.DataFrame:
    """
    Propagate GO annotations to ancestor terms.
    
    Args:
        df: DataFrame with gene-GO associations
        go_dag: GODag object
        
    Returns:
        DataFrame with propagated annotations
    """
    logger.info("Propagating GO annotations to ancestor terms")
    
    rows = []
    
    for _, row in df.iterrows():
        gene_id = row['gene_id']
        go_id = row['go_id']
        
        # Add original annotation
        rows.append(row.to_dict())
        
        # Add ancestor annotations
        ancestors = get_go_ancestors(go_dag, go_id)
        for ancestor_id in ancestors:
            ancestor_row = row.to_dict()
            ancestor_row['go_id'] = ancestor_id
            ancestor_row['evidence'] = 'IEA'  # Inferred from Electronic Annotation
            rows.append(ancestor_row)
    
    propagated_df = pd.DataFrame(rows)
    propagated_df = propagated_df.drop_duplicates(['gene_id', 'go_id'])
    
    logger.info(f"Propagated to {len(propagated_df)} total associations")
    return propagated_df


def create_gene2go_mapping(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Create gene-to-GO mapping dictionary.
    
    Args:
        df: DataFrame with gene-GO associations
        
    Returns:
        Dictionary mapping gene IDs to sets of GO term IDs
    """
    gene2go = {}
    
    for _, row in df.iterrows():
        gene_id = row['gene_id']
        go_id = row['go_id']
        
        if gene_id not in gene2go:
            gene2go[gene_id] = set()
        
        gene2go[gene_id].add(go_id)
    
    return gene2go


def get_go_stats(df: pd.DataFrame) -> Dict[str, int]:
    """
    Get basic statistics about GO annotations.
    
    Args:
        df: DataFrame with gene-GO associations
        
    Returns:
        Dictionary with GO statistics
    """
    if df.empty:
        return {}
    
    stats = {
        'total_associations': len(df),
        'unique_genes': df['gene_id'].nunique(),
        'unique_go_terms': df['go_id'].nunique(),
        'avg_terms_per_gene': len(df) / df['gene_id'].nunique()
    }
    
    # Count by aspect
    if 'aspect' in df.columns:
        aspect_counts = df['aspect'].value_counts().to_dict()
        stats.update({f'aspect_{k}': v for k, v in aspect_counts.items()})
    
    return stats