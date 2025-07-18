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


def _preprocess_gaf_file(gaf_file: str) -> str:
    """
    Preprocess GAF file to handle extra empty columns.
    
    Args:
        gaf_file: Path to original GAF file
        
    Returns:
        Path to preprocessed GAF file (or original if no preprocessing needed)
    """
    import tempfile
    import shutil
    
    # First, check if the file needs preprocessing
    needs_preprocessing = False
    max_fields = 0
    
    with open(gaf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('!') or not line.strip():
                continue
            
            fields = line.split('\t')
            max_fields = max(max_fields, len(fields))
            
            # If we have more than 17 fields, less than 16 fields, or trailing empty fields, preprocess
            if len(fields) > 17 or len(fields) < 16 or (len(fields) > 16 and not fields[-1].strip()):
                needs_preprocessing = True
                break
            
            # Only check first 100 data lines for performance
            if line_num > 100:
                break
    
    if not needs_preprocessing:
        logger.debug(f"GAF file appears well-formatted ({max_fields} max fields), no preprocessing needed")
        return gaf_file
    
    # Create preprocessed file
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.gaf', delete=False)
    logger.info(f"Preprocessing GAF file to handle extra columns (max {max_fields} fields found)")
    
    try:
        with open(gaf_file, 'r') as infile:
            for line in infile:
                if line.startswith('!') or not line.strip():
                    temp_file.write(line)
                    continue
                
                fields = line.rstrip('\n\r').split('\t')
                
                # Keep only the first 17 fields
                if len(fields) > 17:
                    fields = fields[:17]
                
                # Remove trailing empty fields, but ensure we have at least 16 fields for GAF 2.0
                while len(fields) > 16 and not fields[-1].strip():
                    fields.pop()
                
                # Ensure we have exactly 17 fields (GAF 2.0 standard + one extra for safety)
                while len(fields) < 17:
                    fields.append('')
                
                temp_file.write('\t'.join(fields) + '\n')
        
        temp_file.close()
        logger.info(f"Preprocessed GAF file saved to: {temp_file.name}")
        return temp_file.name
        
    except Exception as e:
        temp_file.close()
        os.unlink(temp_file.name)
        logger.warning(f"Failed to preprocess GAF file: {e}")
        return gaf_file


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
    
    # Preprocess the file to handle extra columns
    processed_gaf_file = _preprocess_gaf_file(gaf_file)
    cleanup_temp_file = processed_gaf_file != gaf_file
    
    try:
        # Read GAF file using GOATOOLS
        associations = read_gaf(processed_gaf_file)
        
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
        logger.info(f"Parsed {len(df)} gene-GO associations for {df['gene_id'].nunique()} genes using GOATOOLS")
        return df
        
    except (IndexError, ValueError) as e:
        logger.warning(f"GOATOOLS GAF parsing failed (likely due to extra columns): {e}")
        logger.info("Falling back to manual GAF parsing...")
        return _parse_gaf_manual(gaf_file)
    except Exception as e:
        logger.error(f"Error parsing GAF file with GOATOOLS: {e}")
        logger.info("Falling back to manual GAF parsing...")
        return _parse_gaf_manual(gaf_file)
    finally:
        # Clean up temporary file if created
        if cleanup_temp_file and os.path.exists(processed_gaf_file):
            try:
                os.unlink(processed_gaf_file)
                logger.debug(f"Cleaned up temporary file: {processed_gaf_file}")
            except Exception as e:
                logger.warning(f"Failed to clean up temporary file {processed_gaf_file}: {e}")


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
        for line_num, line in enumerate(f, 1):
            # Only strip newlines, not tabs (to preserve trailing empty fields)
            line = line.rstrip('\n\r')
            
            # Skip comments and empty lines
            if line.startswith('!') or not line:
                continue
            
            fields = line.split('\t')
            
            # GAF format requires at least 15 fields (GAF 1.0) or 16 fields (GAF 2.0)
            if len(fields) < 15:
                logger.warning(f"Skipping line {line_num}: insufficient fields ({len(fields)} < 15)")
                continue
            
            try:
                # Extract required fields (0-indexed)
                gene_id = fields[1].strip()   # DB_Object_ID
                go_id = fields[4].strip()     # GO_ID
                evidence = fields[6].strip()  # Evidence_Code
                aspect = fields[8].strip()    # Aspect
                
                # Validate required fields are not empty
                if not all([gene_id, go_id, evidence, aspect]):
                    logger.warning(f"Skipping line {line_num}: empty required fields")
                    continue
                
                # Validate GO ID format
                if not go_id.startswith('GO:'):
                    logger.warning(f"Skipping line {line_num}: invalid GO ID format: {go_id}")
                    continue
                
                rows.append({
                    'gene_id': gene_id,
                    'go_id': go_id,
                    'evidence': evidence,
                    'aspect': aspect
                })
                
            except IndexError as e:
                logger.warning(f"Skipping line {line_num}: IndexError: {e}")
                continue
            except Exception as e:
                logger.warning(f"Skipping line {line_num}: Error: {e}")
                continue
    
    df = pd.DataFrame(rows)
    logger.info(f"Manually parsed {len(df)} gene-GO associations for {df['gene_id'].nunique() if len(df) > 0 else 0} genes")
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