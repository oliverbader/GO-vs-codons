"""
Classical GO enrichment analysis module.

Implements traditional GO enrichment using Fisher's exact test,
complementing the codon usage distribution analysis.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import logging

logger = logging.getLogger(__name__)


def get_go_term_descriptions(go_terms: Set[str], go_obo_path: str) -> Dict[str, str]:
    """
    Get GO term descriptions from the ontology.
    
    Args:
        go_terms: Set of GO term IDs
        go_obo_path: Path to GO ontology file
        
    Returns:
        Dictionary mapping GO term ID to description
    """
    descriptions = {}
    
    try:
        from goatools.obo_parser import GODag
        
        logger.info(f"Loading GO ontology for descriptions from {go_obo_path}")
        go_dag = GODag(go_obo_path)
        
        for go_term in go_terms:
            if go_term in go_dag:
                descriptions[go_term] = go_dag[go_term].name
            else:
                descriptions[go_term] = go_term  # Fallback to ID
                
    except ImportError:
        logger.warning("GOATOOLS not available, using GO term IDs as descriptions")
        descriptions = {go_term: go_term for go_term in go_terms}
    
    return descriptions


def classify_go_terms_by_category(gene2go_dict: Dict[str, Set[str]], 
                                 go_obo_path: str) -> Dict[str, Dict[str, Set[str]]]:
    """
    Classify GO terms by category (BP, MF, CC).
    
    Args:
        gene2go_dict: Dictionary mapping gene IDs to GO term sets
        go_obo_path: Path to GO ontology file
        
    Returns:
        Dictionary with structure: {category: {gene_id: {go_terms}}}
    """
    try:
        from goatools.obo_parser import GODag
        
        # Load GO ontology
        logger.info(f"Loading GO ontology from {go_obo_path}")
        go_dag = GODag(go_obo_path)
        
        # Initialize category dictionaries
        categorized = {
            'BP': {},  # Biological Process
            'MF': {},  # Molecular Function  
            'CC': {}   # Cellular Component
        }
        
        # Classify each gene's GO terms
        for gene_id, go_terms in gene2go_dict.items():
            categorized['BP'][gene_id] = set()
            categorized['MF'][gene_id] = set()
            categorized['CC'][gene_id] = set()
            
            for go_term in go_terms:
                if go_term in go_dag:
                    namespace = go_dag[go_term].namespace
                    if namespace == 'biological_process':
                        categorized['BP'][gene_id].add(go_term)
                    elif namespace == 'molecular_function':
                        categorized['MF'][gene_id].add(go_term)
                    elif namespace == 'cellular_component':
                        categorized['CC'][gene_id].add(go_term)
        
        # Log statistics
        for category in ['BP', 'MF', 'CC']:
            total_terms = sum(len(terms) for terms in categorized[category].values())
            genes_with_terms = sum(1 for terms in categorized[category].values() if terms)
            logger.info(f"Category {category}: {total_terms} total annotations, {genes_with_terms} genes annotated")
        
        return categorized
        
    except ImportError:
        logger.warning("GOATOOLS not available, using simple GO term classification")
        return _classify_go_terms_simple(gene2go_dict)


def _classify_go_terms_simple(gene2go_dict: Dict[str, Set[str]]) -> Dict[str, Dict[str, Set[str]]]:
    """
    Simple GO term classification based on GO term prefixes.
    Fallback when GOATOOLS is not available.
    """
    categorized = {'BP': {}, 'MF': {}, 'CC': {}}
    
    for gene_id, go_terms in gene2go_dict.items():
        categorized['BP'][gene_id] = set()
        categorized['MF'][gene_id] = set()
        categorized['CC'][gene_id] = set()
        
        # Simple heuristic: assume all terms are BP for now
        # In practice, would need GO ontology or term mapping
        categorized['BP'][gene_id] = go_terms.copy()
    
    return categorized


def classical_go_enrichment(
    target_genes: Set[str],
    background_genes: Set[str],
    gene2go_dict: Dict[str, Set[str]],
    go_obo_path: Optional[str] = None,
    min_genes: int = 5,
    max_genes: int = 500
) -> pd.DataFrame:
    """
    Perform classical GO enrichment analysis using Fisher's exact test.
    
    Args:
        target_genes: Set of genes of interest (e.g., high codon usage)
        background_genes: Set of all genes in the analysis
        gene2go_dict: Dictionary mapping gene IDs to GO term sets
        min_genes: Minimum genes required for a GO term to be tested
        max_genes: Maximum genes allowed for a GO term to be tested
        
    Returns:
        DataFrame with enrichment results
    """
    logger.info(f"Classical GO enrichment: {len(target_genes)} target genes, {len(background_genes)} background genes")
    
    # Get all GO terms and their gene counts
    go_term_genes = {}
    for gene_id, go_terms in gene2go_dict.items():
        if gene_id in background_genes:
            for go_term in go_terms:
                if go_term not in go_term_genes:
                    go_term_genes[go_term] = set()
                go_term_genes[go_term].add(gene_id)
    
    # Filter GO terms by gene count
    filtered_terms = {
        go_term: genes for go_term, genes in go_term_genes.items()
        if min_genes <= len(genes) <= max_genes
    }
    
    logger.info(f"Testing {len(filtered_terms)} GO terms (after filtering by gene count)")
    
    results = []
    
    for go_term, go_genes in filtered_terms.items():
        # Create 2x2 contingency table
        # |              | Has GO term | No GO term | Total |
        # |--------------|-------------|------------|-------|
        # | Target genes | a           | b          | a+b   |
        # | Other genes  | c           | d          | c+d   |
        # | Total        | a+c         | b+d        | n     |
        
        a = len(target_genes & go_genes)  # Target genes with GO term
        b = len(target_genes - go_genes)  # Target genes without GO term
        c = len(go_genes - target_genes)  # Non-target genes with GO term
        d = len(background_genes - target_genes - go_genes)  # Non-target genes without GO term
        
        # Skip if no overlap
        if a == 0:
            continue
        
        # Fisher's exact test
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        
        # Calculate additional statistics
        total_with_term = a + c
        total_target = a + b
        expected = (total_target * total_with_term) / len(background_genes)
        fold_enrichment = (a / total_target) / (total_with_term / len(background_genes)) if total_with_term > 0 else 0
        gene_ratio = a / total_target if total_target > 0 else 0
        bg_ratio = total_with_term / len(background_genes) if len(background_genes) > 0 else 0
        
        results.append({
            'go_id': go_term,
            'target_genes_with_term': a,
            'total_target_genes': total_target,
            'bg_genes_with_term': total_with_term,
            'total_bg_genes': len(background_genes),
            'expected': expected,
            'fold_enrichment': fold_enrichment,
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'gene_ratio': gene_ratio,
            'bg_ratio': bg_ratio,
            'gene_count': a,
            'description': f"{a}/{total_target} vs {total_with_term}/{len(background_genes)}"
        })
    
    if not results:
        logger.warning("No enriched GO terms found")
        return pd.DataFrame()
    
    # Convert to DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df['adj_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
    
    # Add GO term descriptions
    if go_obo_path and not results_df.empty:
        go_terms_in_results = set(results_df['go_id'])
        descriptions = get_go_term_descriptions(go_terms_in_results, go_obo_path)
        results_df['description'] = results_df['go_id'].map(descriptions)
    else:
        results_df['description'] = results_df['go_id'] if not results_df.empty else []
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('adj_p_value')
    
    logger.info(f"Found {len(results_df)} GO terms with classical enrichment results")
    
    return results_df


def combined_codon_go_analysis(
    codon_usage_df: pd.DataFrame,
    gene2go_dict: Dict[str, Set[str]],
    codon: str,
    thresholds: List[float],
    go_obo_path: str,
    min_genes: int = 5
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Perform combined codon-GO analysis with both classical and distribution-based methods.
    
    Args:
        codon_usage_df: DataFrame with codon usage data
        gene2go_dict: Dictionary mapping gene IDs to GO term sets
        codon: Codon to analyze
        thresholds: List of relative usage thresholds
        go_obo_path: Path to GO ontology file
        min_genes: Minimum genes required for analysis
        
    Returns:
        Dictionary with structure: {method: {category: results_df}}
    """
    from .stats import adaptive_go_analysis_by_codon
    
    logger.info(f"Combined codon-GO analysis for {codon}")
    
    # Classify GO terms by category
    categorized_go = classify_go_terms_by_category(gene2go_dict, go_obo_path)
    
    results = {
        'classical': {'BP': pd.DataFrame(), 'MF': pd.DataFrame(), 'CC': pd.DataFrame()},
        'distribution': {'BP': pd.DataFrame(), 'MF': pd.DataFrame(), 'CC': pd.DataFrame()}
    }
    
    # Get codon usage data
    codon_data = codon_usage_df[codon_usage_df['codon'] == codon].copy()
    if codon_data.empty:
        logger.warning(f"No data found for codon {codon}")
        return results
    
    background_genes = set(codon_data['gene_id'].unique())
    
    # Analyze each category
    for category in ['BP', 'MF', 'CC']:
        logger.info(f"Analyzing category {category}")
        category_gene2go = categorized_go[category]
        
        # Filter to genes with data
        filtered_gene2go = {
            gene_id: go_terms for gene_id, go_terms in category_gene2go.items()
            if gene_id in background_genes and go_terms
        }
        
        if not filtered_gene2go:
            logger.warning(f"No GO annotations found for category {category}")
            continue
        
        # Classical enrichment for each threshold
        classical_results = []
        for threshold in thresholds:
            # Get high-usage genes at this threshold
            threshold_percentile = np.percentile(codon_data['rel_usage'], threshold)
            high_usage_genes = set(
                codon_data[codon_data['rel_usage'] >= threshold_percentile]['gene_id']
            )
            
            if len(high_usage_genes) < min_genes:
                continue
            
            # Classical enrichment
            enrichment_df = classical_go_enrichment(
                target_genes=high_usage_genes,
                background_genes=background_genes,
                gene2go_dict=filtered_gene2go,
                go_obo_path=go_obo_path,
                min_genes=min_genes
            )
            
            if not enrichment_df.empty:
                enrichment_df['threshold'] = threshold
                enrichment_df['codon'] = codon
                enrichment_df['category'] = category
                enrichment_df['method'] = 'classical'
                classical_results.append(enrichment_df)
        
        if classical_results:
            results['classical'][category] = pd.concat(classical_results, ignore_index=True)
        
        # Distribution-based analysis (existing method)
        try:
            # Create temporary gene2go DataFrame for compatibility
            gene2go_df = []
            for gene_id, go_terms in filtered_gene2go.items():
                for go_term in go_terms:
                    gene2go_df.append({'gene_id': gene_id, 'go_id': go_term})
            
            if gene2go_df:
                gene2go_df = pd.DataFrame(gene2go_df)
                
                # Run distribution-based analysis
                dist_results, _ = adaptive_go_analysis_by_codon(
                    df_rel=codon_usage_df,
                    gene2go=gene2go_df,
                    start_pct=max(thresholds),
                    step_pct=5,  # Will be calculated from thresholds
                    rounds=len(thresholds),
                    min_genes=min_genes
                )
                
                if not dist_results.empty:
                    # Filter to current codon and category
                    codon_results = dist_results[dist_results['codon'] == codon].copy()
                    codon_results['category'] = category
                    codon_results['method'] = 'distribution'
                    results['distribution'][category] = codon_results
                    
        except Exception as e:
            logger.error(f"Error in distribution analysis for {category}: {e}")
    
    return results