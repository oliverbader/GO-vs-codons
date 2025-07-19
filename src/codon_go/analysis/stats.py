"""
Statistical analysis module for adaptive GO enrichment analysis.
"""

from typing import Dict, List, Optional, Set, Tuple
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact, ks_2samp
from statsmodels.stats.multitest import multipletests
import logging

logger = logging.getLogger(__name__)


def adaptive_go_analysis_by_codon(
    df_rel: pd.DataFrame,
    gene2go: pd.DataFrame,
    start_pct: float,
    step_pct: float,
    rounds: int,
    wobble_aas: Optional[Set[str]] = None,
    min_genes: int = 5,
    test_method: str = 'mannwhitneyu'
) -> pd.DataFrame:
    """
    Perform adaptive GO-term enrichment at descending thresholds for individual codons.
    
    Args:
        df_rel: DataFrame with relative codon usage data
        gene2go: DataFrame with gene-GO mappings
        start_pct: Starting relative usage threshold (%)
        step_pct: Threshold decrement per round (%)
        rounds: Number of adaptive iterations
        wobble_aas: Optional set of wobble-modified amino acids
        min_genes: Minimum number of genes required for testing
        test_method: Statistical test method ('mannwhitneyu', 'ks_2samp')
        
    Returns:
        DataFrame with enrichment results by codon
    """
    logger.info(f"Starting adaptive GO analysis by codon with {rounds} rounds")
    logger.info(f"Threshold range: {start_pct}% to {start_pct - (rounds-1)*step_pct}%")
    
    # Filter for wobble amino acids if specified
    if wobble_aas:
        # Import the conversion function from codon_usage module
        from .codon_usage import _convert_aa_codes_to_single_letter
        single_letter_aas = _convert_aa_codes_to_single_letter(wobble_aas)
        df_rel = df_rel[df_rel['AA'].isin(single_letter_aas)].copy()
        logger.info(f"Filtered to wobble AAs: {wobble_aas} -> {single_letter_aas}")
    
    # Create gene-GO mapping dictionary
    gene2go_dict = create_gene2go_dict(gene2go)
    
    # Get unique codons to analyze
    unique_codons = df_rel['codon'].unique()
    logger.info(f"Analyzing {len(unique_codons)} individual codons: {sorted(unique_codons)}")
    
    results = []
    
    # Analyze each codon separately
    for codon in sorted(unique_codons):
        logger.info(f"Analyzing codon: {codon}")
        codon_data = df_rel[df_rel['codon'] == codon].copy()
        
        if len(codon_data) == 0:
            logger.warning(f"No data for codon {codon}, skipping")
            continue
        
        # Get amino acid for this codon
        aa = codon_data['AA'].iloc[0]
        
        for round_num in range(rounds):
            threshold = (start_pct - round_num * step_pct) / 100
            logger.debug(f"  Round {round_num + 1}: threshold = {threshold:.3f}")
            
            # Filter genes above threshold for this specific codon
            high_usage_genes = codon_data[codon_data['rel_usage'] >= threshold]['gene_id'].unique()
            logger.info(f"  Found {len(high_usage_genes)} genes above threshold for {codon}")
            
            if len(high_usage_genes) < min_genes:
                logger.debug(f"  Too few genes ({len(high_usage_genes)}) above threshold for {codon}, skipping")
                continue
            
            # Test each GO term for this codon
            go_results = _test_go_terms_for_codon(
                codon_data, gene2go_dict, high_usage_genes, 
                threshold, test_method, min_genes, codon, aa
            )
            
            # Add round and codon information
            for result in go_results:
                result['round'] = round_num + 1
                result['threshold_pct'] = threshold * 100
                result['codon'] = codon
                result['amino_acid'] = aa
                results.append(result)
    
    if not results:
        logger.warning("No results generated from adaptive codon analysis")
        return pd.DataFrame()
    
    # Convert to DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction
    if len(results_df) > 0:
        results_df['adj_p_value'] = multipletests(
            results_df['p_value'], 
            method='fdr_bh'
        )[1]
    
    # Sort by codon, then by adjusted p-value
    results_df = results_df.sort_values(['codon', 'adj_p_value'])
    
    logger.info(f"Completed adaptive codon analysis: {len(results_df)} total tests across {len(unique_codons)} codons")
    return results_df


def adaptive_go_analysis(
    df_rel: pd.DataFrame,
    gene2go: pd.DataFrame,
    start_pct: float,
    step_pct: float,
    rounds: int,
    wobble_aas: Optional[Set[str]] = None,
    min_genes: int = 5,
    test_method: str = 'mannwhitneyu'
) -> pd.DataFrame:
    """
    Perform adaptive GO-term enrichment at descending thresholds.
    
    Args:
        df_rel: DataFrame with relative codon usage data
        gene2go: DataFrame with gene-GO mappings
        start_pct: Starting relative usage threshold (%)
        step_pct: Threshold decrement per round (%)
        rounds: Number of adaptive iterations
        wobble_aas: Optional set of wobble-modified amino acids
        min_genes: Minimum number of genes required for testing
        test_method: Statistical test method ('mannwhitneyu', 'ks_2samp')
        
    Returns:
        DataFrame with enrichment results
    """
    logger.info(f"Starting adaptive GO analysis with {rounds} rounds")
    logger.info(f"Threshold range: {start_pct}% to {start_pct - (rounds-1)*step_pct}%")
    
    # Filter for wobble amino acids if specified
    if wobble_aas:
        # Import the conversion function from codon_usage module
        from .codon_usage import _convert_aa_codes_to_single_letter
        single_letter_aas = _convert_aa_codes_to_single_letter(wobble_aas)
        df_rel = df_rel[df_rel['AA'].isin(single_letter_aas)].copy()
        logger.info(f"Filtered to wobble AAs: {wobble_aas} -> {single_letter_aas}")
    
    # Create gene-GO mapping dictionary
    gene2go_dict = create_gene2go_dict(gene2go)
    
    results = []
    
    for round_num in range(rounds):
        threshold = (start_pct - round_num * step_pct) / 100
        logger.info(f"Round {round_num + 1}: threshold = {threshold:.3f}")
        
        # Filter genes above threshold
        high_usage_genes = df_rel[df_rel['rel_usage'] >= threshold]['gene_id'].unique()
        logger.info(f"Found {len(high_usage_genes)} genes above threshold")
        
        if len(high_usage_genes) < min_genes:
            logger.warning(f"Too few genes ({len(high_usage_genes)}) above threshold, skipping")
            continue
        
        # Test each GO term
        go_results = _test_go_terms(
            df_rel, gene2go_dict, high_usage_genes, 
            threshold, test_method, min_genes
        )
        
        # Add round information
        for result in go_results:
            result['round'] = round_num + 1
            result['threshold_pct'] = threshold * 100
            results.append(result)
    
    if not results:
        logger.warning("No results generated from adaptive analysis")
        return pd.DataFrame()
    
    # Convert to DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction
    if len(results_df) > 0:
        results_df['adj_p_value'] = multipletests(
            results_df['p_value'], 
            method='fdr_bh'
        )[1]
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('adj_p_value')
    
    logger.info(f"Completed adaptive analysis: {len(results_df)} total tests")
    return results_df


def _test_go_terms_for_codon(
    codon_data: pd.DataFrame,
    gene2go_dict: Dict[str, Set[str]],
    high_usage_genes: List[str],
    threshold: float,
    test_method: str,
    min_genes: int,
    codon: str,
    aa: str
) -> List[Dict]:
    """
    Test GO terms for enrichment in high-usage genes for a specific codon.
    
    Args:
        codon_data: DataFrame with relative codon usage data for this codon only
        gene2go_dict: Dictionary mapping genes to GO terms
        high_usage_genes: List of genes with high codon usage
        threshold: Current threshold value
        test_method: Statistical test method
        min_genes: Minimum number of genes required
        codon: The specific codon being analyzed
        aa: The amino acid this codon codes for
        
    Returns:
        List of dictionaries with test results
    """
    results = []
    
    # Get all GO terms associated with high-usage genes
    go_terms = set()
    genes_with_go = 0
    for gene_id in high_usage_genes:
        if gene_id in gene2go_dict:
            go_terms.update(gene2go_dict[gene_id])
            genes_with_go += 1
    
    logger.info(f"  High-usage genes: {len(high_usage_genes)}, with GO annotations: {genes_with_go}")
    logger.info(f"  Testing {len(go_terms)} GO terms for codon {codon}")
    
    if len(go_terms) == 0:
        logger.warning(f"  No GO terms found for high-usage genes of codon {codon}")
        if len(high_usage_genes) > 0:
            sample_genes = list(high_usage_genes)[:5]
            logger.debug(f"  Sample high-usage genes: {sample_genes}")
            logger.debug(f"  Total genes in gene2go_dict: {len(gene2go_dict)}")
            # Check if any of the sample genes are in the dictionary
            for gene in sample_genes:
                if gene in gene2go_dict:
                    logger.debug(f"    {gene}: {len(gene2go_dict[gene])} GO terms")
                else:
                    logger.debug(f"    {gene}: NOT FOUND in gene2go_dict")
    
    for go_id in go_terms:
        # Get genes associated with this GO term
        go_genes = {gene for gene, terms in gene2go_dict.items() if go_id in terms}
        
        # Intersection with high-usage genes
        overlap_genes = set(high_usage_genes) & go_genes
        
        if len(overlap_genes) < min_genes:
            continue
        
        # Get relative usage values for statistical test (only for this codon)
        in_group = codon_data[codon_data['gene_id'].isin(overlap_genes)]['rel_usage']
        out_group = codon_data[~codon_data['gene_id'].isin(overlap_genes)]['rel_usage']
        
        if len(in_group) == 0 or len(out_group) == 0:
            continue
        
        # Perform statistical test
        if test_method == 'mannwhitneyu':
            try:
                stat, p_value = mannwhitneyu(
                    in_group, out_group, 
                    alternative='greater'
                )
            except ValueError:
                continue
        elif test_method == 'ks_2samp':
            try:
                stat, p_value = ks_2samp(in_group, out_group)
            except ValueError:
                continue
        else:
            raise ValueError(f"Unknown test method: {test_method}")
        
        # Calculate effect size
        effect_size = in_group.mean() - out_group.mean()
        
        results.append({
            'go_id': go_id,
            'n_genes': len(overlap_genes),
            'n_total_go_genes': len(go_genes),
            'mean_usage_in': in_group.mean(),
            'mean_usage_out': out_group.mean(),
            'effect_size': effect_size,
            'test_statistic': stat,
            'p_value': p_value,
            'threshold': threshold
        })
    
    return results


def _test_go_terms(
    df_rel: pd.DataFrame,
    gene2go_dict: Dict[str, Set[str]],
    high_usage_genes: List[str],
    threshold: float,
    test_method: str,
    min_genes: int
) -> List[Dict]:
    """
    Test GO terms for enrichment in high-usage genes.
    
    Args:
        df_rel: DataFrame with relative codon usage data
        gene2go_dict: Dictionary mapping genes to GO terms
        high_usage_genes: List of genes with high codon usage
        threshold: Current threshold value
        test_method: Statistical test method
        min_genes: Minimum number of genes required
        
    Returns:
        List of dictionaries with test results
    """
    results = []
    
    # Get all GO terms associated with high-usage genes
    go_terms = set()
    for gene_id in high_usage_genes:
        if gene_id in gene2go_dict:
            go_terms.update(gene2go_dict[gene_id])
    
    logger.info(f"Testing {len(go_terms)} GO terms")
    
    for go_id in go_terms:
        # Get genes associated with this GO term
        go_genes = {gene for gene, terms in gene2go_dict.items() if go_id in terms}
        
        # Intersection with high-usage genes
        overlap_genes = set(high_usage_genes) & go_genes
        
        if len(overlap_genes) < min_genes:
            continue
        
        # Get relative usage values for statistical test
        in_group = df_rel[df_rel['gene_id'].isin(overlap_genes)]['rel_usage']
        out_group = df_rel[~df_rel['gene_id'].isin(overlap_genes)]['rel_usage']
        
        if len(in_group) == 0 or len(out_group) == 0:
            continue
        
        # Perform statistical test
        if test_method == 'mannwhitneyu':
            try:
                stat, p_value = mannwhitneyu(
                    in_group, out_group, 
                    alternative='greater'
                )
            except ValueError:
                continue
        elif test_method == 'ks_2samp':
            try:
                stat, p_value = ks_2samp(in_group, out_group)
            except ValueError:
                continue
        else:
            raise ValueError(f"Unknown test method: {test_method}")
        
        # Calculate effect size
        effect_size = in_group.mean() - out_group.mean()
        
        results.append({
            'go_id': go_id,
            'n_genes': len(overlap_genes),
            'n_total_go_genes': len(go_genes),
            'mean_usage_in': in_group.mean(),
            'mean_usage_out': out_group.mean(),
            'effect_size': effect_size,
            'test_statistic': stat,
            'p_value': p_value,
            'threshold': threshold
        })
    
    return results


def perform_enrichment_test(
    target_genes: Set[str],
    background_genes: Set[str],
    gene2go_dict: Dict[str, Set[str]],
    min_genes: int = 5
) -> pd.DataFrame:
    """
    Perform Fisher's exact test for GO term enrichment.
    
    Args:
        target_genes: Set of target genes
        background_genes: Set of background genes
        gene2go_dict: Dictionary mapping genes to GO terms
        min_genes: Minimum number of genes required
        
    Returns:
        DataFrame with enrichment results
    """
    logger.info(f"Testing enrichment for {len(target_genes)} target genes")
    
    # Get all GO terms
    all_go_terms = set()
    for terms in gene2go_dict.values():
        all_go_terms.update(terms)
    
    results = []
    
    for go_id in all_go_terms:
        # Get genes associated with this GO term
        go_genes = {gene for gene, terms in gene2go_dict.items() if go_id in terms}
        
        # Create contingency table
        in_target_in_go = len(target_genes & go_genes)
        in_target_not_go = len(target_genes - go_genes)
        not_target_in_go = len(background_genes & go_genes) - in_target_in_go
        not_target_not_go = len(background_genes - go_genes) - in_target_not_go
        
        if in_target_in_go < min_genes:
            continue
        
        # Fisher's exact test
        contingency_table = [
            [in_target_in_go, in_target_not_go],
            [not_target_in_go, not_target_not_go]
        ]
        
        try:
            odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
        except ValueError:
            continue
        
        # Calculate enrichment fold
        expected = len(target_genes) * len(go_genes) / len(background_genes)
        fold_enrichment = in_target_in_go / expected if expected > 0 else 0
        
        results.append({
            'go_id': go_id,
            'n_target_genes': in_target_in_go,
            'n_total_go_genes': len(go_genes),
            'n_target_total': len(target_genes),
            'n_background_total': len(background_genes),
            'odds_ratio': odds_ratio,
            'fold_enrichment': fold_enrichment,
            'p_value': p_value
        })
    
    if not results:
        return pd.DataFrame()
    
    # Convert to DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df['adj_p_value'] = multipletests(
        results_df['p_value'], 
        method='fdr_bh'
    )[1]
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('adj_p_value')
    
    logger.info(f"Completed enrichment test: {len(results_df)} GO terms tested")
    return results_df


def create_gene2go_dict(gene2go_df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Create gene-to-GO mapping dictionary from DataFrame.
    
    Args:
        gene2go_df: DataFrame with gene-GO mappings
        
    Returns:
        Dictionary mapping gene IDs to sets of GO term IDs
    """
    gene2go_dict = {}
    
    for _, row in gene2go_df.iterrows():
        gene_id = row['gene_id']
        go_id = row['go_id']
        
        if gene_id not in gene2go_dict:
            gene2go_dict[gene_id] = set()
        
        gene2go_dict[gene_id].add(go_id)
    
    return gene2go_dict


def filter_significant_results(
    results_df: pd.DataFrame,
    p_threshold: float = 0.05,
    min_genes: int = 5,
    min_fold_enrichment: float = 1.5
) -> pd.DataFrame:
    """
    Filter results for significant GO terms.
    
    Args:
        results_df: DataFrame with enrichment results
        p_threshold: Adjusted p-value threshold
        min_genes: Minimum number of genes
        min_fold_enrichment: Minimum fold enrichment
        
    Returns:
        Filtered DataFrame
    """
    if results_df.empty:
        return results_df
    
    filtered_df = results_df[
        (results_df['adj_p_value'] < p_threshold) &
        (results_df['n_genes'] >= min_genes)
    ].copy()
    
    if 'fold_enrichment' in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df['fold_enrichment'] >= min_fold_enrichment
        ]
    
    logger.info(f"Filtered to {len(filtered_df)} significant results")
    return filtered_df


def compute_enrichment_scores(
    results_df: pd.DataFrame,
    score_method: str = 'combined'
) -> pd.DataFrame:
    """
    Compute enrichment scores for ranking GO terms.
    
    Args:
        results_df: DataFrame with enrichment results
        score_method: Method for computing scores ('combined', 'p_value', 'fold')
        
    Returns:
        DataFrame with enrichment scores
    """
    if results_df.empty:
        return results_df
    
    results_df = results_df.copy()
    
    if score_method == 'combined':
        # Combined score: -log10(p) * fold_enrichment
        results_df['enrichment_score'] = (
            -np.log10(results_df['adj_p_value'] + 1e-10) *
            results_df.get('fold_enrichment', 1)
        )
    elif score_method == 'p_value':
        # Score based on p-value only
        results_df['enrichment_score'] = -np.log10(results_df['adj_p_value'] + 1e-10)
    elif score_method == 'fold':
        # Score based on fold enrichment only
        results_df['enrichment_score'] = results_df.get('fold_enrichment', 1)
    else:
        raise ValueError(f"Unknown score method: {score_method}")
    
    # Sort by enrichment score
    results_df = results_df.sort_values('enrichment_score', ascending=False)
    
    return results_df


def compare_thresholds(
    results_df: pd.DataFrame,
    threshold_col: str = 'threshold_pct'
) -> pd.DataFrame:
    """
    Compare GO term enrichment across different thresholds.
    
    Args:
        results_df: DataFrame with adaptive analysis results
        threshold_col: Column name for threshold values
        
    Returns:
        DataFrame with threshold comparison
    """
    if results_df.empty or threshold_col not in results_df.columns:
        return pd.DataFrame()
    
    # Pivot to get GO terms vs thresholds
    pivot_df = results_df.pivot_table(
        index='go_id',
        columns=threshold_col,
        values='adj_p_value',
        fill_value=1.0
    )
    
    # Calculate stability metrics
    stability_metrics = []
    
    for go_id in pivot_df.index:
        p_values = pivot_df.loc[go_id].values
        
        # Count significant thresholds
        sig_count = np.sum(p_values < 0.05)
        
        # Calculate coefficient of variation
        log_p_values = -np.log10(p_values + 1e-10)
        cv = np.std(log_p_values) / np.mean(log_p_values) if np.mean(log_p_values) > 0 else np.inf
        
        stability_metrics.append({
            'go_id': go_id,
            'n_significant_thresholds': sig_count,
            'min_p_value': np.min(p_values),
            'max_p_value': np.max(p_values),
            'mean_p_value': np.mean(p_values),
            'cv_log_p': cv,
            'stability_score': sig_count / len(p_values)
        })
    
    stability_df = pd.DataFrame(stability_metrics)
    stability_df = stability_df.sort_values('stability_score', ascending=False)
    
    return stability_df


def get_enrichment_stats(results_df: pd.DataFrame) -> Dict[str, float]:
    """
    Get summary statistics for enrichment results.
    
    Args:
        results_df: DataFrame with enrichment results
        
    Returns:
        Dictionary with summary statistics
    """
    if results_df.empty:
        return {}
    
    stats = {
        'total_tests': len(results_df),
        'significant_tests': np.sum(results_df['adj_p_value'] < 0.05),
        'mean_p_value': results_df['p_value'].mean(),
        'median_p_value': results_df['p_value'].median(),
        'min_p_value': results_df['p_value'].min(),
        'max_p_value': results_df['p_value'].max()
    }
    
    if 'fold_enrichment' in results_df.columns:
        stats.update({
            'mean_fold_enrichment': results_df['fold_enrichment'].mean(),
            'median_fold_enrichment': results_df['fold_enrichment'].median(),
            'max_fold_enrichment': results_df['fold_enrichment'].max()
        })
    
    if 'n_genes' in results_df.columns:
        stats.update({
            'mean_genes_per_term': results_df['n_genes'].mean(),
            'median_genes_per_term': results_df['n_genes'].median(),
            'max_genes_per_term': results_df['n_genes'].max()
        })
    
    return stats