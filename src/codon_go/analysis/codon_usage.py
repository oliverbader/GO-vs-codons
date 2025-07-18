"""
Codon usage analysis module for computing relative usage by amino acid.
"""

import json
import os
from typing import Dict, Optional, Set
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import logging

logger = logging.getLogger(__name__)


def compute_relative_usage_by_aa(records: Dict[str, SeqRecord], 
                                codon_table_file: Optional[str] = None,
                                cug_clade: bool = False) -> pd.DataFrame:
    """
    For each gene, compute relative usage of synonymous codons per AA.
    
    Args:
        records: Dictionary of gene_id → SeqRecord
        codon_table_file: Optional path to custom codon table JSON
        cug_clade: If True, use CUG-clade genetic code (CUG codes for Serine)
        
    Returns:
        DataFrame with columns: gene_id, AA, codon, count, rel_usage
    """
    logger.info(f"Computing codon usage for {len(records)} genes")
    
    # Load codon table
    if codon_table_file and os.path.exists(codon_table_file):
        with open(codon_table_file, 'r') as f:
            custom_table = json.load(f)
        codon_to_aa = custom_table['codon_table']
        logger.info(f"Using custom codon table from {codon_table_file}")
    else:
        # Use standard genetic code
        table = CodonTable.unambiguous_dna_by_id[1]
        codon_to_aa = table.forward_table.copy()
        logger.info("Using standard genetic code")
    
    # Apply CUG-clade modification if requested
    if cug_clade:
        codon_to_aa['CTG'] = 'S'  # CUG codes for Serine instead of Leucine
        logger.info("Applied CUG-clade genetic code modification: CTG → Serine")
    
    rows = []
    
    for gene_id, record in records.items():
        seq = str(record.seq).upper()
        
        # Skip if sequence length is not multiple of 3
        if len(seq) % 3 != 0:
            logger.warning(f"Skipping {gene_id}: sequence length not multiple of 3")
            continue
        
        # Count codons
        codon_counts = {}
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            
            # Skip if codon contains ambiguous bases
            if any(base not in 'ATCG' for base in codon):
                continue
            
            if codon in codon_to_aa:
                aa = codon_to_aa[codon]
                
                # Skip stop codons for relative usage calculation
                if aa == '*':
                    continue
                
                key = (aa, codon)
                codon_counts[key] = codon_counts.get(key, 0) + 1
        
        # Compute relative usage per amino acid
        aa_totals = {}
        for (aa, codon), count in codon_counts.items():
            aa_totals[aa] = aa_totals.get(aa, 0) + count
        
        # Add rows for this gene
        for (aa, codon), count in codon_counts.items():
            if aa_totals[aa] > 0:  # Avoid division by zero
                rel_usage = count / aa_totals[aa]
                rows.append({
                    'gene_id': gene_id,
                    'AA': aa,
                    'codon': codon,
                    'count': count,
                    'rel_usage': rel_usage
                })
    
    df = pd.DataFrame(rows)
    logger.info(f"Computed codon usage for {df['gene_id'].nunique()} genes, "
                f"{df['AA'].nunique()} amino acids, {len(df)} total observations")
    
    return df


def get_synonymous_codons(aa: str, 
                         codon_table_file: Optional[str] = None,
                         cug_clade: bool = False) -> Set[str]:
    """
    Get all synonymous codons for a given amino acid.
    
    Args:
        aa: Single letter amino acid code
        codon_table_file: Optional path to custom codon table JSON
        cug_clade: If True, use CUG-clade genetic code
        
    Returns:
        Set of synonymous codons
    """
    # Load codon table
    if codon_table_file and os.path.exists(codon_table_file):
        with open(codon_table_file, 'r') as f:
            custom_table = json.load(f)
        codon_to_aa = custom_table['codon_table']
    else:
        table = CodonTable.unambiguous_dna_by_id[1]
        codon_to_aa = table.forward_table.copy()
    
    # Apply CUG-clade modification if requested
    if cug_clade:
        codon_to_aa['CTG'] = 'S'
    
    # Find all codons that code for this amino acid
    synonymous = set()
    for codon, amino_acid in codon_to_aa.items():
        if amino_acid == aa:
            synonymous.add(codon)
    
    return synonymous


def get_cug_clade_info() -> Dict[str, str]:
    """
    Get information about the CUG-clade genetic code differences.
    
    Returns:
        Dictionary with genetic code differences
    """
    return {
        'CTG': {
            'standard': 'L',  # Leucine in standard code
            'cug_clade': 'S'  # Serine in CUG-clade
        }
    }


def validate_cug_clade_usage(df: pd.DataFrame, cug_clade: bool = False) -> Dict[str, any]:
    """
    Validate and analyze CUG codon usage patterns.
    
    Args:
        df: DataFrame with codon usage data
        cug_clade: Whether CUG-clade genetic code was used
        
    Returns:
        Dictionary with validation results
    """
    results = {
        'cug_clade_mode': cug_clade,
        'ctg_present': 'CTG' in df['codon'].values,
        'ctg_count': 0,
        'ctg_assigned_to': None,
        'leucine_codons': [],
        'serine_codons': []
    }
    
    if 'CTG' in df['codon'].values:
        ctg_data = df[df['codon'] == 'CTG']
        results['ctg_count'] = ctg_data['count'].sum()
        
        # Check which amino acid CTG is assigned to
        ctg_aa = ctg_data['AA'].iloc[0] if len(ctg_data) > 0 else None
        results['ctg_assigned_to'] = ctg_aa
        
        # Validate assignment matches expectation
        if cug_clade and ctg_aa != 'S':
            logger.warning(f"CUG-clade mode enabled but CTG assigned to {ctg_aa}, expected S")
        elif not cug_clade and ctg_aa != 'L':
            logger.warning(f"Standard mode but CTG assigned to {ctg_aa}, expected L")
    
    # Get codon lists for leucine and serine
    leucine_codons = df[df['AA'] == 'L']['codon'].unique().tolist()
    serine_codons = df[df['AA'] == 'S']['codon'].unique().tolist()
    
    results['leucine_codons'] = sorted(leucine_codons)
    results['serine_codons'] = sorted(serine_codons)
    
    return results


def create_genetic_code_comparison(standard_df: pd.DataFrame, 
                                  cug_clade_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compare codon usage between standard and CUG-clade genetic codes.
    
    Args:
        standard_df: DataFrame with standard genetic code usage
        cug_clade_df: DataFrame with CUG-clade genetic code usage
        
    Returns:
        DataFrame with comparison results
    """
    comparison_rows = []
    
    # Focus on CTG codon specifically
    for df, code_type in [(standard_df, 'standard'), (cug_clade_df, 'cug_clade')]:
        ctg_data = df[df['codon'] == 'CTG']
        
        if len(ctg_data) > 0:
            comparison_rows.append({
                'genetic_code': code_type,
                'codon': 'CTG',
                'amino_acid': ctg_data['AA'].iloc[0],
                'total_usage': ctg_data['count'].sum(),
                'mean_rel_usage': ctg_data['rel_usage'].mean(),
                'genes_using': ctg_data['gene_id'].nunique()
            })
    
    # Compare leucine and serine codon families
    for aa in ['L', 'S']:
        for df, code_type in [(standard_df, 'standard'), (cug_clade_df, 'cug_clade')]:
            aa_data = df[df['AA'] == aa]
            
            if len(aa_data) > 0:
                comparison_rows.append({
                    'genetic_code': code_type,
                    'codon': f'{aa}_family',
                    'amino_acid': aa,
                    'total_usage': aa_data['count'].sum(),
                    'mean_rel_usage': aa_data['rel_usage'].mean(),
                    'genes_using': aa_data['gene_id'].nunique()
                })
    
    return pd.DataFrame(comparison_rows)


def _convert_aa_codes_to_single_letter(wobble_aas: Set[str]) -> Set[str]:
    """
    Convert amino acid codes to single-letter format.
    
    Args:
        wobble_aas: Set of amino acid codes (single or three-letter)
        
    Returns:
        Set of single-letter amino acid codes
    """
    # Mapping from three-letter to single-letter codes
    three_to_one = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    
    converted = set()
    for aa in wobble_aas:
        if len(aa) == 3 and aa in three_to_one:
            # Three-letter code
            converted.add(three_to_one[aa])
        elif len(aa) == 1 and aa in three_to_one.values():
            # Single-letter code
            converted.add(aa)
        else:
            logger.warning(f"Unknown amino acid code: {aa}")
    
    return converted


def get_amino_acid_name(single_letter_code: str) -> str:
    """
    Convert single-letter amino acid code to three-letter name.
    
    Args:
        single_letter_code: Single letter amino acid code (e.g., 'A')
        
    Returns:
        Three-letter amino acid name (e.g., 'Ala')
    """
    # Mapping from single-letter to three-letter codes
    one_to_three = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
    }
    
    return one_to_three.get(single_letter_code, single_letter_code)


def filter_wobble_codons(df: pd.DataFrame, wobble_aas: Set[str]) -> pd.DataFrame:
    """
    Filter codon usage data to include only wobble-modified amino acids.
    
    Args:
        df: DataFrame with codon usage data
        wobble_aas: Set of amino acids with wobble modifications (single or three-letter)
        
    Returns:
        Filtered DataFrame
    """
    if not wobble_aas:
        return df
    
    # Convert wobble_aas to single-letter codes to match the AA column
    single_letter_aas = _convert_aa_codes_to_single_letter(wobble_aas)
    logger.debug(f"Converted wobble AAs to single-letter: {wobble_aas} -> {single_letter_aas}")
    
    filtered_df = df[df['AA'].isin(single_letter_aas)].copy()
    logger.info(f"Filtered to wobble AAs: {len(filtered_df)} observations for "
                f"{filtered_df['AA'].nunique()} amino acids")
    
    return filtered_df


def compute_codon_bias_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute codon bias metrics for each gene.
    
    Args:
        df: DataFrame with codon usage data
        
    Returns:
        DataFrame with bias metrics per gene
    """
    metrics = []
    
    for gene_id in df['gene_id'].unique():
        gene_data = df[df['gene_id'] == gene_id]
        
        # Compute effective number of codons (ENC)
        enc = _compute_enc(gene_data)
        
        # Compute codon adaptation index (CAI) - simplified version
        cai = _compute_cai_simple(gene_data)
        
        # Compute relative synonymous codon usage (RSCU)
        rscu_values = _compute_rscu(gene_data)
        
        metrics.append({
            'gene_id': gene_id,
            'enc': enc,
            'cai': cai,
            'mean_rscu': np.mean(list(rscu_values.values())) if rscu_values else 0,
            'std_rscu': np.std(list(rscu_values.values())) if rscu_values else 0
        })
    
    return pd.DataFrame(metrics)


def _compute_enc(gene_data: pd.DataFrame) -> float:
    """
    Compute Effective Number of Codons (ENC) for a gene.
    
    Args:
        gene_data: DataFrame with codon usage for one gene
        
    Returns:
        ENC value
    """
    aa_groups = gene_data.groupby('AA')
    enc_sum = 0
    
    for aa, group in aa_groups:
        n = len(group)  # Number of synonymous codons
        if n > 1:
            # Compute homozygosity
            total_codons = group['count'].sum()
            if total_codons > 0:
                p_squares = [(count / total_codons) ** 2 for count in group['count']]
                homozygosity = sum(p_squares)
                enc_sum += n * homozygosity
    
    return enc_sum if enc_sum > 0 else 0


def _compute_cai_simple(gene_data: pd.DataFrame) -> float:
    """
    Compute simplified Codon Adaptation Index (CAI) for a gene.
    
    Args:
        gene_data: DataFrame with codon usage for one gene
        
    Returns:
        Simplified CAI value
    """
    # Simplified CAI: geometric mean of relative usage values
    rel_usage_values = gene_data['rel_usage'].values
    rel_usage_values = rel_usage_values[rel_usage_values > 0]  # Remove zeros
    
    if len(rel_usage_values) == 0:
        return 0
    
    # Geometric mean
    log_sum = np.sum(np.log(rel_usage_values))
    return np.exp(log_sum / len(rel_usage_values))


def _compute_rscu(gene_data: pd.DataFrame) -> Dict[str, float]:
    """
    Compute Relative Synonymous Codon Usage (RSCU) for a gene.
    
    Args:
        gene_data: DataFrame with codon usage for one gene
        
    Returns:
        Dictionary mapping codons to RSCU values
    """
    rscu_values = {}
    
    for aa in gene_data['AA'].unique():
        aa_data = gene_data[gene_data['AA'] == aa]
        n_synonymous = len(aa_data)
        
        if n_synonymous > 1:
            for _, row in aa_data.iterrows():
                codon = row['codon']
                rel_usage = row['rel_usage']
                rscu_values[codon] = rel_usage * n_synonymous
    
    return rscu_values


def get_codon_usage_stats(df: pd.DataFrame) -> Dict[str, float]:
    """
    Get summary statistics for codon usage data.
    
    Args:
        df: DataFrame with codon usage data
        
    Returns:
        Dictionary with summary statistics
    """
    if df.empty:
        return {}
    
    stats = {
        'total_observations': len(df),
        'unique_genes': df['gene_id'].nunique(),
        'unique_amino_acids': df['AA'].nunique(),
        'unique_codons': df['codon'].nunique(),
        'mean_rel_usage': df['rel_usage'].mean(),
        'std_rel_usage': df['rel_usage'].std(),
        'min_rel_usage': df['rel_usage'].min(),
        'max_rel_usage': df['rel_usage'].max()
    }
    
    # Add per-AA statistics
    aa_stats = df.groupby('AA')['rel_usage'].agg(['mean', 'std', 'count'])
    for aa in aa_stats.index:
        stats[f'mean_rel_usage_{aa}'] = aa_stats.loc[aa, 'mean']
        stats[f'std_rel_usage_{aa}'] = aa_stats.loc[aa, 'std']
        stats[f'count_{aa}'] = aa_stats.loc[aa, 'count']
    
    return stats


def create_codon_frequency_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a gene x codon frequency matrix.
    
    Args:
        df: DataFrame with codon usage data
        
    Returns:
        DataFrame with genes as rows and codons as columns
    """
    # Pivot to create matrix
    matrix = df.pivot_table(
        index='gene_id',
        columns='codon',
        values='rel_usage',
        fill_value=0
    )
    
    logger.info(f"Created codon frequency matrix: {matrix.shape[0]} genes x {matrix.shape[1]} codons")
    return matrix


def identify_optimal_codons(df: pd.DataFrame, top_percentile: float = 0.1) -> Dict[str, str]:
    """
    Identify optimal codons for each amino acid based on high-expression genes.
    
    Args:
        df: DataFrame with codon usage data
        top_percentile: Percentile of genes to consider as highly expressed
        
    Returns:
        Dictionary mapping amino acids to optimal codons
    """
    # For simplicity, use genes with highest total codon counts as proxy for expression
    gene_totals = df.groupby('gene_id')['count'].sum()
    top_genes = gene_totals.quantile(1 - top_percentile)
    high_expr_genes = gene_totals[gene_totals >= top_genes].index
    
    high_expr_data = df[df['gene_id'].isin(high_expr_genes)]
    
    optimal_codons = {}
    
    for aa in high_expr_data['AA'].unique():
        aa_data = high_expr_data[high_expr_data['AA'] == aa]
        
        # Find codon with highest mean relative usage
        codon_means = aa_data.groupby('codon')['rel_usage'].mean()
        optimal_codon = codon_means.idxmax()
        optimal_codons[aa] = optimal_codon
    
    logger.info(f"Identified optimal codons for {len(optimal_codons)} amino acids")
    return optimal_codons