"""
Boxplot visualization module for codon usage distributions.
"""

import os
from typing import Dict, List, Optional, Set
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import logging

logger = logging.getLogger(__name__)

# Set style
plt.style.use('default')
sns.set_palette("husl")


def create_codon_boxplot(
    df: pd.DataFrame,
    amino_acid: str,
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (10, 6),
    format: str = 'svg'
) -> None:
    """
    Create boxplot showing codon usage distribution for a specific amino acid.
    
    Args:
        df: DataFrame with codon usage data
        amino_acid: Single letter amino acid code
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating codon boxplot for amino acid {amino_acid}")
    
    # Filter data for the specific amino acid
    aa_data = df[df['AA'] == amino_acid].copy()
    
    if aa_data.empty:
        logger.warning(f"No data found for amino acid {amino_acid}")
        return
    
    # Get synonymous codons
    codons = sorted(aa_data['codon'].unique())
    
    if len(codons) < 2:
        logger.warning(f"Only {len(codons)} codon(s) found for {amino_acid}, skipping boxplot")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create boxplot
    box_data = [aa_data[aa_data['codon'] == codon]['rel_usage'].values for codon in codons]
    
    bp = ax.boxplot(
        box_data,
        labels=codons,
        patch_artist=True,
        showfliers=True,
        flierprops=dict(marker='o', markersize=4, alpha=0.6)
    )
    
    # Color boxes
    colors = sns.color_palette("husl", len(codons))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Customize plot
    ax.set_xlabel('Codon', fontsize=12)
    ax.set_ylabel('Relative Usage', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'Codon Usage Distribution for {amino_acid}', 
                    fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    _add_boxplot_stats(ax, aa_data, codons)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved codon boxplot to {output_path}")


def create_go_term_boxplot(
    df_codon: pd.DataFrame,
    df_go: pd.DataFrame,
    go_id: str,
    amino_acid: str,
    codon: str,
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (8, 6),
    format: str = 'svg'
) -> None:
    """
    Create boxplot comparing codon usage between GO term and background genes.
    
    Args:
        df_codon: DataFrame with codon usage data
        df_go: DataFrame with gene-GO mappings
        go_id: GO term ID
        amino_acid: Single letter amino acid code
        codon: Codon sequence
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating GO term boxplot for {go_id}, {amino_acid}, {codon}")
    
    # Get genes associated with GO term
    go_genes = set(df_go[df_go['go_id'] == go_id]['gene_id'].unique())
    
    if len(go_genes) == 0:
        logger.warning(f"No genes found for GO term {go_id}")
        return
    
    # Filter codon usage data
    codon_data = df_codon[
        (df_codon['AA'] == amino_acid) & 
        (df_codon['codon'] == codon)
    ].copy()
    
    if codon_data.empty:
        logger.warning(f"No codon usage data found for {amino_acid} {codon}")
        return
    
    # Split into GO term and background
    go_usage = codon_data[codon_data['gene_id'].isin(go_genes)]['rel_usage']
    bg_usage = codon_data[~codon_data['gene_id'].isin(go_genes)]['rel_usage']
    
    if len(go_usage) == 0 or len(bg_usage) == 0:
        logger.warning(f"Insufficient data for comparison")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create boxplot
    box_data = [go_usage.values, bg_usage.values]
    labels = [f'GO:{go_id}\n(n={len(go_usage)})', f'Background\n(n={len(bg_usage)})']
    
    bp = ax.boxplot(
        box_data,
        labels=labels,
        patch_artist=True,
        showfliers=True,
        flierprops=dict(marker='o', markersize=4, alpha=0.6)
    )
    
    # Color boxes
    colors = ['lightcoral', 'lightblue']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Customize plot
    ax.set_ylabel('Relative Usage', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'Codon Usage: {amino_acid} {codon}', 
                    fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add statistical comparison
    _add_comparison_stats(ax, go_usage, bg_usage)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved GO term boxplot to {output_path}")


def create_multi_codon_boxplot(
    df: pd.DataFrame,
    amino_acids: List[str],
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (15, 10),
    format: str = 'svg'
) -> None:
    """
    Create multi-panel boxplot for multiple amino acids.
    
    Args:
        df: DataFrame with codon usage data
        amino_acids: List of amino acid codes
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating multi-codon boxplot for {len(amino_acids)} amino acids")
    
    # Filter data
    aa_data = df[df['AA'].isin(amino_acids)].copy()
    
    if aa_data.empty:
        logger.warning("No data found for specified amino acids")
        return
    
    # Calculate subplot dimensions
    n_aas = len(amino_acids)
    n_cols = min(3, n_aas)
    n_rows = (n_aas + n_cols - 1) // n_cols
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    
    for i, aa in enumerate(amino_acids):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        # Get data for this amino acid
        aa_subset = aa_data[aa_data['AA'] == aa]
        
        if aa_subset.empty:
            ax.text(0.5, 0.5, f'No data for {aa}', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{aa}', fontsize=12, fontweight='bold')
            continue
        
        # Get codons
        codons = sorted(aa_subset['codon'].unique())
        
        if len(codons) < 2:
            ax.text(0.5, 0.5, f'Only {len(codons)} codon(s)', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{aa}', fontsize=12, fontweight='bold')
            continue
        
        # Create boxplot
        box_data = [aa_subset[aa_subset['codon'] == codon]['rel_usage'].values 
                   for codon in codons]
        
        bp = ax.boxplot(
            box_data,
            labels=codons,
            patch_artist=True,
            showfliers=False
        )
        
        # Color boxes
        colors = sns.color_palette("husl", len(codons))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Customize subplot
        ax.set_title(f'{aa}', fontsize=12, fontweight='bold')
        ax.set_ylabel('Relative Usage', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Rotate x-axis labels if needed
        if len(codons) > 3:
            ax.tick_params(axis='x', rotation=45)
    
    # Hide empty subplots
    for i in range(n_aas, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    # Add overall title
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold')
    else:
        fig.suptitle('Codon Usage Distributions', fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved multi-codon boxplot to {output_path}")


def create_comparative_boxplot(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    amino_acid: str,
    codon: str,
    labels: List[str],
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (8, 6),
    format: str = 'svg'
) -> None:
    """
    Create comparative boxplot between two datasets.
    
    Args:
        df1: First DataFrame with codon usage data
        df2: Second DataFrame with codon usage data
        amino_acid: Single letter amino acid code
        codon: Codon sequence
        labels: Labels for the two datasets
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating comparative boxplot for {amino_acid} {codon}")
    
    # Filter data
    data1 = df1[(df1['AA'] == amino_acid) & (df1['codon'] == codon)]['rel_usage']
    data2 = df2[(df2['AA'] == amino_acid) & (df2['codon'] == codon)]['rel_usage']
    
    if len(data1) == 0 or len(data2) == 0:
        logger.warning("Insufficient data for comparison")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create boxplot
    box_data = [data1.values, data2.values]
    plot_labels = [f'{labels[0]}\n(n={len(data1)})', f'{labels[1]}\n(n={len(data2)})']
    
    bp = ax.boxplot(
        box_data,
        labels=plot_labels,
        patch_artist=True,
        showfliers=True,
        flierprops=dict(marker='o', markersize=4, alpha=0.6)
    )
    
    # Color boxes
    colors = ['lightcoral', 'lightblue']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Customize plot
    ax.set_ylabel('Relative Usage', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'Comparative Usage: {amino_acid} {codon}', 
                    fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add statistical comparison
    _add_comparison_stats(ax, data1, data2)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved comparative boxplot to {output_path}")


def _add_boxplot_stats(ax, data: pd.DataFrame, codons: List[str]) -> None:
    """Add statistical annotations to boxplot."""
    # Add mean values as text
    for i, codon in enumerate(codons):
        codon_data = data[data['codon'] == codon]['rel_usage']
        mean_val = codon_data.mean()
        ax.text(i + 1, mean_val, f'μ={mean_val:.3f}', 
               ha='center', va='bottom', fontsize=9, 
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))


def _add_comparison_stats(ax, data1: pd.Series, data2: pd.Series) -> None:
    """Add statistical comparison to boxplot."""
    from scipy.stats import mannwhitneyu
    
    try:
        stat, p_value = mannwhitneyu(data1, data2, alternative='two-sided')
        
        # Add p-value annotation
        y_max = max(data1.max(), data2.max())
        ax.text(0.5, 0.95, f'p = {p_value:.4f}', 
               transform=ax.transAxes, ha='center', va='top',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
        
        # Add means
        ax.text(0.02, 0.95, f'Mean 1: {data1.mean():.3f}', 
               transform=ax.transAxes, ha='left', va='top', fontsize=9)
        ax.text(0.02, 0.90, f'Mean 2: {data2.mean():.3f}', 
               transform=ax.transAxes, ha='left', va='top', fontsize=9)
        
    except Exception as e:
        logger.warning(f"Could not compute statistics: {e}")


def _save_figure(fig, output_path: str, format: str) -> None:
    """Save figure in specified format."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    if format.lower() == 'pdf':
        fig.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
    elif format.lower() == 'png':
        fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    else:  # default to SVG
        fig.savefig(output_path, format='svg', bbox_inches='tight')


def create_batch_boxplots(
    df: pd.DataFrame,
    output_dir: str,
    amino_acids: Optional[List[str]] = None,
    format: str = 'svg'
) -> None:
    """
    Create boxplots for all amino acids in batch.
    
    Args:
        df: DataFrame with codon usage data
        output_dir: Output directory
        amino_acids: Optional list of amino acids (default: all)
        format: Output format
    """
    if amino_acids is None:
        amino_acids = sorted(df['AA'].unique())
    
    logger.info(f"Creating batch boxplots for {len(amino_acids)} amino acids")
    
    for aa in amino_acids:
        output_path = os.path.join(output_dir, f'boxplot_{aa}.{format}')
        create_codon_boxplot(df, aa, output_path, format=format)