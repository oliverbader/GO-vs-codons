"""
Heatmap visualization module for adaptive GO analysis results.
"""

import os
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
import logging

logger = logging.getLogger(__name__)

# Set style
plt.style.use('default')


def create_adaptive_heatmap(
    results_df: pd.DataFrame,
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (12, 8),
    format: str = 'svg',
    max_terms: int = 50
) -> None:
    """
    Create heatmap of -log10(adj_p) vs thresholds for adaptive GO analysis.
    
    Args:
        results_df: DataFrame with adaptive analysis results
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
        max_terms: Maximum number of GO terms to display
    """
    logger.info("Creating adaptive GO analysis heatmap")
    
    if results_df.empty:
        logger.warning("No results data provided")
        return
    
    # Check if this is codon-specific analysis
    is_codon_analysis = 'codon' in results_df.columns
    
    if is_codon_analysis:
        logger.info("Creating codon-specific heatmap")
        # For codon analysis, create separate heatmaps for each codon
        _create_codon_specific_heatmaps(results_df, output_path, title, figsize, format, max_terms)
    else:
        logger.info("Creating amino acid-level heatmap")
        # Original amino acid-level analysis
        _create_single_heatmap(results_df, output_path, title, figsize, format, max_terms)


def _create_single_heatmap(
    results_df: pd.DataFrame,
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (12, 8),
    format: str = 'svg',
    max_terms: int = 50
) -> None:
    """Create a single heatmap for amino acid-level analysis."""
    
    # Filter to most significant terms
    if len(results_df['go_id'].unique()) > max_terms:
        top_terms = (results_df.groupby('go_id')['adj_p_value'].min()
                    .sort_values().head(max_terms).index)
        results_df = results_df[results_df['go_id'].isin(top_terms)]
    
    # Pivot data for heatmap
    heatmap_data = results_df.pivot_table(
        index='go_id',
        columns='threshold_pct',
        values='adj_p_value',
        fill_value=1.0
    )
    
    # Convert to -log10 scale
    log_p_data = -np.log10(heatmap_data + 1e-10)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    sns.heatmap(
        log_p_data,
        annot=False,
        cmap='viridis',
        cbar_kws={'label': '-log₁₀(adjusted p-value)'},
        ax=ax,
        xticklabels=True,
        yticklabels=True
    )
    
    # Customize plot
    ax.set_xlabel('Threshold (%)', fontsize=12)
    ax.set_ylabel('GO Terms', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title('Adaptive GO Enrichment Analysis', fontsize=14, fontweight='bold')
    
    # Rotate x-axis labels
    ax.tick_params(axis='x', rotation=45)
    ax.tick_params(axis='y', rotation=0)
    
    # Add significance line
    significance_line = -np.log10(0.05)
    ax.axhline(y=significance_line, color='red', linestyle='--', alpha=0.5)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved adaptive heatmap to {output_path}")


def _create_codon_specific_heatmaps(
    results_df: pd.DataFrame,
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (12, 8),
    format: str = 'svg',
    max_terms: int = 50
) -> None:
    """Create separate heatmaps for each codon."""
    
    # Get unique codons
    codons = sorted(results_df['codon'].unique())
    
    # Create a combined heatmap with codons as subplots
    n_codons = len(codons)
    cols = min(3, n_codons)
    rows = (n_codons + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(figsize[0] * cols, figsize[1] * rows))
    if n_codons == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.reshape(1, -1)
    
    for i, codon in enumerate(codons):
        row = i // cols
        col = i % cols
        ax = axes[row, col] if rows > 1 else axes[col]
        
        # Filter data for this codon
        codon_data = results_df[results_df['codon'] == codon].copy()
        
        if codon_data.empty:
            ax.set_visible(False)
            continue
        
        # Filter to most significant terms for this codon
        if len(codon_data['go_id'].unique()) > max_terms:
            top_terms = (codon_data.groupby('go_id')['adj_p_value'].min()
                        .sort_values().head(max_terms).index)
            codon_data = codon_data[codon_data['go_id'].isin(top_terms)]
        
        # Pivot data for heatmap
        heatmap_data = codon_data.pivot_table(
            index='go_id',
            columns='threshold_pct',
            values='adj_p_value',
            fill_value=1.0
        )
        
        if heatmap_data.empty:
            ax.set_visible(False)
            continue
        
        # Convert to -log10 scale
        log_p_data = -np.log10(heatmap_data + 1e-10)
        
        # Create heatmap
        sns.heatmap(
            log_p_data,
            annot=False,
            cmap='viridis',
            cbar=True,
            ax=ax,
            xticklabels=True,
            yticklabels=True
        )
        
        # Get amino acid for this codon
        aa = codon_data['amino_acid'].iloc[0] if 'amino_acid' in codon_data.columns else 'Unknown'
        
        # Customize subplot
        ax.set_title(f'{codon} ({aa})', fontsize=12, fontweight='bold')
        ax.set_xlabel('Threshold (%)', fontsize=10)
        ax.set_ylabel('GO Terms', fontsize=10)
        ax.tick_params(axis='x', rotation=45, labelsize=8)
        ax.tick_params(axis='y', rotation=0, labelsize=8)
    
    # Hide unused subplots
    for i in range(n_codons, rows * cols):
        if rows > 1:
            axes[i // cols, i % cols].set_visible(False)
        else:
            axes[i].set_visible(False)
    
    # Overall title
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold')
    else:
        fig.suptitle('Codon-Specific GO Enrichment Analysis', fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved codon-specific heatmaps to {output_path}")


def create_codon_usage_heatmap(
    df: pd.DataFrame,
    output_path: str,
    amino_acids: Optional[List[str]] = None,
    title: Optional[str] = None,
    figsize: tuple = (12, 8),
    format: str = 'svg'
) -> None:
    """
    Create heatmap of codon usage patterns across genes.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        amino_acids: Optional list of amino acids to include
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info("Creating codon usage heatmap")
    
    if df.empty:
        logger.warning("No codon usage data provided")
        return
    
    # Filter amino acids if specified
    if amino_acids:
        df = df[df['AA'].isin(amino_acids)]
    
    # Create gene x codon matrix
    matrix = df.pivot_table(
        index='gene_id',
        columns=['AA', 'codon'],
        values='rel_usage',
        fill_value=0
    )
    
    # Sample genes if too many
    if len(matrix) > 100:
        matrix = matrix.sample(n=100, random_state=42)
        logger.info("Sampled 100 genes for visualization")
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    sns.heatmap(
        matrix,
        annot=False,
        cmap='viridis',
        cbar_kws={'label': 'Relative Usage'},
        ax=ax,
        xticklabels=True,
        yticklabels=False
    )
    
    # Customize plot
    ax.set_xlabel('Amino Acid - Codon', fontsize=12)
    ax.set_ylabel('Genes', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title('Codon Usage Patterns', fontsize=14, fontweight='bold')
    
    # Rotate x-axis labels
    ax.tick_params(axis='x', rotation=90)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved codon usage heatmap to {output_path}")


def create_correlation_heatmap(
    df: pd.DataFrame,
    output_path: str,
    method: str = 'pearson',
    title: Optional[str] = None,
    figsize: tuple = (10, 8),
    format: str = 'svg'
) -> None:
    """
    Create correlation heatmap between codon usage patterns.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        method: Correlation method ('pearson', 'spearman')
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating correlation heatmap using {method} correlation")
    
    if df.empty:
        logger.warning("No codon usage data provided")
        return
    
    # Create gene x codon matrix
    matrix = df.pivot_table(
        index='gene_id',
        columns=['AA', 'codon'],
        values='rel_usage',
        fill_value=0
    )
    
    # Compute correlation matrix
    corr_matrix = matrix.corr(method=method)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    sns.heatmap(
        corr_matrix,
        annot=False,
        cmap='RdBu_r',
        center=0,
        cbar_kws={'label': f'{method.capitalize()} Correlation'},
        ax=ax,
        xticklabels=True,
        yticklabels=True
    )
    
    # Customize plot
    ax.set_xlabel('Amino Acid - Codon', fontsize=12)
    ax.set_ylabel('Amino Acid - Codon', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'Codon Usage Correlation ({method.capitalize()})', 
                    fontsize=14, fontweight='bold')
    
    # Rotate labels
    ax.tick_params(axis='x', rotation=90)
    ax.tick_params(axis='y', rotation=0)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved correlation heatmap to {output_path}")


def create_go_term_heatmap(
    results_df: pd.DataFrame,
    go_info: Dict[str, Dict[str, str]],
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (12, 10),
    format: str = 'svg',
    max_terms: int = 30
) -> None:
    """
    Create heatmap of GO term enrichment with term names.
    
    Args:
        results_df: DataFrame with enrichment results
        go_info: Dictionary with GO term information
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
        max_terms: Maximum number of terms to display
    """
    logger.info("Creating GO term enrichment heatmap")
    
    if results_df.empty:
        logger.warning("No results data provided")
        return
    
    # Filter to most significant terms
    top_results = results_df.nsmallest(max_terms, 'adj_p_value')
    
    # Prepare data
    plot_data = []
    for _, row in top_results.iterrows():
        go_id = row['go_id']
        go_name = go_info.get(go_id, {}).get('name', go_id)
        
        # Truncate long names
        if len(go_name) > 40:
            go_name = go_name[:37] + '...'
        
        plot_data.append({
            'go_term': f"{go_id}: {go_name}",
            'log_p': -np.log10(row['adj_p_value'] + 1e-10),
            'n_genes': row.get('n_genes', 0),
            'fold_enrichment': row.get('fold_enrichment', 1)
        })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, 
                                   gridspec_kw={'width_ratios': [3, 1]})
    
    # Main heatmap: -log10(p-value)
    heatmap_data = plot_df.set_index('go_term')[['log_p']].T
    
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt='.1f',
        cmap='viridis',
        cbar_kws={'label': '-log₁₀(adjusted p-value)'},
        ax=ax1,
        xticklabels=True,
        yticklabels=False
    )
    
    ax1.set_xlabel('GO Terms', fontsize=12)
    ax1.tick_params(axis='x', rotation=90)
    
    # Side heatmap: fold enrichment
    if 'fold_enrichment' in plot_df.columns:
        fold_data = plot_df.set_index('go_term')[['fold_enrichment']].T
        
        sns.heatmap(
            fold_data,
            annot=True,
            fmt='.1f',
            cmap='Reds',
            cbar_kws={'label': 'Fold Enrichment'},
            ax=ax2,
            xticklabels=True,
            yticklabels=False
        )
        
        ax2.set_xlabel('GO Terms', fontsize=12)
        ax2.tick_params(axis='x', rotation=90)
    
    # Overall title
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold')
    else:
        fig.suptitle('GO Term Enrichment Analysis', fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved GO term heatmap to {output_path}")


def create_threshold_comparison_heatmap(
    results_df: pd.DataFrame,
    output_path: str,
    metric: str = 'adj_p_value',
    title: Optional[str] = None,
    figsize: tuple = (12, 8),
    format: str = 'svg'
) -> None:
    """
    Create heatmap comparing GO terms across different thresholds.
    
    Args:
        results_df: DataFrame with adaptive analysis results
        output_path: Path to save the figure
        metric: Metric to display ('adj_p_value', 'fold_enrichment', 'n_genes')
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating threshold comparison heatmap for {metric}")
    
    if results_df.empty:
        logger.warning("No results data provided")
        return
    
    # Check if metric exists
    if metric not in results_df.columns:
        logger.error(f"Metric {metric} not found in results")
        return
    
    # Pivot data
    pivot_data = results_df.pivot_table(
        index='go_id',
        columns='threshold_pct',
        values=metric,
        fill_value=np.nan
    )
    
    # Transform data based on metric
    if metric == 'adj_p_value':
        # Convert to -log10 scale
        plot_data = -np.log10(pivot_data + 1e-10)
        cbar_label = '-log₁₀(adjusted p-value)'
        cmap = 'viridis'
    elif metric == 'fold_enrichment':
        plot_data = pivot_data
        cbar_label = 'Fold Enrichment'
        cmap = 'Reds'
    else:
        plot_data = pivot_data
        cbar_label = metric.replace('_', ' ').title()
        cmap = 'Blues'
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    sns.heatmap(
        plot_data,
        annot=False,
        cmap=cmap,
        cbar_kws={'label': cbar_label},
        ax=ax,
        xticklabels=True,
        yticklabels=True
    )
    
    # Customize plot
    ax.set_xlabel('Threshold (%)', fontsize=12)
    ax.set_ylabel('GO Terms', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title(f'Threshold Comparison: {cbar_label}', 
                    fontsize=14, fontweight='bold')
    
    # Rotate labels
    ax.tick_params(axis='x', rotation=45)
    ax.tick_params(axis='y', rotation=0)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved threshold comparison heatmap to {output_path}")


def create_clustered_heatmap(
    df: pd.DataFrame,
    output_path: str,
    cluster_method: str = 'ward',
    distance_metric: str = 'euclidean',
    title: Optional[str] = None,
    figsize: tuple = (12, 8),
    format: str = 'svg'
) -> None:
    """
    Create clustered heatmap of codon usage patterns.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        cluster_method: Clustering method for linkage
        distance_metric: Distance metric for clustering
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating clustered heatmap using {cluster_method} clustering")
    
    if df.empty:
        logger.warning("No codon usage data provided")
        return
    
    # Create gene x codon matrix
    matrix = df.pivot_table(
        index='gene_id',
        columns=['AA', 'codon'],
        values='rel_usage',
        fill_value=0
    )
    
    # Sample if too many genes
    if len(matrix) > 100:
        matrix = matrix.sample(n=100, random_state=42)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create clustered heatmap
    sns.clustermap(
        matrix,
        method=cluster_method,
        metric=distance_metric,
        cmap='viridis',
        cbar_kws={'label': 'Relative Usage'},
        figsize=figsize,
        xticklabels=True,
        yticklabels=False
    )
    
    # Save figure
    plt.savefig(output_path, format=format.lower(), bbox_inches='tight', dpi=300)
    plt.close()
    
    logger.info(f"Saved clustered heatmap to {output_path}")


def _save_figure(fig, output_path: str, format: str) -> None:
    """Save figure in specified format."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    if format.lower() == 'pdf':
        fig.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
    elif format.lower() == 'png':
        fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    else:  # default to SVG
        fig.savefig(output_path, format='svg', bbox_inches='tight')


def create_batch_heatmaps(
    results_df: pd.DataFrame,
    output_dir: str,
    format: str = 'svg'
) -> None:
    """
    Create multiple heatmaps for different aspects of the analysis.
    
    Args:
        results_df: DataFrame with analysis results
        output_dir: Output directory
        format: Output format
    """
    logger.info("Creating batch heatmaps")
    
    # Adaptive analysis heatmap
    adaptive_path = os.path.join(output_dir, f'heatmap_adaptive.{format}')
    create_adaptive_heatmap(results_df, adaptive_path, format=format)
    
    # Threshold comparison heatmaps
    for metric in ['adj_p_value', 'n_genes']:
        if metric in results_df.columns:
            threshold_path = os.path.join(output_dir, f'heatmap_threshold_{metric}.{format}')
            create_threshold_comparison_heatmap(results_df, threshold_path, 
                                              metric=metric, format=format)