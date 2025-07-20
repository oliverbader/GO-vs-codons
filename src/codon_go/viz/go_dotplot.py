"""
GO enrichment dot plot visualization module.

Creates dot plots showing GO term enrichment results with gene ratios,
gene counts, and significance levels.
"""

import os
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging

logger = logging.getLogger(__name__)

# Set style
plt.style.use('default')


def create_go_dotplot(
    enrichment_df: pd.DataFrame,
    output_path: str,
    title: str = "GO Enrichment Analysis",
    max_terms: int = 25,
    figsize: Tuple[float, float] = (10, 12),
    format: str = 'svg'
) -> None:
    """
    Create a GO enrichment dot plot.
    
    Args:
        enrichment_df: DataFrame with GO enrichment results
        output_path: Path to save the figure
        title: Plot title
        max_terms: Maximum number of terms to display
        figsize: Figure size (width, height)
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info(f"Creating GO dot plot with {len(enrichment_df)} terms")
    
    if enrichment_df.empty:
        logger.warning("No enrichment data provided")
        return
    
    # Prepare data
    plot_data = enrichment_df.copy()
    
    # Add GO term descriptions if available
    if 'description' not in plot_data.columns:
        plot_data['description'] = plot_data['go_id']
    
    # Calculate -log10(adjusted p-value)
    plot_data['neg_log_padj'] = -np.log10(plot_data['adj_p_value'] + 1e-300)
    
    # Sort by significance and take top terms
    plot_data = plot_data.sort_values('adj_p_value').head(max_terms)
    
    # Reverse order for plotting (most significant at top)
    plot_data = plot_data.iloc[::-1]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create scatter plot
    scatter = ax.scatter(
        plot_data['gene_count'],
        range(len(plot_data)),
        s=plot_data['gene_ratio'] * 1000,  # Size by gene ratio
        c=plot_data['neg_log_padj'],       # Color by significance
        cmap='viridis',
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5
    )
    
    # Customize plot
    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels(plot_data['description'], fontsize=10)
    ax.set_xlabel('Gene Count', fontsize=12)
    ax.set_ylabel('GO Terms', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add colorbar for significance
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-Log₁₀(Adjusted p-value)', rotation=270, labelpad=20)
    
    # Create size legend for gene ratio
    sizes = [0.01, 0.02, 0.03, 0.04]
    size_labels = ['1%', '2%', '3%', '4%']
    legend_elements = [
        plt.scatter([], [], s=s*1000, c='gray', alpha=0.7, edgecolors='black', linewidth=0.5)
        for s in sizes
    ]
    
    legend = ax.legend(
        legend_elements, size_labels,
        title='Gene Ratio',
        loc='lower right',
        frameon=True,
        title_fontsize=10,
        fontsize=9
    )
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.8)
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='x')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved GO dot plot to {output_path}")


def create_go_dotplot_by_category(
    enrichment_results: Dict[str, pd.DataFrame],
    output_dir: str,
    species_name: str = "Species",
    codon: str = "Codon",
    threshold: Optional[float] = None,
    max_terms_per_category: int = 20,
    format: str = 'svg'
) -> None:
    """
    Create separate GO dot plots for each category (BP, MF, CC).
    
    Args:
        enrichment_results: Dictionary with category -> enrichment DataFrame
        output_dir: Output directory
        species_name: Species name for titles
        codon: Codon name for titles
        threshold: Threshold percentage for titles
        max_terms_per_category: Maximum terms per category
        format: Output format
    """
    logger.info(f"Creating category-specific GO dot plots for {codon}")
    
    category_names = {
        'BP': 'Biological Process',
        'MF': 'Molecular Function', 
        'CC': 'Cellular Component'
    }
    
    for category, results_df in enrichment_results.items():
        if results_df.empty:
            logger.info(f"No results for category {category}")
            continue
        
        # Create title
        full_category_name = category_names.get(category, category)
        if threshold is not None:
            title = f"{species_name} - {codon} ({threshold}%)\n{full_category_name}"
        else:
            title = f"{species_name} - {codon}\n{full_category_name}"
        
        # Create output path
        output_path = os.path.join(
            output_dir, 
            f"go_dotplot_{codon}_{category}_{species_name.lower().replace(' ', '_')}.{format}"
        )
        
        # Create plot
        create_go_dotplot(
            enrichment_df=results_df,
            output_path=output_path,
            title=title,
            max_terms=max_terms_per_category,
            format=format
        )


def create_combined_go_dotplot(
    enrichment_results: Dict[str, pd.DataFrame],
    output_path: str,
    title: str = "GO Enrichment Analysis",
    max_terms_per_category: int = 15,
    figsize: Tuple[float, float] = (12, 16),
    format: str = 'svg'
) -> None:
    """
    Create a combined GO dot plot with all categories in one figure.
    
    Args:
        enrichment_results: Dictionary with category -> enrichment DataFrame
        output_path: Path to save the figure
        title: Plot title
        max_terms_per_category: Maximum terms per category
        figsize: Figure size (width, height)
        format: Output format
    """
    logger.info("Creating combined GO dot plot")
    
    # Filter and prepare data for each category
    plot_data_list = []
    category_colors = {'BP': 'circle', 'MF': 'triangle', 'CC': 'square'}
    category_names = {
        'BP': 'Biological Process',
        'MF': 'Molecular Function',
        'CC': 'Cellular Component'
    }
    
    for category, results_df in enrichment_results.items():
        if results_df.empty:
            continue
        
        # Take top terms for this category
        category_data = results_df.copy()
        category_data['neg_log_padj'] = -np.log10(category_data['adj_p_value'] + 1e-300)
        category_data = category_data.sort_values('adj_p_value').head(max_terms_per_category)
        category_data['category'] = category
        category_data['category_name'] = category_names.get(category, category)
        
        plot_data_list.append(category_data)
    
    if not plot_data_list:
        logger.warning("No data to plot")
        return
    
    # Combine all data
    combined_data = pd.concat(plot_data_list, ignore_index=True)
    combined_data = combined_data.sort_values(['category', 'adj_p_value'])
    combined_data = combined_data.iloc[::-1]  # Reverse for plotting
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot each category with different markers
    for category in ['BP', 'MF', 'CC']:
        cat_data = combined_data[combined_data['category'] == category]
        if cat_data.empty:
            continue
        
        # Get positions for this category
        positions = [i for i, cat in enumerate(combined_data['category']) if cat == category]
        
        # Choose marker style
        marker = 'o' if category == 'BP' else ('^' if category == 'MF' else 's')
        
        scatter = ax.scatter(
            cat_data['gene_count'],
            positions,
            s=cat_data['gene_ratio'] * 1000,
            c=cat_data['neg_log_padj'],
            cmap='viridis',
            marker=marker,
            alpha=0.7,
            edgecolors='black',
            linewidth=0.5,
            label=category_names[category]
        )
    
    # Customize plot
    ax.set_yticks(range(len(combined_data)))
    ax.set_yticklabels(combined_data['description'], fontsize=9)
    ax.set_xlabel('Gene Count', fontsize=12)
    ax.set_ylabel('GO Terms', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-Log₁₀(Adjusted p-value)', rotation=270, labelpad=20)
    
    # Add category legend
    ax.legend(title='Category', loc='lower right', frameon=True)
    
    # Add size legend
    sizes = [0.01, 0.02, 0.03]
    size_labels = ['1%', '2%', '3%']
    size_legend_elements = [
        plt.scatter([], [], s=s*1000, c='gray', alpha=0.7, edgecolors='black', linewidth=0.5)
        for s in sizes
    ]
    
    size_legend = ax.legend(
        size_legend_elements, size_labels,
        title='Gene Ratio',
        loc='upper right',
        frameon=True
    )
    ax.add_artist(size_legend)  # Keep both legends
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='x')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved combined GO dot plot to {output_path}")


def _save_figure(fig, output_path: str, format: str) -> None:
    """Save figure in the specified format."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    if format.lower() == 'pdf':
        fig.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
    elif format.lower() == 'png':
        fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    else:  # default to SVG
        fig.savefig(output_path, format='svg', bbox_inches='tight')