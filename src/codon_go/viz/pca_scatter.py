"""
PCA scatter plot visualization module for codon usage analysis.
"""

import os
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import logging

logger = logging.getLogger(__name__)

# Set style
plt.style.use('default')
sns.set_palette("husl")


def create_pca_plot(
    df: pd.DataFrame,
    output_path: str,
    color_by: Optional[str] = None,
    title: Optional[str] = None,
    figsize: tuple = (10, 8),
    format: str = 'svg',
    n_components: int = 2
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create PCA scatter plot of codon usage patterns.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        color_by: Column name to color points by
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
        n_components: Number of PCA components
        
    Returns:
        Tuple of (transformed_data, explained_variance_ratio)
    """
    logger.info("Creating PCA scatter plot")
    
    if df.empty:
        logger.warning("No codon usage data provided")
        return np.array([]), np.array([])
    
    # Create gene x codon matrix
    matrix = df.pivot_table(
        index='gene_id',
        columns=['AA', 'codon'],
        values='rel_usage',
        fill_value=0
    )
    
    if matrix.empty:
        logger.warning("No data after pivoting")
        return np.array([]), np.array([])
    
    # Standardize data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    transformed_data = pca.fit_transform(scaled_data)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create scatter plot
    if color_by and color_by in df.columns:
        # Create color mapping
        color_data = df.set_index('gene_id')[color_by].reindex(matrix.index)
        
        if color_data.dtype == 'object':
            # Categorical coloring
            unique_vals = color_data.unique()
            colors = sns.color_palette("husl", len(unique_vals))
            color_map = dict(zip(unique_vals, colors))
            
            for val in unique_vals:
                mask = color_data == val
                ax.scatter(
                    transformed_data[mask, 0],
                    transformed_data[mask, 1],
                    c=[color_map[val]],
                    label=val,
                    alpha=0.7,
                    s=50
                )
            ax.legend()
        else:
            # Continuous coloring
            scatter = ax.scatter(
                transformed_data[:, 0],
                transformed_data[:, 1],
                c=color_data,
                cmap='viridis',
                alpha=0.7,
                s=50
            )
            plt.colorbar(scatter, ax=ax, label=color_by)
    else:
        # Default coloring
        ax.scatter(
            transformed_data[:, 0],
            transformed_data[:, 1],
            alpha=0.7,
            s=50
        )
    
    # Customize plot
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title('PCA of Codon Usage Patterns', fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved PCA plot to {output_path}")
    logger.info(f"Explained variance: {pca.explained_variance_ratio_}")
    
    return transformed_data, pca.explained_variance_ratio_


def create_codon_pca(
    df: pd.DataFrame,
    output_path: str,
    amino_acids: Optional[List[str]] = None,
    title: Optional[str] = None,
    figsize: tuple = (12, 10),
    format: str = 'svg'
) -> None:
    """
    Create PCA plot with codon loadings overlay.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        amino_acids: Optional list of amino acids to include
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info("Creating PCA plot with codon loadings")
    
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
    
    if matrix.empty:
        logger.warning("No data after pivoting")
        return
    
    # Standardize data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)
    
    # Perform PCA
    pca = PCA(n_components=2)
    transformed_data = pca.fit_transform(scaled_data)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Left plot: Gene scatter
    ax1.scatter(
        transformed_data[:, 0],
        transformed_data[:, 1],
        alpha=0.6,
        s=30
    )
    
    ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=12)
    ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', fontsize=12)
    ax1.set_title('Gene Scatter', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Right plot: Codon loadings
    loadings = pca.components_.T
    
    # Color by amino acid
    aa_colors = {}
    unique_aas = [col[0] for col in matrix.columns]
    colors = sns.color_palette("husl", len(set(unique_aas)))
    
    for i, aa in enumerate(set(unique_aas)):
        aa_colors[aa] = colors[i % len(colors)]
    
    for i, (aa, codon) in enumerate(matrix.columns):
        ax2.arrow(
            0, 0,
            loadings[i, 0],
            loadings[i, 1],
            head_width=0.02,
            head_length=0.02,
            fc=aa_colors[aa],
            ec=aa_colors[aa],
            alpha=0.7
        )
        
        # Add codon labels for significant loadings
        if np.sqrt(loadings[i, 0]**2 + loadings[i, 1]**2) > 0.1:
            ax2.text(
                loadings[i, 0] * 1.1,
                loadings[i, 1] * 1.1,
                codon,
                fontsize=8,
                ha='center',
                va='center'
            )
    
    ax2.set_xlabel('PC1 Loadings', fontsize=12)
    ax2.set_ylabel('PC2 Loadings', fontsize=12)
    ax2.set_title('Codon Loadings', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Create legend for amino acids
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=aa_colors[aa], markersize=8, label=aa)
                      for aa in sorted(set(unique_aas))]
    ax2.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.3, 1))
    
    # Overall title
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold')
    else:
        fig.suptitle('PCA Analysis of Codon Usage', fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved codon PCA plot to {output_path}")


def create_3d_pca_plot(
    df: pd.DataFrame,
    output_path: str,
    color_by: Optional[str] = None,
    title: Optional[str] = None,
    figsize: tuple = (10, 8),
    format: str = 'svg'
) -> None:
    """
    Create 3D PCA scatter plot.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        color_by: Column name to color points by
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info("Creating 3D PCA scatter plot")
    
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
    
    if matrix.empty:
        logger.warning("No data after pivoting")
        return
    
    # Standardize data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)
    
    # Perform PCA
    pca = PCA(n_components=3)
    transformed_data = pca.fit_transform(scaled_data)
    
    # Create 3D figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    # Create scatter plot
    if color_by and color_by in df.columns:
        color_data = df.set_index('gene_id')[color_by].reindex(matrix.index)
        
        if color_data.dtype == 'object':
            # Categorical coloring
            unique_vals = color_data.unique()
            colors = sns.color_palette("husl", len(unique_vals))
            color_map = dict(zip(unique_vals, colors))
            
            for val in unique_vals:
                mask = color_data == val
                ax.scatter(
                    transformed_data[mask, 0],
                    transformed_data[mask, 1],
                    transformed_data[mask, 2],
                    c=[color_map[val]],
                    label=val,
                    alpha=0.7,
                    s=50
                )
            ax.legend()
        else:
            # Continuous coloring
            scatter = ax.scatter(
                transformed_data[:, 0],
                transformed_data[:, 1],
                transformed_data[:, 2],
                c=color_data,
                cmap='viridis',
                alpha=0.7,
                s=50
            )
            plt.colorbar(scatter, ax=ax, label=color_by)
    else:
        # Default coloring
        ax.scatter(
            transformed_data[:, 0],
            transformed_data[:, 1],
            transformed_data[:, 2],
            alpha=0.7,
            s=50
        )
    
    # Customize plot
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=12)
    ax.set_zlabel(f'PC3 ({pca.explained_variance_ratio_[2]:.1%})', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title('3D PCA of Codon Usage', fontsize=14, fontweight='bold')
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved 3D PCA plot to {output_path}")


def create_pca_variance_plot(
    df: pd.DataFrame,
    output_path: str,
    max_components: int = 10,
    title: Optional[str] = None,
    figsize: tuple = (10, 6),
    format: str = 'svg'
) -> None:
    """
    Create PCA variance explained plot.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        max_components: Maximum number of components to show
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
    """
    logger.info("Creating PCA variance explained plot")
    
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
    
    if matrix.empty:
        logger.warning("No data after pivoting")
        return
    
    # Standardize data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)
    
    # Perform PCA
    n_components = min(max_components, scaled_data.shape[1])
    pca = PCA(n_components=n_components)
    pca.fit(scaled_data)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Left plot: Individual variance explained
    components = range(1, n_components + 1)
    ax1.bar(components, pca.explained_variance_ratio_, alpha=0.7)
    ax1.set_xlabel('Principal Component', fontsize=12)
    ax1.set_ylabel('Variance Explained', fontsize=12)
    ax1.set_title('Individual Variance Explained', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Right plot: Cumulative variance explained
    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)
    ax2.plot(components, cumulative_variance, 'bo-', alpha=0.7)
    ax2.axhline(y=0.8, color='r', linestyle='--', alpha=0.5, label='80% threshold')
    ax2.axhline(y=0.95, color='r', linestyle='--', alpha=0.5, label='95% threshold')
    ax2.set_xlabel('Principal Component', fontsize=12)
    ax2.set_ylabel('Cumulative Variance Explained', fontsize=12)
    ax2.set_title('Cumulative Variance Explained', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Overall title
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold')
    else:
        fig.suptitle('PCA Variance Analysis', fontsize=16, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved PCA variance plot to {output_path}")


def create_biplot(
    df: pd.DataFrame,
    output_path: str,
    title: Optional[str] = None,
    figsize: tuple = (12, 10),
    format: str = 'svg',
    max_loadings: int = 20
) -> None:
    """
    Create PCA biplot showing both samples and variable loadings.
    
    Args:
        df: DataFrame with codon usage data
        output_path: Path to save the figure
        title: Optional custom title
        figsize: Figure size tuple
        format: Output format ('svg', 'pdf', 'png')
        max_loadings: Maximum number of loading vectors to show
    """
    logger.info("Creating PCA biplot")
    
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
    
    if matrix.empty:
        logger.warning("No data after pivoting")
        return
    
    # Standardize data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)
    
    # Perform PCA
    pca = PCA(n_components=2)
    transformed_data = pca.fit_transform(scaled_data)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot samples
    ax.scatter(
        transformed_data[:, 0],
        transformed_data[:, 1],
        alpha=0.6,
        s=30,
        color='blue',
        label='Genes'
    )
    
    # Plot loadings
    loadings = pca.components_.T
    
    # Select top loadings by magnitude
    loading_magnitudes = np.sqrt(loadings[:, 0]**2 + loadings[:, 1]**2)
    top_indices = np.argsort(loading_magnitudes)[-max_loadings:]
    
    for i in top_indices:
        aa, codon = matrix.columns[i]
        
        # Scale loadings for visualization
        scale_factor = 3
        
        ax.arrow(
            0, 0,
            loadings[i, 0] * scale_factor,
            loadings[i, 1] * scale_factor,
            head_width=0.1,
            head_length=0.1,
            fc='red',
            ec='red',
            alpha=0.7
        )
        
        ax.text(
            loadings[i, 0] * scale_factor * 1.1,
            loadings[i, 1] * scale_factor * 1.1,
            f'{aa}-{codon}',
            fontsize=8,
            ha='center',
            va='center',
            color='red'
        )
    
    # Customize plot
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)', fontsize=12)
    
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold')
    else:
        ax.set_title('PCA Biplot', fontsize=14, fontweight='bold')
    
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    _save_figure(fig, output_path, format)
    plt.close()
    
    logger.info(f"Saved PCA biplot to {output_path}")


def _save_figure(fig, output_path: str, format: str) -> None:
    """Save figure in specified format."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    if format.lower() == 'pdf':
        fig.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
    elif format.lower() == 'png':
        fig.savefig(output_path, format='png', bbox_inches='tight', dpi=300)
    else:  # default to SVG
        fig.savefig(output_path, format='svg', bbox_inches='tight')


def create_batch_pca_plots(
    df: pd.DataFrame,
    output_dir: str,
    format: str = 'svg'
) -> None:
    """
    Create multiple PCA plots for comprehensive analysis.
    
    Args:
        df: DataFrame with codon usage data
        output_dir: Output directory
        format: Output format
    """
    logger.info("Creating batch PCA plots")
    
    # Standard PCA plot
    pca_path = os.path.join(output_dir, f'pca_scatter.{format}')
    create_pca_plot(df, pca_path, format=format)
    
    # PCA with loadings
    loadings_path = os.path.join(output_dir, f'pca_loadings.{format}')
    create_codon_pca(df, loadings_path, format=format)
    
    # Variance explained plot
    variance_path = os.path.join(output_dir, f'pca_variance.{format}')
    create_pca_variance_plot(df, variance_path, format=format)
    
    # Biplot
    biplot_path = os.path.join(output_dir, f'pca_biplot.{format}')
    create_biplot(df, biplot_path, format=format)