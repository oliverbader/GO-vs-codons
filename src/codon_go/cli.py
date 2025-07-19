"""
Command-line interface for the Codon-GO analysis pipeline.
"""

import os
import sys
import logging
from typing import Dict, List, Optional, Set
import click
import pandas as pd

from .utils.config_loader import load_config, validate_config, expand_paths
from .utils.file_utils import ensure_directory, save_dataframe
from .utils.warnings_config import configure_warnings
from .parsers.genome_parser import load_genome_annotations, validate_cds_sequences
from .parsers.go_parser import load_go_ontology, parse_gaf_file, propagate_go_annotations
from .analysis.codon_usage import (compute_relative_usage_by_aa, filter_wobble_codons,
                                   validate_cug_clade_usage, get_cug_clade_info)
from .analysis.stats import adaptive_go_analysis_by_codon, create_gene2go_dict
from .viz.boxplots import create_batch_boxplots
from .viz.heatmap import create_batch_heatmaps
from .viz.pca_scatter import create_batch_pca_plots

# Configure logging (will be enhanced in main functions)
logger = logging.getLogger(__name__)


@click.command()
@click.option('--config', 
              type=click.Path(exists=True),
              help='Path to YAML config file')
@click.option('--genome-dir', 
              type=click.Path(exists=True),
              help='Genome annotation directory (per-contig files)')
@click.option('--go-obo', 
              type=click.Path(exists=True),
              help='GO ontology (.obo)')
@click.option('--go-gaf', 
              type=click.Path(exists=True),
              help='GO annotations (.gaf)')
@click.option('--adaptive-start', 
              type=click.IntRange(1, 100),
              default=75,
              help='Starting relative-usage threshold (%)')
@click.option('--adaptive-step', 
              type=click.IntRange(1, 100),
              default=10,
              help='Threshold decrement per round (%)')
@click.option('--adaptive-rounds', 
              type=click.IntRange(1, 20),
              default=3,
              help='Number of adaptive iterations')
@click.option('--wobble-only', 
              is_flag=True,
              help='Restrict analysis to wobble-modified AAs')
@click.option('--wobble-list', 
              type=click.Path(exists=True),
              help='JSON/YAML list of AAs for wobble filter')
@click.option('--cug-clade', 
              is_flag=True,
              help='Use CUG-clade genetic code (CTG codes for Serine)')
@click.option('--outdir', 
              type=click.Path(),
              help='Output directory (overrides config)')
@click.option('--species', 
              type=str,
              help='Species code to analyze (from config)')
@click.option('--validate-only', 
              is_flag=True,
              help='Only validate CDS sequences, do not run analysis')
@click.option('--skip-validation', 
              is_flag=True,
              help='Skip CDS sequence validation')
@click.option('--figure-format', 
              type=click.Choice(['svg', 'pdf', 'png']),
              default='svg',
              help='Output format for figures')
@click.option('--verbose', '-v', 
              is_flag=True,
              help='Enable verbose logging')
@click.option('--quiet', '-q', 
              is_flag=True,
              help='Enable quiet mode (errors only)')
def main(config: Optional[str],
         genome_dir: Optional[str],
         go_obo: Optional[str],
         go_gaf: Optional[str],
         adaptive_start: int,
         adaptive_step: int,
         adaptive_rounds: int,
         wobble_only: bool,
         wobble_list: Optional[str],
         cug_clade: bool,
         outdir: Optional[str],
         species: Optional[str],
         validate_only: bool,
         skip_validation: bool,
         figure_format: str,
         verbose: bool,
         quiet: bool) -> None:
    """
    Codon-GO Analysis Pipeline
    
    A modular Python pipeline for analyzing codon usage and GO term enrichment
    in eukaryotic genomes, with support for CUG-clade fungi.
    """
    # Load configuration first to get output directory
    try:
        if config:
            config_data = load_config(config)
        else:
            logger.warning("No configuration file provided, using CLI options only")
            config_data = _create_cli_config(
                genome_dir, go_obo, go_gaf, adaptive_start, 
                adaptive_step, adaptive_rounds, wobble_only, cug_clade, outdir
            )
        
        # Set up logging to file
        output_dir = config_data.get('output_dir', 'results')
        log_file = os.path.join(output_dir, 'codon_go_analysis.log')
        _setup_logging(log_file, verbose)
        
    except Exception as e:
        # If config loading fails, set up basic console logging
        logging.basicConfig(
            level=logging.DEBUG if verbose else logging.ERROR if quiet else logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        logger.error(f"Failed to load configuration: {e}")
        return
    
    # Configure warnings based on verbosity
    configure_warnings(verbose=verbose)
    
    logger.info("Starting Codon-GO Analysis Pipeline")
    
    # Display CUG-clade information if requested
    if cug_clade:
        cug_info = get_cug_clade_info()
        logger.info("CUG-clade genetic code enabled:")
        for codon, info in cug_info.items():
            logger.info(f"  {codon}: {info['standard']} (standard) → {info['cug_clade']} (CUG-clade)")
    
    try:
        # Expand paths for file-based config
        if config:
            config_data = expand_paths(config_data, os.path.dirname(config))
        
        # Override config with command line options
        config_data = _override_config(
            config_data, genome_dir, go_obo, go_gaf, 
            adaptive_start, adaptive_step, adaptive_rounds, 
            wobble_only, cug_clade, outdir
        )
        
        # Load wobble list if provided
        wobble_aas = None
        if wobble_list:
            from .utils.config_loader import load_wobble_list
            wobble_aas = set(load_wobble_list(wobble_list))
        elif config_data.get('wobble_only', False):
            wobble_aas = set(config_data.get('wobble_aas', []))
        
        # Determine species to analyze
        species_list = config_data.get('species', [])
        if species:
            species_list = [s for s in species_list if s.get('code') == species]
            if not species_list:
                raise ValueError(f"Species '{species}' not found in configuration")
        
        # Process each species
        for species_config in species_list:
            logger.info(f"Processing species: {species_config['name']} ({species_config['code']})")
            
            # Check if this species is CUG-clade
            species_cug_clade = species_config.get('cug_clade', False) or cug_clade
            if species_cug_clade:
                logger.info(f"Species {species_config['code']} configured as CUG-clade")
            
            success = process_species(
                species_config=species_config,
                go_obo_path=config_data['go_obo'],
                adaptive_config=config_data['adaptive'],
                wobble_aas=wobble_aas,
                output_dir=config_data['output_dir'],
                validate_only=validate_only,
                skip_validation=skip_validation,
                figure_format=figure_format,
                cug_clade=species_cug_clade
            )
            
            if success:
                logger.info(f"Successfully processed {species_config['code']}")
            else:
                logger.error(f"Failed to process {species_config['code']}")
                sys.exit(1)
        
        logger.info("Codon-GO Analysis Pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def process_species(species_config: Dict,
                   go_obo_path: str,
                   adaptive_config: Dict,
                   wobble_aas: Optional[Set[str]],
                   output_dir: str,
                   validate_only: bool,
                   skip_validation: bool,
                   figure_format: str,
                   cug_clade: bool = False) -> bool:
    """
    Process a single species through the complete pipeline.
    
    Args:
        species_config: Species configuration dictionary
        go_obo_path: Path to GO ontology file
        adaptive_config: Adaptive analysis configuration
        wobble_aas: Set of wobble amino acids
        output_dir: Output directory
        validate_only: Only validate sequences
        skip_validation: Skip sequence validation
        figure_format: Output format for figures
        cug_clade: Use CUG-clade genetic code
        
    Returns:
        True if successful, False otherwise
    """
    species_code = species_config['code']
    species_name = species_config['name']
    
    try:
        # Create output directories
        processed_dir = os.path.join(output_dir, 'processed')
        figures_dir = os.path.join(output_dir, 'figures')
        ensure_directory(processed_dir)
        ensure_directory(figures_dir)
        
        # Step 1: Load genome annotations
        logger.info("Loading genome annotations")
        genome_records = load_genome_annotations(species_config['genome_dir'])
        
        if not genome_records:
            logger.error("No genome records loaded")
            return False
        
        # Step 2: Validate CDS sequences
        if not skip_validation:
            logger.info("Validating CDS sequences")
            genome_records = validate_cds_sequences(genome_records)
        
        if validate_only:
            logger.info("Validation complete, exiting as requested")
            return True
        
        # Step 3: Compute codon usage
        logger.info("Computing codon usage")
        codon_usage_df = compute_relative_usage_by_aa(genome_records, cug_clade=cug_clade)
        
        if codon_usage_df.empty:
            logger.error("No codon usage data computed")
            return False
        
        # Validate CUG-clade usage if enabled
        if cug_clade:
            validation_results = validate_cug_clade_usage(codon_usage_df, cug_clade=True)
            logger.info(f"CUG-clade validation results: {validation_results}")
            
            # Save validation results
            validation_path = os.path.join(processed_dir, f'{species_code}_cug_validation.tsv')
            validation_df = pd.DataFrame([validation_results])
            save_dataframe(validation_df, validation_path)
        
        # Save codon usage data
        codon_usage_path = os.path.join(processed_dir, f'{species_code}_codon_usage.tsv')
        save_dataframe(codon_usage_df, codon_usage_path)
        
        # Step 4: Load GO data
        logger.info("Loading GO ontology and annotations")
        go_dag = load_go_ontology(go_obo_path)
        gene2go_df = parse_gaf_file(species_config['gaf'])
        
        if gene2go_df.empty:
            logger.error("No GO annotations loaded")
            return False
        
        # Propagate GO annotations
        logger.info("Propagating GO annotations")
        gene2go_df = propagate_go_annotations(gene2go_df, go_dag)
        
        # Save gene2go data
        gene2go_path = os.path.join(processed_dir, f'{species_code}_gene2go.tsv')
        save_dataframe(gene2go_df, gene2go_path)
        
        # Step 5: Filter for wobble amino acids if requested
        analysis_df = codon_usage_df.copy()
        if wobble_aas:
            logger.info(f"Filtering for wobble amino acids: {wobble_aas}")
            analysis_df = filter_wobble_codons(analysis_df, wobble_aas)
        
        # Step 6: Perform adaptive GO analysis
        logger.info("Performing adaptive GO analysis")
        adaptive_results, diagnostic_data = adaptive_go_analysis_by_codon(
            df_rel=analysis_df,
            gene2go=gene2go_df,
            start_pct=adaptive_config['start_pct'],
            step_pct=adaptive_config['step_pct'],
            rounds=adaptive_config['rounds'],
            wobble_aas=wobble_aas
        )
        
        # Always save diagnostic data (even if no significant results)
        diagnostic_path = os.path.join(processed_dir, f'{species_code}_diagnostic.tsv')
        save_dataframe(diagnostic_data, diagnostic_path)
        logger.info(f"Saved diagnostic data: {len(diagnostic_data)} codon/threshold combinations")
        
        # Create and save summary table
        summary_data = _create_summary_table(diagnostic_data)
        summary_path = os.path.join(processed_dir, f'{species_code}_summary.tsv')
        save_dataframe(summary_data, summary_path)
        logger.info(f"Saved summary table: {len(summary_data)} rows")
        
        if not adaptive_results.empty:
            # Save adaptive results
            adaptive_path = os.path.join(processed_dir, f'{species_code}_adaptive.tsv')
            save_dataframe(adaptive_results, adaptive_path)
            logger.info(f"Saved {len(adaptive_results)} significant GO enrichment results")
        else:
            logger.warning("No significant results from adaptive analysis")
        
        # Step 7: Create visualizations (always create, even with empty results)
        logger.info("Creating visualizations")
        
        try:
            # Boxplots - show codon usage patterns
            create_batch_boxplots(
                df=analysis_df,
                output_dir=figures_dir,
                species_name=species_name,
                format=figure_format
            )
            logger.info("Created boxplots")
        except Exception as e:
            logger.error(f"Error creating boxplots: {e}")
        
        try:
            # Heatmaps - show results or diagnostic data
            if not adaptive_results.empty:
                create_batch_heatmaps(
                    results_df=adaptive_results,
                    output_dir=figures_dir,
                    format=figure_format,
                    diagnostic_df=diagnostic_data
                )
                logger.info("Created heatmaps with significant results")
            else:
                # Create diagnostic heatmap showing gene counts
                _create_diagnostic_heatmap(
                    diagnostic_data,
                    output_dir=figures_dir,
                    species_code=species_code,
                    format=figure_format
                )
                logger.info("Created diagnostic heatmap showing gene counts")
        except Exception as e:
            logger.error(f"Error creating heatmaps: {e}")
        
        try:
            # PCA plots - show codon usage patterns
            create_batch_pca_plots(
                df=analysis_df,
                output_dir=figures_dir,
                format=figure_format
            )
            logger.info("Created PCA plots")
        except Exception as e:
            logger.error(f"Error creating PCA plots: {e}")
        
        logger.info("Visualization creation completed")
        
        return True
        
    except Exception as e:
        logger.error(f"Error processing species {species_code}: {e}")
        return False


def _setup_logging(log_file: str, verbose: bool = False) -> None:
    """Set up logging to both console and file."""
    
    # Clear any existing handlers
    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    
    # Set log level
    log_level = logging.DEBUG if verbose else logging.INFO
    root_logger.setLevel(log_level)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    
    # File handler
    ensure_directory(os.path.dirname(log_file))
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)
    
    logger.info(f"Logging configured - console and file: {log_file}")


def _create_summary_table(diagnostic_data: pd.DataFrame) -> pd.DataFrame:
    """Create a summary table from diagnostic data."""
    if diagnostic_data.empty:
        return pd.DataFrame()
    
    # Group by codon and calculate summary statistics
    summary = diagnostic_data.groupby('codon').agg({
        'amino_acid': 'first',
        'total_genes_above_threshold': ['min', 'max', 'mean'],
        'genes_with_go_annotations': ['min', 'max', 'mean', 'sum'],
        'unique_go_terms': ['min', 'max', 'mean'],
        'passed_min_genes_filter': 'sum'
    }).round(2)
    
    # Flatten column names
    summary.columns = [f'{col[0]}_{col[1]}' if col[1] else col[0] for col in summary.columns]
    summary = summary.reset_index()
    
    # Add percentage of rounds that passed filter
    total_rounds = diagnostic_data['round'].nunique()
    summary['pct_rounds_passed'] = (summary['passed_min_genes_filter_sum'] / total_rounds * 100).round(1)
    
    return summary


def _create_diagnostic_heatmap(
    diagnostic_data: pd.DataFrame,
    output_dir: str,
    species_code: str,
    format: str = 'svg'
) -> None:
    """Create diagnostic heatmap showing gene counts and GO term availability."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from ..analysis.codon_usage import get_amino_acid_name
    
    if diagnostic_data.empty:
        logger.warning("No diagnostic data to visualize")
        return
    
    # Get all possible codons and thresholds for consistent dimensions
    all_codons = sorted(diagnostic_data['codon'].unique())
    all_thresholds = sorted(diagnostic_data['threshold_pct'].unique())
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 10))
    
    # Heatmap 1: Genes above threshold
    pivot1 = diagnostic_data.pivot_table(
        index='codon',
        columns='threshold_pct',
        values='total_genes_above_threshold',
        fill_value=0
    )
    
    # Ensure consistent dimensions by reindexing
    pivot1 = pivot1.reindex(index=all_codons, columns=all_thresholds, fill_value=0)
    
    # Get maximum value for consistent scaling
    max_genes = diagnostic_data['total_genes_above_threshold'].max()
    
    sns.heatmap(
        pivot1,
        annot=True,
        fmt='d',
        cmap='Blues',
        vmin=0,
        vmax=max_genes,
        cbar_kws={'label': 'Number of genes'},
        ax=ax1
    )
    
    # Add N=xyz labels above columns for gene counts
    for i, threshold in enumerate(all_thresholds):
        total_genes = diagnostic_data[diagnostic_data['threshold_pct'] == threshold]['total_genes_above_threshold'].sum()
        ax1.text(i + 0.5, -0.5, f'N={total_genes}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Add amino acid names to y-axis labels
    y_labels = []
    for codon in all_codons:
        aa_data = diagnostic_data[diagnostic_data['codon'] == codon]
        if not aa_data.empty:
            aa = aa_data['amino_acid'].iloc[0]
            aa_name = get_amino_acid_name(aa)
            y_labels.append(f'{codon}\n({aa_name})')
        else:
            y_labels.append(codon)
    
    ax1.set_yticklabels(y_labels, rotation=0, ha='right')
    ax1.set_title(f'{species_code}: Genes Above Threshold by Codon')
    ax1.set_xlabel('Threshold (%)')
    ax1.set_ylabel('Codon')
    
    # Heatmap 2: Genes with GO annotations
    pivot2 = diagnostic_data.pivot_table(
        index='codon',
        columns='threshold_pct',
        values='genes_with_go_annotations',
        fill_value=0
    )
    
    # Ensure consistent dimensions
    pivot2 = pivot2.reindex(index=all_codons, columns=all_thresholds, fill_value=0)
    
    # Get maximum value for consistent scaling
    max_go_genes = diagnostic_data['genes_with_go_annotations'].max()
    
    sns.heatmap(
        pivot2,
        annot=True,
        fmt='d',
        cmap='Reds',
        vmin=0,
        vmax=max_go_genes,
        cbar_kws={'label': 'Number of genes with GO'},
        ax=ax2
    )
    
    # Add N=xyz labels above columns for GO annotations
    for i, threshold in enumerate(all_thresholds):
        total_go_genes = diagnostic_data[diagnostic_data['threshold_pct'] == threshold]['genes_with_go_annotations'].sum()
        ax2.text(i + 0.5, -0.5, f'N={total_go_genes}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax2.set_yticklabels(y_labels, rotation=0, ha='right')
    ax2.set_title(f'{species_code}: Genes with GO Annotations by Codon')
    ax2.set_xlabel('Threshold (%)')
    ax2.set_ylabel('Codon')
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, f'{species_code}_diagnostic_heatmap.{format}')
    plt.savefig(output_path, format=format, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved diagnostic heatmap: {output_path}")


def _create_cli_config(genome_dir: Optional[str],
                      go_obo: Optional[str],
                      go_gaf: Optional[str],
                      adaptive_start: int,
                      adaptive_step: int,
                      adaptive_rounds: int,
                      wobble_only: bool,
                      cug_clade: bool,
                      outdir: Optional[str]) -> Dict:
    """Create configuration from command line arguments."""
    if not all([genome_dir, go_obo, go_gaf]):
        raise ValueError("Must provide either --config or all of --genome-dir, --go-obo, --go-gaf")
    
    config = {
        'species': [{
            'code': 'cli_species',
            'name': 'CLI Species',
            'genome_dir': genome_dir,
            'gaf': go_gaf,
            'cug_clade': cug_clade
        }],
        'go_obo': go_obo,
        'adaptive': {
            'start_pct': adaptive_start,
            'step_pct': adaptive_step,
            'rounds': adaptive_rounds
        },
        'wobble_only': wobble_only,
        'output_dir': outdir or 'results'
    }
    
    return config


def _override_config(config: Dict,
                    genome_dir: Optional[str],
                    go_obo: Optional[str],
                    go_gaf: Optional[str],
                    adaptive_start: int,
                    adaptive_step: int,
                    adaptive_rounds: int,
                    wobble_only: bool,
                    cug_clade: bool,
                    outdir: Optional[str]) -> Dict:
    """Override configuration with command line arguments."""
    if outdir:
        config['output_dir'] = outdir
    
    if go_obo:
        config['go_obo'] = go_obo
    
    # Override adaptive settings
    config['adaptive']['start_pct'] = adaptive_start
    config['adaptive']['step_pct'] = adaptive_step
    config['adaptive']['rounds'] = adaptive_rounds
    
    if wobble_only:
        config['wobble_only'] = True
    
    # Override species settings if provided
    if genome_dir or go_gaf or cug_clade:
        for species in config.get('species', []):
            if genome_dir:
                species['genome_dir'] = genome_dir
            if go_gaf:
                species['gaf'] = go_gaf
            if cug_clade:
                species['cug_clade'] = cug_clade
    
    return config


@click.command()
@click.option('--output', '-o',
              type=click.Path(),
              default='config/example_species.yaml',
              help='Output path for example configuration')
def create_config(output: str) -> None:
    """Create an example configuration file."""
    from .utils.config_loader import create_example_config
    
    logger.info(f"Creating example configuration at {output}")
    create_example_config(output)
    logger.info("Example configuration created successfully")


@click.command()
@click.argument('config_path', type=click.Path(exists=True))
def validate_config_cmd(config_path: str) -> None:
    """Validate a configuration file."""
    logger.info(f"Validating configuration: {config_path}")
    
    try:
        config = load_config(config_path)
        logger.info("Configuration is valid")
        
        # Check file paths
        from .utils.config_loader import validate_file_paths
        missing_paths = validate_file_paths(config)
        
        if missing_paths:
            logger.warning(f"Missing file paths: {missing_paths}")
        else:
            logger.info("All file paths exist")
            
    except Exception as e:
        logger.error(f"Configuration validation failed: {e}")
        sys.exit(1)


@click.command()
def show_cug_info() -> None:
    """Show information about CUG-clade genetic code differences."""
    cug_info = get_cug_clade_info()
    
    click.echo("CUG-Clade Genetic Code Information")
    click.echo("=" * 40)
    click.echo()
    click.echo("CUG-clade fungi use a non-standard genetic code where:")
    
    for codon, info in cug_info.items():
        click.echo(f"  {codon}: {info['standard']} (standard) → {info['cug_clade']} (CUG-clade)")
    
    click.echo()
    click.echo("Common CUG-clade species include:")
    click.echo("  - Candida albicans")
    click.echo("  - Candida tropicalis")
    click.echo("  - Candida parapsilosis")
    click.echo("  - Candida dubliniensis")
    click.echo("  - Debaryomyces hansenii")
    click.echo("  - Lodderomyces elongisporus")
    click.echo()
    click.echo("Use --cug-clade flag or set cug_clade: true in config for these species.")


@click.group()
def cli():
    """Codon-GO Analysis Pipeline CLI."""
    pass


# Add subcommands
cli.add_command(main, name='run')
cli.add_command(create_config, name='create-config')
cli.add_command(validate_config_cmd, name='validate-config')
cli.add_command(show_cug_info, name='show-cug-info')


if __name__ == '__main__':
    cli()