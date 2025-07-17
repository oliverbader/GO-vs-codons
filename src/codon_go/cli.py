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
from .analysis.stats import adaptive_go_analysis, create_gene2go_dict
from .viz.boxplots import create_batch_boxplots
from .viz.heatmap import create_batch_heatmaps
from .viz.pca_scatter import create_batch_pca_plots

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
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
    # Configure logging level and warnings
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif quiet:
        logging.getLogger().setLevel(logging.ERROR)
    
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
        # Load and validate configuration
        if config:
            config_data = load_config(config)
            config_data = expand_paths(config_data, os.path.dirname(config))
        else:
            # Create minimal config from command line options
            config_data = _create_cli_config(
                genome_dir, go_obo, go_gaf, adaptive_start, 
                adaptive_step, adaptive_rounds, wobble_only, cug_clade, outdir
            )
        
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
        adaptive_results = adaptive_go_analysis(
            df_rel=analysis_df,
            gene2go=gene2go_df,
            start_pct=adaptive_config['start_pct'],
            step_pct=adaptive_config['step_pct'],
            rounds=adaptive_config['rounds'],
            wobble_aas=wobble_aas
        )
        
        if not adaptive_results.empty:
            # Save adaptive results
            adaptive_path = os.path.join(processed_dir, f'{species_code}_adaptive.tsv')
            save_dataframe(adaptive_results, adaptive_path)
            
            # Step 7: Create visualizations
            logger.info("Creating visualizations")
            
            # Boxplots
            create_batch_boxplots(
                df=analysis_df,
                output_dir=figures_dir,
                format=figure_format
            )
            
            # Heatmaps
            create_batch_heatmaps(
                results_df=adaptive_results,
                output_dir=figures_dir,
                format=figure_format
            )
            
            # PCA plots
            create_batch_pca_plots(
                df=analysis_df,
                output_dir=figures_dir,
                format=figure_format
            )
            
            logger.info("Visualizations created successfully")
        else:
            logger.warning("No significant results from adaptive analysis")
        
        return True
        
    except Exception as e:
        logger.error(f"Error processing species {species_code}: {e}")
        return False


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