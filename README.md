# Codon-GO Analysis Pipeline

A modular Python pipeline for analyzing codon usage and GO term enrichment in eukaryotic genomes, with specialized support for CUG-clade fungi.

## Overview

The Codon-GO Analysis Pipeline is a comprehensive tool designed to:

1. **Parse** eukaryotic genome annotations (multi-file EMBL/GenBank) and GO data
2. **Compute** per-gene relative codon usage by amino acid
3. **Perform** adaptive GO-term enrichment at descending thresholds, with optional wobble-modification filtering
4. **Visualize** results in publication-quality figures
5. **Support CUG-clade fungi** that use a non-standard genetic code (CTG codes for serine instead of leucine)

## Features

- **Modular Architecture**: Well-organized modules for parsing, analysis, and visualization
- **Multiple Input Formats**: Support for EMBL, GenBank, and GAF files
- **Adaptive Analysis**: Threshold-based GO enrichment analysis
- **Wobble Filtering**: Focus on wobble-position modified amino acids
- **CUG-Clade Support**: Specialized handling for fungi with non-standard genetic codes
- **Rich Visualizations**: Boxplots, heatmaps, and PCA scatter plots
- **Flexible Configuration**: YAML-based configuration with CLI overrides
- **Publication Ready**: SVG/PDF output for high-quality figures

## CUG-Clade Fungi Support

This pipeline includes specialized support for CUG-clade fungi, which use a non-standard genetic code where **CTG codes for serine instead of leucine**. This is crucial for accurate codon usage analysis in species such as:

- **Candida albicans**
- **Candida tropicalis**
- **Candida parapsilosis**
- **Candida dubliniensis**
- **Debaryomyces hansenii**
- **Lodderomyces elongisporus**

### Using CUG-Clade Mode

#### Via Configuration File
```yaml
species:
  - code: candida_albicans
    name: Candida albicans
    genome_dir: data/genomes/candida_albicans/
    gaf: data/go/mappings/candida_albicans.gaf
    cug_clade: true  # Enable CUG-clade genetic code
```

#### Via Command Line
```bash
python codon_go.py run --config config/species.yaml --cug-clade
```

#### Show CUG-Clade Information
```bash
python codon_go.py show-cug-info
```

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Install Package (Optional)

```bash
pip install -e .
```

## Quick Start

### 1. Create Configuration

```bash
python codon_go.py create-config --output config/my_species.yaml
```

### 2. Edit Configuration

Edit the generated YAML file to specify your species data:

```yaml
species:
  - code: afumigatus
    name: Aspergillus fumigatus
    genome_dir: data/genomes/afumigatus/
    gaf: data/go/mappings/afumigatus.gaf
    cug_clade: false
  - code: candida_albicans
    name: Candida albicans
    genome_dir: data/genomes/candida_albicans/
    gaf: data/go/mappings/candida_albicans.gaf
    cug_clade: true  # CUG-clade species
go_obo: data/go/go.obo
adaptive:
  start_pct: 75
  step_pct: 10
  rounds: 3
wobble_only: true
wobble_aas:
  - Leu
  - Lys
  - Gln
  - Glu
  - Phe
  - Trp
  - Ser  # Important for CUG-clade analysis
output_dir: results
```

### 3. Run Analysis

```bash
python codon_go.py run --config config/my_species.yaml
```

## Usage

### Command Line Interface

The pipeline provides a comprehensive CLI with multiple subcommands:

```bash
# Run complete analysis
python codon_go.py run --config config/species.yaml

# Run with CUG-clade support
python codon_go.py run --config config/species.yaml --cug-clade

# Run with CLI options (override config)
python codon_go.py run \
    --config config/species.yaml \
    --adaptive-start 80 \
    --adaptive-rounds 5 \
    --wobble-only \
    --cug-clade \
    --outdir results_custom

# Run without config file
python codon_go.py run \
    --genome-dir data/genomes/species/ \
    --go-obo data/go/go.obo \
    --go-gaf data/go/mappings/species.gaf \
    --cug-clade \
    --outdir results

# Show CUG-clade information
python codon_go.py show-cug-info

# Validate configuration
python codon_go.py validate-config config/species.yaml

# Create example configuration
python codon_go.py create-config --output config/example.yaml
```

### Available Options

- `--config`: Path to YAML configuration file
- `--genome-dir`: Directory containing EMBL/GenBank files
- `--go-obo`: GO ontology file (.obo)
- `--go-gaf`: GO annotation file (.gaf)
- `--adaptive-start`: Starting threshold percentage (default: 75)
- `--adaptive-step`: Step size for threshold reduction (default: 10)
- `--adaptive-rounds`: Number of analysis rounds (default: 3)
- `--wobble-only`: Restrict analysis to wobble-modified amino acids
- `--wobble-list`: Custom wobble amino acid list (JSON/YAML)
- `--cug-clade`: Use CUG-clade genetic code (CTG codes for serine)
- `--outdir`: Output directory
- `--species`: Specific species code to analyze
- `--validate-only`: Only validate CDS sequences
- `--skip-validation`: Skip CDS validation
- `--figure-format`: Output format (svg/pdf/png)
- `--verbose`: Enable verbose logging
- `--quiet`: Quiet mode (errors only)

## Directory Structure

```
project_root/
├── data/
│   ├── genomes/           # per-species dirs, each with annotation files
│   │   ├── <species_code>/
│   │   │   ├── chr1.embl
│   │   │   ├── chr2.gb
│   │   │   └── ...
│   │   └── ...
│   ├── go/
│   │   ├── go.obo         # GO ontology file
│   │   └── mappings/
│   │       ├── <species_code>.gaf
│   │       └── ...
│   └── references/
│       └── codon_table.json  # optional custom codon→AA mapping
├── config/
│   └── species.yaml       # configuration file
├── results/
│   ├── processed/         # TSV output tables
│   └── figures/           # visualization outputs
└── src/
    └── codon_go/          # main package
        ├── cli.py         # command-line interface
        ├── parsers/       # data parsing modules
        ├── analysis/      # analysis modules
        ├── viz/           # visualization modules
        └── utils/         # utility modules
```

## Output Files

### Processed Data (TSV format)

- `<species>_codon_usage.tsv`: Per-gene codon usage statistics
- `<species>_gene2go.tsv`: Gene-to-GO term mappings
- `<species>_adaptive.tsv`: Adaptive enrichment results
- `<species>_cug_validation.tsv`: CUG-clade validation results (when applicable)

### Figures (SVG/PDF/PNG)

- `heatmap_adaptive.svg`: Adaptive analysis heatmap
- `boxplot_<AA>.svg`: Codon usage boxplots per amino acid
- `pca_scatter.svg`: PCA scatter plot of codon usage
- `pca_loadings.svg`: PCA with codon loadings
- `pca_variance.svg`: PCA variance explained
- `cug_clade_comparison.svg`: CUG-clade specific comparisons (when applicable)

## Configuration

### Species Configuration

```yaml
species:
  - code: species_code          # Short identifier
    name: Species Name          # Full species name
    genome_dir: path/to/genome/ # Directory with annotation files
    gaf: path/to/species.gaf    # GO annotation file
    cug_clade: false           # Use CUG-clade genetic code (default: false)
```

### Analysis Parameters

```yaml
adaptive:
  start_pct: 75     # Starting threshold (75%)
  step_pct: 10      # Decrement per round (10%)
  rounds: 3         # Number of rounds

wobble_only: true   # Restrict to wobble amino acids
wobble_aas:         # List of wobble amino acids
  - Leu
  - Lys
  - Gln
  - Glu
  - Phe
  - Trp
  - Ser  # Important for CUG-clade species
```

## CUG-Clade Analysis Features

### Genetic Code Handling
- Automatic detection and handling of CTG → Serine assignment
- Validation of genetic code consistency
- Comparison tools for standard vs. CUG-clade analysis

### Specialized Visualizations
- CTG codon usage comparisons
- Leucine vs. Serine codon family analysis
- Side-by-side genetic code comparisons

### Validation Tools
- CUG-clade genetic code validation
- Codon assignment verification
- Statistical comparison between genetic codes

## Data Requirements

### Genome Annotations

- **Format**: EMBL (.embl) or GenBank (.gb, .gbk, .genbank)
- **Content**: CDS features with proper gene identifiers
- **Structure**: One file per chromosome/contig in species directory

### GO Data

- **Ontology**: GO ontology file in OBO format
- **Annotations**: Gene Association File (GAF) format
- **Source**: Download from [Gene Ontology Consortium](http://geneontology.org/)

## Analysis Workflow

1. **Genome Parsing**: Extract CDS sequences from annotation files
2. **Sequence Validation**: Validate CDS sequences for proper start/stop codons
3. **Genetic Code Selection**: Apply standard or CUG-clade genetic code
4. **Codon Usage**: Calculate relative usage of synonymous codons
5. **GO Processing**: Load ontology and propagate annotations
6. **Adaptive Analysis**: Perform threshold-based enrichment analysis
7. **Visualization**: Generate publication-quality figures

## Visualization Types

### Boxplots
- Distribution of codon usage for each amino acid
- Comparison between GO terms and background
- CUG-clade specific comparisons (CTG usage, leucine vs. serine families)

### Heatmaps
- Adaptive analysis results across thresholds
- Codon usage patterns across genes
- GO term enrichment with significance levels

### PCA Plots
- Principal component analysis of codon usage
- Gene clustering based on codon preferences
- Codon loading vectors and variance explained

## Dependencies

- **Biopython**: Sequence parsing and analysis
- **GOATOOLS**: GO ontology processing
- **pandas**: Data manipulation and analysis
- **NumPy**: Numerical computations
- **SciPy**: Statistical functions
- **matplotlib**: Basic plotting
- **seaborn**: Statistical visualizations
- **scikit-learn**: Machine learning (PCA)
- **Click**: Command-line interface
- **PyYAML**: Configuration file parsing
- **statsmodels**: Statistical modeling

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Install all required packages from requirements.txt
2. **File Format Errors**: Ensure EMBL/GenBank files are properly formatted
3. **Memory Issues**: Large genomes may require more RAM for analysis
4. **GO File Errors**: Verify GO ontology and GAF files are current
5. **CUG-Clade Issues**: Ensure species are correctly marked as CUG-clade

### Dependency Warnings

The pipeline automatically suppresses common deprecation warnings from dependencies (such as the `pkg_resources` warning from older packages) to keep the output clean. If you need to see all warnings for debugging:

```bash
python codon_go.py run --config config/species.yaml --verbose
```

### Testing Warning Suppression

You can test that warning suppression is working correctly:

```bash
python check_warnings.py
```

### Logging

Enable verbose logging for detailed information:

```bash
python codon_go.py run --config config/species.yaml --verbose
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this pipeline in your research, please cite:

```
Codon-GO Analysis Pipeline: A modular Python tool for analyzing codon usage 
and GO term enrichment in eukaryotic genomes, with support for CUG-clade fungi. (2024)
```

## Contact

For questions, issues, or contributions, please contact:
- Email: codon-go@example.com
- GitHub: https://github.com/example/codon-go

## Acknowledgments

- Gene Ontology Consortium for GO data and tools
- Biopython developers for sequence analysis tools
- Scientific Python community for analysis libraries
- CUG-clade research community for genetic code insights