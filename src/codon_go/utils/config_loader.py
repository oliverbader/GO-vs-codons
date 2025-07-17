"""
Configuration loader utility module.
"""

import os
import yaml
import json
from typing import Dict, Any, List, Optional
import logging

logger = logging.getLogger(__name__)


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Dictionary with configuration data
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    logger.info(f"Loading configuration from {config_path}")
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Validate configuration
        validate_config(config)
        
        logger.info("Configuration loaded successfully")
        return config
        
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML configuration: {e}")
        raise
    except Exception as e:
        logger.error(f"Error loading configuration: {e}")
        raise


def validate_config(config: Dict[str, Any]) -> None:
    """
    Validate configuration structure and required fields.
    
    Args:
        config: Configuration dictionary
        
    Raises:
        ValueError: If configuration is invalid
    """
    logger.info("Validating configuration")
    
    # Required top-level keys
    required_keys = ['species', 'go_obo', 'adaptive', 'output_dir']
    
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required configuration key: {key}")
    
    # Validate species configuration
    if not isinstance(config['species'], list) or len(config['species']) == 0:
        raise ValueError("Species configuration must be a non-empty list")
    
    for i, species in enumerate(config['species']):
        if not isinstance(species, dict):
            raise ValueError(f"Species {i} must be a dictionary")
        
        required_species_keys = ['code', 'name', 'genome_dir', 'gaf']
        for key in required_species_keys:
            if key not in species:
                raise ValueError(f"Missing required species key: {key}")
    
    # Validate adaptive configuration
    adaptive_config = config['adaptive']
    if not isinstance(adaptive_config, dict):
        raise ValueError("Adaptive configuration must be a dictionary")
    
    required_adaptive_keys = ['start_pct', 'step_pct', 'rounds']
    for key in required_adaptive_keys:
        if key not in adaptive_config:
            raise ValueError(f"Missing required adaptive key: {key}")
    
    # Validate numeric values
    try:
        start_pct = float(adaptive_config['start_pct'])
        step_pct = float(adaptive_config['step_pct'])
        rounds = int(adaptive_config['rounds'])
        
        if not (0 < start_pct <= 100):
            raise ValueError("start_pct must be between 0 and 100")
        
        if not (0 < step_pct <= 100):
            raise ValueError("step_pct must be between 0 and 100")
        
        if rounds < 1:
            raise ValueError("rounds must be at least 1")
        
        if start_pct - (rounds - 1) * step_pct <= 0:
            raise ValueError("Threshold would become negative or zero")
            
    except (ValueError, TypeError) as e:
        raise ValueError(f"Invalid numeric values in adaptive config: {e}")
    
    # Validate wobble configuration if present
    if 'wobble_aas' in config:
        wobble_aas = config['wobble_aas']
        if not isinstance(wobble_aas, list):
            raise ValueError("wobble_aas must be a list")
        
        valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
        for aa in wobble_aas:
            if aa not in valid_aas:
                raise ValueError(f"Invalid amino acid code: {aa}")
    
    logger.info("Configuration validation passed")


def get_default_config() -> Dict[str, Any]:
    """
    Get default configuration template.
    
    Returns:
        Dictionary with default configuration
    """
    return {
        'species': [
            {
                'code': 'example_species',
                'name': 'Example Species',
                'genome_dir': 'data/genomes/example_species/',
                'gaf': 'data/go/mappings/example_species.gaf'
            }
        ],
        'go_obo': 'data/go/go.obo',
        'adaptive': {
            'start_pct': 75,
            'step_pct': 10,
            'rounds': 3
        },
        'wobble_only': False,
        'wobble_aas': ['Leu', 'Lys', 'Gln', 'Glu', 'Phe', 'Trp'],
        'output_dir': 'results'
    }


def save_config(config: Dict[str, Any], output_path: str) -> None:
    """
    Save configuration to YAML file.
    
    Args:
        config: Configuration dictionary
        output_path: Path to save configuration
    """
    logger.info(f"Saving configuration to {output_path}")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    try:
        with open(output_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, indent=2)
        
        logger.info("Configuration saved successfully")
        
    except Exception as e:
        logger.error(f"Error saving configuration: {e}")
        raise


def merge_configs(base_config: Dict[str, Any], 
                 override_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge two configuration dictionaries.
    
    Args:
        base_config: Base configuration
        override_config: Override configuration
        
    Returns:
        Merged configuration dictionary
    """
    logger.info("Merging configurations")
    
    merged = base_config.copy()
    
    def _merge_dict(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
        result = base.copy()
        
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = _merge_dict(result[key], value)
            else:
                result[key] = value
        
        return result
    
    merged = _merge_dict(merged, override_config)
    
    # Validate merged configuration
    validate_config(merged)
    
    return merged


def expand_paths(config: Dict[str, Any], base_dir: str = '.') -> Dict[str, Any]:
    """
    Expand relative paths in configuration to absolute paths.
    
    Args:
        config: Configuration dictionary
        base_dir: Base directory for relative paths
        
    Returns:
        Configuration with expanded paths
    """
    logger.info("Expanding configuration paths")
    
    expanded = config.copy()
    
    # Expand GO ontology path
    if 'go_obo' in expanded:
        expanded['go_obo'] = os.path.abspath(
            os.path.join(base_dir, expanded['go_obo'])
        )
    
    # Expand output directory
    if 'output_dir' in expanded:
        expanded['output_dir'] = os.path.abspath(
            os.path.join(base_dir, expanded['output_dir'])
        )
    
    # Expand species paths
    if 'species' in expanded:
        for species in expanded['species']:
            if 'genome_dir' in species:
                species['genome_dir'] = os.path.abspath(
                    os.path.join(base_dir, species['genome_dir'])
                )
            if 'gaf' in species:
                species['gaf'] = os.path.abspath(
                    os.path.join(base_dir, species['gaf'])
                )
    
    return expanded


def load_wobble_list(file_path: str) -> List[str]:
    """
    Load wobble amino acid list from JSON or YAML file.
    
    Args:
        file_path: Path to wobble list file
        
    Returns:
        List of amino acid codes
    """
    logger.info(f"Loading wobble list from {file_path}")
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Wobble list file not found: {file_path}")
    
    try:
        with open(file_path, 'r') as f:
            if file_path.endswith('.json'):
                data = json.load(f)
            else:
                data = yaml.safe_load(f)
        
        if isinstance(data, list):
            wobble_aas = data
        elif isinstance(data, dict) and 'wobble_aas' in data:
            wobble_aas = data['wobble_aas']
        else:
            raise ValueError("Invalid wobble list format")
        
        # Validate amino acid codes
        valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
        for aa in wobble_aas:
            if aa not in valid_aas:
                raise ValueError(f"Invalid amino acid code: {aa}")
        
        logger.info(f"Loaded {len(wobble_aas)} wobble amino acids")
        return wobble_aas
        
    except Exception as e:
        logger.error(f"Error loading wobble list: {e}")
        raise


def create_example_config(output_path: str) -> None:
    """
    Create an example configuration file.
    
    Args:
        output_path: Path to save example configuration
    """
    logger.info(f"Creating example configuration at {output_path}")
    
    example_config = get_default_config()
    save_config(example_config, output_path)
    
    logger.info("Example configuration created successfully")


def get_species_config(config: Dict[str, Any], species_code: str) -> Optional[Dict[str, Any]]:
    """
    Get configuration for a specific species.
    
    Args:
        config: Full configuration dictionary
        species_code: Species code to look up
        
    Returns:
        Species configuration dictionary or None if not found
    """
    for species in config.get('species', []):
        if species.get('code') == species_code:
            return species
    
    return None


def list_species_codes(config: Dict[str, Any]) -> List[str]:
    """
    Get list of species codes from configuration.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        List of species codes
    """
    return [species.get('code') for species in config.get('species', [])]


def validate_file_paths(config: Dict[str, Any]) -> List[str]:
    """
    Validate that all file paths in configuration exist.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        List of missing file paths
    """
    logger.info("Validating file paths in configuration")
    
    missing_paths = []
    
    # Check GO ontology file
    go_obo = config.get('go_obo')
    if go_obo and not os.path.exists(go_obo):
        missing_paths.append(go_obo)
    
    # Check species files
    for species in config.get('species', []):
        genome_dir = species.get('genome_dir')
        if genome_dir and not os.path.exists(genome_dir):
            missing_paths.append(genome_dir)
        
        gaf_file = species.get('gaf')
        if gaf_file and not os.path.exists(gaf_file):
            missing_paths.append(gaf_file)
    
    if missing_paths:
        logger.warning(f"Missing file paths: {missing_paths}")
    else:
        logger.info("All file paths validated successfully")
    
    return missing_paths