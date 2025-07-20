"""
Cache management and dependency tracking for the Codon-GO pipeline.

Provides intelligent caching to avoid recomputing expensive operations
and enables incremental analysis updates.
"""

import os
import pickle
import hashlib
import json
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
import pandas as pd
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


class CacheManager:
    """Manages caching and dependency tracking for pipeline stages."""
    
    def __init__(self, base_dir: str, species_code: str):
        """
        Initialize cache manager.
        
        Args:
            base_dir: Base directory for cache files
            species_code: Species identifier
        """
        self.base_dir = Path(base_dir)
        self.species_code = species_code
        self.cache_dir = self.base_dir / 'cache' / species_code
        self.processed_dir = self.base_dir / 'processed'
        self.figures_dir = self.base_dir / 'figures'
        
        # Create directories
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.processed_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
        # Cache metadata file
        self.metadata_file = self.cache_dir / 'cache_metadata.json'
        self.metadata = self._load_metadata()
    
    def _load_metadata(self) -> Dict:
        """Load cache metadata."""
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        return {}
    
    def _save_metadata(self):
        """Save cache metadata."""
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
    
    def _get_file_hash(self, filepath: str) -> str:
        """Get MD5 hash of a file."""
        hash_md5 = hashlib.md5()
        try:
            with open(filepath, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hash_md5.update(chunk)
            return hash_md5.hexdigest()
        except FileNotFoundError:
            return ""
    
    def _get_files_hash(self, filepaths: List[str]) -> str:
        """Get combined hash of multiple files."""
        combined_hash = hashlib.md5()
        for filepath in sorted(filepaths):  # Sort for consistent ordering
            file_hash = self._get_file_hash(filepath)
            combined_hash.update(file_hash.encode())
        return combined_hash.hexdigest()
    
    def check_stage_cache(self, stage_name: str, input_files: List[str], 
                         output_files: List[str], force: bool = False) -> bool:
        """
        Check if a pipeline stage can be skipped based on cache.
        
        Args:
            stage_name: Name of the pipeline stage
            input_files: List of input file paths
            output_files: List of expected output file paths
            force: Force recomputation even if cache is valid
            
        Returns:
            True if stage can be skipped, False if it needs to run
        """
        if force:
            logger.info(f"Stage {stage_name}: Forced recomputation (--force)")
            return False
        
        # Check if all output files exist
        missing_outputs = []
        for output_file in output_files:
            if not Path(output_file).exists():
                missing_outputs.append(output_file)
        
        if missing_outputs:
            logger.info(f"Stage {stage_name}: Missing outputs {missing_outputs}")
            return False
        
        # Check if input files have changed since last run
        current_input_hash = self._get_files_hash(input_files)
        
        stage_key = f"{self.species_code}_{stage_name}"
        if stage_key in self.metadata:
            cached_input_hash = self.metadata[stage_key].get('input_hash', '')
            cached_timestamp = self.metadata[stage_key].get('timestamp', '')
            
            if current_input_hash == cached_input_hash:
                logger.info(f"Stage {stage_name}: Cache valid, skipping (cached: {cached_timestamp})")
                return True
            else:
                logger.info(f"Stage {stage_name}: Input files changed, recomputing")
                return False
        else:
            logger.info(f"Stage {stage_name}: No cache entry found, computing")
            return False
    
    def update_stage_cache(self, stage_name: str, input_files: List[str], 
                          output_files: List[str]):
        """
        Update cache metadata after successful stage completion.
        
        Args:
            stage_name: Name of the pipeline stage
            input_files: List of input file paths
            output_files: List of output file paths
        """
        input_hash = self._get_files_hash(input_files)
        timestamp = datetime.now().isoformat()
        
        stage_key = f"{self.species_code}_{stage_name}"
        self.metadata[stage_key] = {
            'input_hash': input_hash,
            'input_files': input_files,
            'output_files': output_files,
            'timestamp': timestamp
        }
        
        self._save_metadata()
        logger.info(f"Stage {stage_name}: Cache updated ({timestamp})")
    
    def save_pickle(self, obj: Any, filename: str) -> str:
        """
        Save object to pickle file in cache directory.
        
        Args:
            obj: Object to save
            filename: Filename (without path)
            
        Returns:
            Full path to saved file
        """
        filepath = self.cache_dir / filename
        with open(filepath, 'wb') as f:
            pickle.dump(obj, f)
        logger.debug(f"Saved pickle: {filepath}")
        return str(filepath)
    
    def load_pickle(self, filename: str) -> Any:
        """
        Load object from pickle file in cache directory.
        
        Args:
            filename: Filename (without path)
            
        Returns:
            Loaded object
        """
        filepath = self.cache_dir / filename
        with open(filepath, 'rb') as f:
            obj = pickle.load(f)
        logger.debug(f"Loaded pickle: {filepath}")
        return obj
    
    def get_cache_status(self) -> Dict:
        """
        Get comprehensive cache status for this species.
        
        Returns:
            Dictionary with cache information
        """
        status = {
            'species': self.species_code,
            'cache_dir': str(self.cache_dir),
            'stages': {}
        }
        
        for stage_key, stage_data in self.metadata.items():
            if stage_key.startswith(f"{self.species_code}_"):
                stage_name = stage_key.replace(f"{self.species_code}_", "")
                
                # Check if outputs still exist
                outputs_exist = all(Path(f).exists() for f in stage_data.get('output_files', []))
                
                status['stages'][stage_name] = {
                    'timestamp': stage_data.get('timestamp'),
                    'outputs_exist': outputs_exist,
                    'output_files': stage_data.get('output_files', []),
                    'input_files': stage_data.get('input_files', [])
                }
        
        return status
    
    def clear_cache(self, stage_name: Optional[str] = None):
        """
        Clear cache for specific stage or all stages.
        
        Args:
            stage_name: Specific stage to clear, or None for all stages
        """
        if stage_name:
            stage_key = f"{self.species_code}_{stage_name}"
            if stage_key in self.metadata:
                # Remove output files
                stage_data = self.metadata[stage_key]
                for output_file in stage_data.get('output_files', []):
                    try:
                        Path(output_file).unlink()
                        logger.info(f"Removed cached file: {output_file}")
                    except FileNotFoundError:
                        pass
                
                # Remove metadata entry
                del self.metadata[stage_key]
                self._save_metadata()
                logger.info(f"Cleared cache for stage: {stage_name}")
        else:
            # Clear all caches for this species
            species_keys = [k for k in self.metadata.keys() if k.startswith(f"{self.species_code}_")]
            for stage_key in species_keys:
                stage_data = self.metadata[stage_key]
                for output_file in stage_data.get('output_files', []):
                    try:
                        Path(output_file).unlink()
                    except FileNotFoundError:
                        pass
                del self.metadata[stage_key]
            
            # Remove cache directory
            import shutil
            if self.cache_dir.exists():
                shutil.rmtree(self.cache_dir)
            
            self._save_metadata()
            logger.info(f"Cleared all cache for species: {self.species_code}")


def get_expected_outputs(species_code: str, base_dir: str, config: Dict) -> Dict[str, List[str]]:
    """
    Get expected output files for each pipeline stage.
    
    Args:
        species_code: Species identifier
        base_dir: Base output directory
        config: Pipeline configuration
        
    Returns:
        Dictionary mapping stage names to expected output files
    """
    processed_dir = Path(base_dir) / 'processed'
    figures_dir = Path(base_dir) / 'figures'
    cache_dir = Path(base_dir) / 'cache' / species_code
    
    # Get codons and categories from config
    codons = ['TTT', 'TTC', 'TTA', 'TTG']  # Will be dynamic based on analysis
    categories = ['BP', 'MF', 'CC']
    
    outputs = {
        'genome_loading': [
            str(cache_dir / f'{species_code}_genome_records.pkl')
        ],
        'codon_usage': [
            str(processed_dir / f'{species_code}_codon_usage.tsv')
        ],
        'go_loading': [
            str(processed_dir / f'{species_code}_gene2go.tsv')
        ],
        'adaptive_analysis': [
            str(processed_dir / f'{species_code}_adaptive.tsv'),
            str(processed_dir / f'{species_code}_diagnostic.tsv'),
            str(processed_dir / f'{species_code}_summary.tsv')
        ],
        'classical_enrichment': [],
        'visualizations': [
            str(figures_dir / f'boxplot_comprehensive_{species_code.lower()}.svg'),
            str(figures_dir / f'heatmap_adaptive.svg')
        ]
    }
    
    # Add classical enrichment outputs (will be generated dynamically)
    for codon in codons:
        for category in categories:
            outputs['classical_enrichment'].append(
                str(processed_dir / f'{species_code}_classical_{category}_{codon}.tsv')
            )
            outputs['visualizations'].append(
                str(figures_dir / f'go_dotplot_{codon}_{category}_{species_code.lower()}.svg')
            )
    
    return outputs