"""
File utilities module for common file operations.
"""

import os
import shutil
from typing import Dict, Any, List, Optional
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def ensure_directory(directory: str) -> None:
    """
    Ensure a directory exists, creating it if necessary.
    
    Args:
        directory: Directory path to create
    """
    if not os.path.exists(directory):
        logger.info(f"Creating directory: {directory}")
        os.makedirs(directory, exist_ok=True)
    else:
        logger.debug(f"Directory already exists: {directory}")


def save_dataframe(df: pd.DataFrame, 
                  output_path: str, 
                  format: str = 'tsv',
                  **kwargs) -> None:
    """
    Save DataFrame to file in specified format.
    
    Args:
        df: DataFrame to save
        output_path: Output file path
        format: Output format ('tsv', 'csv', 'xlsx', 'parquet')
        **kwargs: Additional arguments for pandas save methods
    """
    logger.info(f"Saving DataFrame to {output_path} in {format} format")
    
    # Ensure output directory exists
    ensure_directory(os.path.dirname(output_path))
    
    try:
        if format.lower() == 'tsv':
            df.to_csv(output_path, sep='\t', index=False, **kwargs)
        elif format.lower() == 'csv':
            df.to_csv(output_path, index=False, **kwargs)
        elif format.lower() == 'xlsx':
            df.to_excel(output_path, index=False, **kwargs)
        elif format.lower() == 'parquet':
            df.to_parquet(output_path, index=False, **kwargs)
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        logger.info(f"DataFrame saved successfully: {len(df)} rows")
        
    except Exception as e:
        logger.error(f"Error saving DataFrame: {e}")
        raise


def load_dataframe(file_path: str, 
                  format: Optional[str] = None,
                  **kwargs) -> pd.DataFrame:
    """
    Load DataFrame from file, auto-detecting format if not specified.
    
    Args:
        file_path: Path to input file
        format: File format ('tsv', 'csv', 'xlsx', 'parquet')
        **kwargs: Additional arguments for pandas read methods
        
    Returns:
        Loaded DataFrame
    """
    logger.info(f"Loading DataFrame from {file_path}")
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Auto-detect format if not specified
    if format is None:
        _, ext = os.path.splitext(file_path)
        format_map = {
            '.tsv': 'tsv',
            '.csv': 'csv',
            '.xlsx': 'xlsx',
            '.parquet': 'parquet'
        }
        format = format_map.get(ext.lower(), 'tsv')
    
    try:
        if format.lower() == 'tsv':
            df = pd.read_csv(file_path, sep='\t', **kwargs)
        elif format.lower() == 'csv':
            df = pd.read_csv(file_path, **kwargs)
        elif format.lower() == 'xlsx':
            df = pd.read_excel(file_path, **kwargs)
        elif format.lower() == 'parquet':
            df = pd.read_parquet(file_path, **kwargs)
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        logger.info(f"DataFrame loaded successfully: {len(df)} rows, {len(df.columns)} columns")
        return df
        
    except Exception as e:
        logger.error(f"Error loading DataFrame: {e}")
        raise


def backup_file(file_path: str, backup_suffix: str = '.bak') -> str:
    """
    Create a backup copy of a file.
    
    Args:
        file_path: Path to file to backup
        backup_suffix: Suffix to add to backup file
        
    Returns:
        Path to backup file
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    backup_path = file_path + backup_suffix
    logger.info(f"Creating backup: {file_path} -> {backup_path}")
    
    shutil.copy2(file_path, backup_path)
    return backup_path


def clean_directory(directory: str, 
                   pattern: str = '*',
                   keep_subdirs: bool = True) -> int:
    """
    Clean files from a directory matching a pattern.
    
    Args:
        directory: Directory to clean
        pattern: File pattern to match (glob style)
        keep_subdirs: Whether to keep subdirectories
        
    Returns:
        Number of files removed
    """
    import glob
    
    if not os.path.exists(directory):
        logger.warning(f"Directory does not exist: {directory}")
        return 0
    
    logger.info(f"Cleaning directory: {directory} (pattern: {pattern})")
    
    files_removed = 0
    search_pattern = os.path.join(directory, pattern)
    
    for file_path in glob.glob(search_pattern):
        if os.path.isfile(file_path):
            try:
                os.remove(file_path)
                files_removed += 1
                logger.debug(f"Removed file: {file_path}")
            except Exception as e:
                logger.warning(f"Could not remove file {file_path}: {e}")
        elif os.path.isdir(file_path) and not keep_subdirs:
            try:
                shutil.rmtree(file_path)
                files_removed += 1
                logger.debug(f"Removed directory: {file_path}")
            except Exception as e:
                logger.warning(f"Could not remove directory {file_path}: {e}")
    
    logger.info(f"Cleaned {files_removed} items from {directory}")
    return files_removed


def get_file_size(file_path: str) -> int:
    """
    Get file size in bytes.
    
    Args:
        file_path: Path to file
        
    Returns:
        File size in bytes
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    return os.path.getsize(file_path)


def format_file_size(size_bytes: int) -> str:
    """
    Format file size in human-readable format.
    
    Args:
        size_bytes: Size in bytes
        
    Returns:
        Formatted size string
    """
    if size_bytes == 0:
        return "0 B"
    
    size_names = ["B", "KB", "MB", "GB", "TB"]
    i = 0
    while size_bytes >= 1024 and i < len(size_names) - 1:
        size_bytes /= 1024.0
        i += 1
    
    return f"{size_bytes:.1f} {size_names[i]}"


def list_files(directory: str, 
               extensions: Optional[List[str]] = None,
               recursive: bool = False) -> List[str]:
    """
    List files in a directory with optional filtering.
    
    Args:
        directory: Directory to search
        extensions: List of file extensions to include (e.g., ['.txt', '.csv'])
        recursive: Whether to search recursively
        
    Returns:
        List of file paths
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    files = []
    
    if recursive:
        for root, dirs, filenames in os.walk(directory):
            for filename in filenames:
                file_path = os.path.join(root, filename)
                if extensions is None or any(filename.endswith(ext) for ext in extensions):
                    files.append(file_path)
    else:
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                if extensions is None or any(filename.endswith(ext) for ext in extensions):
                    files.append(file_path)
    
    return sorted(files)


def copy_file(src: str, dst: str, overwrite: bool = False) -> None:
    """
    Copy a file from source to destination.
    
    Args:
        src: Source file path
        dst: Destination file path
        overwrite: Whether to overwrite existing destination
    """
    if not os.path.exists(src):
        raise FileNotFoundError(f"Source file not found: {src}")
    
    if os.path.exists(dst) and not overwrite:
        raise FileExistsError(f"Destination file exists: {dst}")
    
    logger.info(f"Copying file: {src} -> {dst}")
    
    # Ensure destination directory exists
    ensure_directory(os.path.dirname(dst))
    
    shutil.copy2(src, dst)
    logger.debug(f"File copied successfully")


def move_file(src: str, dst: str, overwrite: bool = False) -> None:
    """
    Move a file from source to destination.
    
    Args:
        src: Source file path
        dst: Destination file path
        overwrite: Whether to overwrite existing destination
    """
    if not os.path.exists(src):
        raise FileNotFoundError(f"Source file not found: {src}")
    
    if os.path.exists(dst) and not overwrite:
        raise FileExistsError(f"Destination file exists: {dst}")
    
    logger.info(f"Moving file: {src} -> {dst}")
    
    # Ensure destination directory exists
    ensure_directory(os.path.dirname(dst))
    
    shutil.move(src, dst)
    logger.debug(f"File moved successfully")


def create_symlink(src: str, dst: str, overwrite: bool = False) -> None:
    """
    Create a symbolic link from source to destination.
    
    Args:
        src: Source file path
        dst: Destination link path
        overwrite: Whether to overwrite existing destination
    """
    if not os.path.exists(src):
        raise FileNotFoundError(f"Source file not found: {src}")
    
    if os.path.exists(dst) and not overwrite:
        raise FileExistsError(f"Destination link exists: {dst}")
    
    logger.info(f"Creating symlink: {src} -> {dst}")
    
    # Remove existing link if overwrite is True
    if os.path.exists(dst) and overwrite:
        os.remove(dst)
    
    # Ensure destination directory exists
    ensure_directory(os.path.dirname(dst))
    
    os.symlink(src, dst)
    logger.debug(f"Symlink created successfully")


def get_directory_size(directory: str) -> int:
    """
    Get total size of a directory and its contents.
    
    Args:
        directory: Directory path
        
    Returns:
        Total size in bytes
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    total_size = 0
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if os.path.exists(file_path):
                total_size += os.path.getsize(file_path)
    
    return total_size


def compress_directory(directory: str, 
                      output_path: str,
                      format: str = 'zip') -> None:
    """
    Compress a directory to an archive.
    
    Args:
        directory: Directory to compress
        output_path: Output archive path
        format: Archive format ('zip', 'tar', 'gztar', 'bztar', 'xztar')
    """
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    logger.info(f"Compressing directory {directory} to {output_path}")
    
    # Ensure output directory exists
    ensure_directory(os.path.dirname(output_path))
    
    # Remove extension from output path for shutil.make_archive
    base_name = os.path.splitext(output_path)[0]
    
    try:
        shutil.make_archive(base_name, format, directory)
        logger.info(f"Directory compressed successfully")
    except Exception as e:
        logger.error(f"Error compressing directory: {e}")
        raise


def extract_archive(archive_path: str, 
                   extract_to: str,
                   format: Optional[str] = None) -> None:
    """
    Extract an archive to a directory.
    
    Args:
        archive_path: Path to archive file
        extract_to: Directory to extract to
        format: Archive format (auto-detected if None)
    """
    if not os.path.exists(archive_path):
        raise FileNotFoundError(f"Archive not found: {archive_path}")
    
    logger.info(f"Extracting archive {archive_path} to {extract_to}")
    
    # Ensure extraction directory exists
    ensure_directory(extract_to)
    
    try:
        shutil.unpack_archive(archive_path, extract_to, format)
        logger.info(f"Archive extracted successfully")
    except Exception as e:
        logger.error(f"Error extracting archive: {e}")
        raise