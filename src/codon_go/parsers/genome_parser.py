"""
Genome parser module for loading EMBL/GenBank annotations.
"""

import glob
import os
from typing import Dict, List, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import logging

logger = logging.getLogger(__name__)


def load_genome_annotations(genome_dir: str) -> Dict[str, SeqRecord]:
    """
    Glob all EMBL/GenBank files in genome_dir, extract CDS sequences,
    return dict gene_id → SeqRecord.
    
    Args:
        genome_dir: Directory containing EMBL/GenBank files
        
    Returns:
        Dictionary mapping gene IDs to SeqRecord objects
    """
    if not os.path.exists(genome_dir):
        raise FileNotFoundError(f"Genome directory not found: {genome_dir}")
    
    records = {}
    
    # Pattern to match EMBL and GenBank files
    patterns = [
        os.path.join(genome_dir, "*.embl"),
        os.path.join(genome_dir, "*.gb"),
        os.path.join(genome_dir, "*.gbk"),
        os.path.join(genome_dir, "*.genbank")
    ]
    
    files_found = []
    for pattern in patterns:
        files_found.extend(glob.glob(pattern))
    
    if not files_found:
        logger.warning(f"No EMBL/GenBank files found in {genome_dir}")
        return records
    
    logger.info(f"Found {len(files_found)} annotation files")
    
    for file_path in files_found:
        logger.info(f"Processing {file_path}")
        
        # Determine file format
        if file_path.endswith('.embl'):
            file_format = 'embl'
        else:
            file_format = 'genbank'
        
        try:
            for record in SeqIO.parse(file_path, file_format):
                for feature in record.features:
                    if feature.type == 'CDS':
                        gene_id = _extract_gene_id(feature)
                        if gene_id:
                            try:
                                # Extract CDS sequence
                                cds_seq = feature.extract(record.seq)
                                
                                # Create SeqRecord for the CDS
                                cds_record = SeqRecord(
                                    cds_seq,
                                    id=gene_id,
                                    description=f"CDS from {os.path.basename(file_path)}"
                                )
                                
                                records[gene_id] = cds_record
                                
                            except Exception as e:
                                logger.warning(f"Error extracting CDS for {gene_id}: {e}")
                                continue
                                
        except Exception as e:
            logger.error(f"Error parsing {file_path}: {e}")
            continue
    
    logger.info(f"Loaded {len(records)} CDS records")
    return records


def _extract_gene_id(feature: SeqFeature) -> Optional[str]:
    """
    Extract gene ID from a CDS feature.
    
    Tries multiple qualifier fields in order of preference:
    1. locus_tag
    2. gene
    3. protein_id
    4. db_xref (for systematic names, including CGD)
    
    Args:
        feature: Bio.SeqFeature.SeqFeature object
        
    Returns:
        Gene ID string or None if not found
    """
    # Priority order for gene ID extraction
    id_fields = ['locus_tag', 'gene', 'protein_id']
    
    for field in id_fields:
        if field in feature.qualifiers:
            return feature.qualifiers[field][0]
    
    # Check db_xref for systematic names (including CGD format)
    if 'db_xref' in feature.qualifiers:
        for xref in feature.qualifiers['db_xref']:
            if xref.startswith('GeneID:'):
                return xref.split(':')[1]
            elif xref.startswith('CGD:'):
                return xref.split(':')[1]
            elif xref.startswith('CAL'):  # Direct CGD systematic ID
                return xref
    
    return None


def _extract_all_gene_ids(feature: SeqFeature) -> Dict[str, str]:
    """
    Extract all possible gene IDs from a CDS feature for mapping purposes.
    
    Args:
        feature: Bio.SeqFeature.SeqFeature object
        
    Returns:
        Dictionary mapping ID type to ID value
    """
    ids = {}
    
    # Extract standard qualifiers
    id_fields = ['locus_tag', 'gene', 'protein_id']
    for field in id_fields:
        if field in feature.qualifiers:
            ids[field] = feature.qualifiers[field][0]
    
    # Extract db_xref entries
    if 'db_xref' in feature.qualifiers:
        for xref in feature.qualifiers['db_xref']:
            if ':' in xref:
                db_name, db_id = xref.split(':', 1)
                ids[f'db_xref_{db_name}'] = db_id
            else:
                ids['db_xref_other'] = xref
    
    return ids


def create_gene_id_mapping(records: Dict[str, SeqRecord]) -> Dict[str, str]:
    """
    Create a mapping from genome gene IDs to GAF-compatible IDs.
    
    This function analyzes all gene IDs in the genome records and attempts
    to create a mapping to GAF file compatible IDs.
    
    Args:
        records: Dictionary of gene_id → SeqRecord from genome files
        
    Returns:
        Dictionary mapping genome_gene_id → gaf_gene_id
    """
    logger.info("Creating gene ID mapping for GAF compatibility")
    
    # For now, return identity mapping - this can be enhanced based on
    # the specific ID formats found in the files
    mapping = {gene_id: gene_id for gene_id in records.keys()}
    
    # Log some examples for debugging
    sample_ids = list(records.keys())[:5]
    logger.info(f"Sample genome gene IDs: {sample_ids}")
    logger.info(f"Gene ID mapping created for {len(mapping)} genes")
    
    return mapping


def validate_cds_sequences(records: Dict[str, SeqRecord]) -> Dict[str, SeqRecord]:
    """
    Validate CDS sequences for proper length and start/stop codons.
    
    Args:
        records: Dictionary of gene_id → SeqRecord
        
    Returns:
        Filtered dictionary with valid CDS sequences
    """
    valid_records = {}
    
    for gene_id, record in records.items():
        seq = str(record.seq).upper()
        
        # Check if sequence length is multiple of 3
        if len(seq) % 3 != 0:
            logger.warning(f"CDS {gene_id} length not multiple of 3: {len(seq)}")
            continue
        
        # Check minimum length (at least one codon)
        if len(seq) < 3:
            logger.warning(f"CDS {gene_id} too short: {len(seq)} bp")
            continue
        
        # Check for valid start codon (ATG, GTG, TTG)
        start_codons = ['ATG', 'GTG', 'TTG']
        if seq[:3] not in start_codons:
            logger.warning(f"CDS {gene_id} invalid start codon: {seq[:3]}")
            continue
        
        # Check for valid stop codon
        stop_codons = ['TAA', 'TAG', 'TGA']
        if seq[-3:] not in stop_codons:
            logger.warning(f"CDS {gene_id} invalid stop codon: {seq[-3:]}")
            continue
        
        valid_records[gene_id] = record
    
    logger.info(f"Validated {len(valid_records)}/{len(records)} CDS sequences")
    return valid_records


def get_genome_stats(records: Dict[str, SeqRecord]) -> Dict[str, int]:
    """
    Get basic statistics about the genome annotation.
    
    Args:
        records: Dictionary of gene_id → SeqRecord
        
    Returns:
        Dictionary with genome statistics
    """
    if not records:
        return {}
    
    lengths = [len(record.seq) for record in records.values()]
    
    stats = {
        'total_genes': len(records),
        'total_bp': sum(lengths),
        'mean_length': sum(lengths) / len(lengths),
        'min_length': min(lengths),
        'max_length': max(lengths)
    }
    
    return stats