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
    return dict gene_id -> SeqRecord.
    
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
                                
                                # Log if we found a CAL* ID (good for GAF matching)
                                if gene_id.startswith('CAL'):
                                    logger.debug(f"Found CAL* ID: {gene_id}")
                                
                            except Exception as e:
                                logger.warning(f"Error extracting CDS for {gene_id}: {e}")
                                continue
                                
        except Exception as e:
            logger.error(f"Error parsing {file_path}: {e}")
            continue
    
    # Analyze gene ID types for diagnostic purposes
    cal_ids = sum(1 for gene_id in records.keys() if gene_id.startswith('CAL'))
    systematic_ids = sum(1 for gene_id in records.keys() if '_' in gene_id and not gene_id.startswith('CAL'))
    locus_tag_ids = sum(1 for gene_id in records.keys() if any(c.isdigit() for c in gene_id) and not gene_id.startswith('CAL') and '_' not in gene_id)
    other_ids = len(records) - cal_ids - systematic_ids - locus_tag_ids
    
    logger.info(f"Loaded {len(records)} CDS records:")
    logger.info(f"  - CAL* database IDs: {cal_ids} (matches GAF files)")
    logger.info(f"  - Systematic names: {systematic_ids} (e.g., C1_00010W_A)")
    logger.info(f"  - Locus tags: {locus_tag_ids} (e.g., gene001, orf19.123)")
    logger.info(f"  - Other IDs: {other_ids}")
    
    # Show samples of each type
    if cal_ids > 0:
        sample_cal = [gene_id for gene_id in list(records.keys())[:10] if gene_id.startswith('CAL')][:3]
        logger.info(f"  - Sample CAL* IDs: {sample_cal}")
    
    if systematic_ids > 0:
        sample_sys = [gene_id for gene_id in list(records.keys())[:10] if '_' in gene_id and not gene_id.startswith('CAL')][:3]
        logger.info(f"  - Sample systematic names: {sample_sys}")
    
    if locus_tag_ids > 0:
        sample_locus = [gene_id for gene_id in list(records.keys())[:10] if any(c.isdigit() for c in gene_id) and not gene_id.startswith('CAL') and '_' not in gene_id][:3]
        logger.info(f"  - Sample locus tags: {sample_locus}")
    
    # Warn if we have mostly locus tags (might need GAF file with locus_tag format)
    if locus_tag_ids > cal_ids and locus_tag_ids > systematic_ids:
        logger.info("  📍 NOTICE: Most gene IDs are locus_tag format")
        logger.info("  📍 Ensure your GAF file uses corresponding identifiers for matching")
    
    return records


def _extract_gene_id(feature: SeqFeature) -> Optional[str]:
    """
    Extract gene ID from a CDS feature.
    
    Smart priority system for CGD/Candida gene IDs:
    1. id (if it looks like CAL* database identifier - matches GAF files)
    2. locus_tag (systematic name like C1_00010W_A - ALWAYS as fallback)
    3. id (if not CAL* but still useful)
    4. gene
    5. protein_id
    6. db_xref (for systematic names, including CGD)
    
    Args:
        feature: Bio.SeqFeature.SeqFeature object
        
    Returns:
        Gene ID string or None if not found
    """
    # Collect all available IDs
    available_ids = {}
    
    # Standard fields
    id_fields = ['id', 'locus_tag', 'gene', 'protein_id']
    for field in id_fields:
        if field in feature.qualifiers:
            available_ids[field] = feature.qualifiers[field][0]
    
    # db_xref fields
    if 'db_xref' in feature.qualifiers:
        for xref in feature.qualifiers['db_xref']:
            if xref.startswith('GeneID:'):
                available_ids['geneid'] = xref.split(':')[1]
            elif xref.startswith('CGD:'):
                available_ids['cgd'] = xref.split(':')[1]
            elif xref.startswith('CAL'):  # Direct CGD systematic ID
                available_ids['cal_xref'] = xref
    
    # PRIORITY 1: Prefer /id field if it looks like CAL* database ID (matches GAF)
    if 'id' in available_ids:
        id_value = available_ids['id']
        if id_value.startswith('CAL') and len(id_value) > 10:  # CAL identifiers are long
            logger.debug(f"Using /id field (CAL* database ID): {id_value}")
            return id_value
    
    # PRIORITY 2: ALWAYS try locus_tag as fallback (systematic names)
    if 'locus_tag' in available_ids:
        locus_tag = available_ids['locus_tag']
        logger.debug(f"Using locus_tag (systematic name): {locus_tag}")
        return locus_tag
    
    # PRIORITY 3: Use /id field even if not CAL* (might still be useful)
    if 'id' in available_ids:
        id_value = available_ids['id']
        logger.debug(f"Using /id field (non-CAL): {id_value}")
        return id_value
    
    # PRIORITY 4-5: Other standard fields
    for field in ['gene', 'protein_id']:
        if field in available_ids:
            gene_id = available_ids[field]
            logger.debug(f"Using {field}: {gene_id}")
            return gene_id
    
    # PRIORITY 6: db_xref entries
    for xref_field in ['cgd', 'geneid', 'cal_xref']:
        if xref_field in available_ids:
            gene_id = available_ids[xref_field]
            logger.debug(f"Using {xref_field}: {gene_id}")
            return gene_id
    
    logger.warning(f"No suitable gene ID found in qualifiers: {list(available_ids.keys())}")
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
    
    # Extract standard qualifiers including the critical /id field
    id_fields = ['id', 'locus_tag', 'gene', 'protein_id']
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
        records: Dictionary of gene_id -> SeqRecord from genome files
        
    Returns:
        Dictionary mapping genome_gene_id -> gaf_gene_id
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
        records: Dictionary of gene_id -> SeqRecord
        
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
        records: Dictionary of gene_id -> SeqRecord
        
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