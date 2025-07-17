#!/usr/bin/env python3
"""
Test script for the Codon-GO Analysis Pipeline.
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from codon_go.utils.config_loader import create_example_config, load_config
from codon_go.parsers.genome_parser import load_genome_annotations
from codon_go.analysis.codon_usage import compute_relative_usage_by_aa
from codon_go.viz.boxplots import create_codon_boxplot

def create_test_data():
    """Create minimal test data for demonstration."""
    print("Creating test data...")
    
    # Create test directories
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)
    
    genome_dir = test_dir / "genomes" / "test_species"
    genome_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a minimal EMBL file
    embl_content = """ID   TEST_CHR1; SV 1; linear; genomic DNA; STD; UNK; 1000 BP.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..1000
FT                   /organism="Test species"
FT   CDS             1..300
FT                   /locus_tag="TEST_001"
FT                   /gene="test_gene1"
FT                   /product="test protein 1"
FT                   /translation="MKLLVVVGGVGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDG
FT                   ETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDD
FT                   VPMVLVGNKSDLHEHPLLGSKVMSTLSAPLLVQFSSSFYGSIGQAHVAIIIYCRDNKQTM
FT                   QDDFFVRTSSLPEEVSELNEVLNKPRMRRNLVKRPQRQKLQNLFINFCLILICLLLNPVLA"
FT   CDS             400..600
FT                   /locus_tag="TEST_002"
FT                   /gene="test_gene2"
FT                   /product="test protein 2"
FT                   /translation="MKLVVVGGVGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDG
FT                   ETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDD
FT                   VPMVLVGNKSDLHEHPLLGSKVMSTLSAPLLVQFSSSFYGSIGQAHVAIIIYCRDNKQTM"
XX
SQ   Sequence 1000 BP; 250 A; 250 C; 250 G; 250 T; 0 other;
     atgaagctgc tggtggtggt gggcggcgtg ggcgtgggca agagcgcgct gaccatccag        60
     ctgatccaga accacttcgt ggacgagtac gacccgacca tcgaggacag ctaccgcaag       120
     caggtggtga tcgacggcga gacctgcctg ctggacatcc tggacacggc cggccaggag       180
     gagtacagcg cgatgcgcga ccagtacatg cgcaccggcg agggcttcct gtgcgtgttc       240
     gccatcaaca acaccaagag cttcgaggac atccaccagt accgcgagca gatcaagcgc       300
     gtgaaggaca gcgacgacgt gccgatggtg ctggtgggca acaagagcga cctgcacgag       360
     cacccgctgc tgggcagcaa ggtgatgagc accctgagcg cgccgctgct ggtgcagttc       420
     agcagcagct tctacggcag catcggccag gcccacgtgg ccatcatcat ctactgccgc       480
     gacaacaagc agaccatgca ggacgacttc ttcgtgcgca ccagcagcct gccggaggag       540
     gtgagcgagc tgaacgaggt gctgaacaag ccgcgcatgc gccgcaacct ggtgaagcgc       600
     ccgcagcgcc agaagctgca gaacctgttc atcaacttct gcctgatcct gatctgcctg       660
     ctgctgaacc cggtgctggc ctgaggcgag gcgctggtgg tggtgggcgg cgtgggcgtg       720
     ggcaagagcg cgctgaccat ccagctgatc cagaaccact tcgtggacga gtacgacccg       780
     accatcgagg acagctaccg caagcaggtg gtgatcgacg gcgagacctg cctgctggac       840
     atcctggaca cggccggcca ggaggagtac agcgcgatgc gcgaccagta catgcgcacc       900
     ggcgagggct tcctgtgcgt gttcgccatc aacaacacca agagcttcga ggacatccac       960
     cagtaccgcg agcagatcaa gcgcgtgaag gacagcgacg acgtgccgat ggtgctggtg      1020
//
"""
    
    with open(genome_dir / "chr1.embl", "w") as f:
        f.write(embl_content)
    
    # Create a minimal GAF file
    gaf_dir = test_dir / "go" / "mappings"
    gaf_dir.mkdir(parents=True, exist_ok=True)
    
    gaf_content = """!gaf-version: 2.2
TEST	TEST_001	test_gene1		GO:0005737	PMID:12345678	IEA		P	test protein 1		protein	taxon:12345	20240101	TEST
TEST	TEST_002	test_gene2		GO:0016740	PMID:12345678	IEA		F	test protein 2		protein	taxon:12345	20240101	TEST
"""
    
    with open(gaf_dir / "test_species.gaf", "w") as f:
        f.write(gaf_content)
    
    # Create a minimal GO ontology file
    go_obo_content = """format-version: 1.2
data-version: releases/2024-01-01
ontology: go

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component
def: "The contents of a cell excluding the plasma membrane and nucleus." [GOC:go_curators]

[Term]
id: GO:0016740
name: transferase activity
namespace: molecular_function
def: "Catalysis of the transfer of a group." [GOC:go_curators]
"""
    
    with open(test_dir / "go" / "go.obo", "w") as f:
        f.write(go_obo_content)
    
    print(f"Test data created in {test_dir}")
    return test_dir

def test_pipeline():
    """Test the main pipeline components."""
    print("Testing Codon-GO Analysis Pipeline...")
    
    # Create test data
    test_dir = create_test_data()
    
    try:
        # Test 1: Configuration
        print("\n1. Testing configuration...")
        config_path = test_dir / "config.yaml"
        create_example_config(str(config_path))
        
        # Modify config for test data
        config_content = f"""species:
  - code: test_species
    name: Test Species
    genome_dir: {test_dir}/genomes/test_species/
    gaf: {test_dir}/go/mappings/test_species.gaf
go_obo: {test_dir}/go/go.obo
adaptive:
  start_pct: 75
  step_pct: 10
  rounds: 3
wobble_only: false
wobble_aas:
  - Leu
  - Lys
output_dir: {test_dir}/results
"""
        
        with open(config_path, "w") as f:
            f.write(config_content)
        
        config = load_config(str(config_path))
        print("✓ Configuration loaded successfully")
        
        # Test 2: Genome parsing
        print("\n2. Testing genome parsing...")
        genome_dir = test_dir / "genomes" / "test_species"
        records = load_genome_annotations(str(genome_dir))
        print(f"✓ Loaded {len(records)} genome records")
        
        # Test 3: Codon usage analysis
        print("\n3. Testing codon usage analysis...")
        codon_usage_df = compute_relative_usage_by_aa(records)
        print(f"✓ Computed codon usage: {len(codon_usage_df)} observations")
        
        # Test 4: Visualization (if possible)
        print("\n4. Testing visualization...")
        if not codon_usage_df.empty:
            output_dir = test_dir / "results" / "figures"
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Test boxplot for one amino acid
            amino_acids = codon_usage_df['AA'].unique()
            if len(amino_acids) > 0:
                test_aa = amino_acids[0]
                output_path = output_dir / f"test_boxplot_{test_aa}.svg"
                create_codon_boxplot(codon_usage_df, test_aa, str(output_path))
                print(f"✓ Created test boxplot for {test_aa}")
        
        print("\n✅ All tests passed!")
        print(f"Test results saved in: {test_dir}")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

def main():
    """Main test function."""
    print("Codon-GO Analysis Pipeline Test")
    print("=" * 40)
    
    success = test_pipeline()
    
    if success:
        print("\n🎉 Pipeline test completed successfully!")
        print("\nNext steps:")
        print("1. Install dependencies: pip install -r requirements.txt")
        print("2. Prepare your data in the required format")
        print("3. Create configuration: python codon_go.py create-config")
        print("4. Run analysis: python codon_go.py run --config config/species.yaml")
    else:
        print("\n💥 Pipeline test failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()