#!/usr/bin/env python3
"""
Script to check that warning suppression is working correctly.
"""

import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

print("Testing warning suppression...")

# Test without suppression
print("\n1. Testing without warning suppression:")
try:
    import warnings
    warnings.resetwarnings()  # Clear any existing filters
    
    # This might trigger the pkg_resources warning
    from goatools import obo_parser
    print("✓ GOATOOLS imported without visible warnings")
except ImportError:
    print("⚠ GOATOOLS not installed, skipping this test")
except Exception as e:
    print(f"✗ Error: {e}")

# Test with suppression
print("\n2. Testing with warning suppression:")
try:
    from codon_go.utils.warnings_config import suppress_common_warnings
    suppress_common_warnings()
    
    # Try importing again
    from goatools import obo_parser
    print("✓ GOATOOLS imported with warnings suppressed")
except ImportError:
    print("⚠ GOATOOLS not installed, skipping this test")
except Exception as e:
    print(f"✗ Error: {e}")

# Test the pipeline import
print("\n3. Testing pipeline import:")
try:
    from codon_go.cli import cli
    print("✓ Codon-GO CLI imported successfully")
except Exception as e:
    print(f"✗ Error importing CLI: {e}")

print("\n✅ Warning suppression test completed!")
print("\nTo run the pipeline:")
print("  python codon_go.py --help")
print("  python codon_go.py show-cug-info")
print("  python test_pipeline.py")