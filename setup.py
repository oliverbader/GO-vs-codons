#!/usr/bin/env python3
"""
Setup script for the Codon-GO Analysis Pipeline.
"""

from setuptools import setup, find_packages
import os

# Read README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

# Read requirements
def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(requirements_path):
        with open(requirements_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return []

setup(
    name="codon-go",
    version="1.0.0",
    author="Codon-GO Pipeline",
    author_email="codon-go@example.com",
    description="A modular Python pipeline for analyzing codon usage and GO term enrichment in eukaryotic genomes",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/example/codon-go",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    entry_points={
        "console_scripts": [
            "codon-go=codon_go.cli:cli",
        ],
    },
    include_package_data=True,
    package_data={
        "codon_go": ["data/*"],
    },
    zip_safe=False,
)