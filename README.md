# Genomic Tools

## Overview
This repository contains various scripts for genomic data processing and analysis.

## Files

### extract_seq_from_genome_based_on_gff_feature.py
This script extracts sequences from a genome based on GFF features.

### format_genbank.py
This script contains functions to format sequences into GenBank files.

#### Functions:
- `makegb(**kwargs) -> str`: Creates a GenBank file from the provided keyword arguments.
- `help() -> None`: Prints the documentation for the `makegb` function.

### gen_gel_map.py
This script generates fragment sizes for gel electrophoresis maps.

#### Classes:
- `gen_frag_sizes`: Generates fragment sizes for given restriction enzymes.

### gen_seq_from_motifs.py
This script generates sequences from motif files.

#### Functions:
- `gen_seq_from_motifs(motif_files: str|list[str]) -> dict[str,float]`: Generates sequences from motifs.

### generate_CAS9_primers.py
This script helps generate CAS9 primers.

#### Functions:
- `help()`: Prints the documentation for the script.

### gff_genome_to_annfa.py
This script converts GFF genome annotations to ANNFA format.

### guide_helper.py
This script contains helper functions for guide RNA design.

#### Classes:
- `CAS`: Contains methods for CAS system operations.

#### Functions:
- `help()`: Prints the documentation for the script.

### plot_gel.py
This script plots gel electrophoresis results.

#### Functions:
- `plot_gel(ax, file_path, restriction_enzymes, index=0)`: Plots the gel electrophoresis results.

### select_genes.py
This script selects genes based on specific features from a genome.

### smith_waterman.py
This script implements the Smith-Waterman algorithm for sequence alignment.

#### Functions:
- `parse_fasta(file)`: Parses a FASTA file and returns the sequences.

### __init__.py
This file initializes the package.

### __pycache__
This directory contains the compiled bytecode of the Python files.

### .gitignore
This file specifies the files and directories to be ignored by Git.

### README.md
This file contains the documentation for the project.