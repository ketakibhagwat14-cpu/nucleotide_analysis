# Nucleotide Composition Analyser

This script analyses the coding sequence (CDS) of multiple genes of interest 
fetched directly from NCBI. For each gene, it computes the absolute count of 
each nucleotide (A, T, G, C) and the A/T and G/C ratios, returning the results 
as a clean table.

## Requirements

- Python 3
- Biopython
- pandas

Install dependencies with:

pip install biopython pandas

## How to Use

1. Open the script and add your NCBI accession IDs to the `ACCESSION_IDS` list:

ACCESSION_IDS = ["NM_002232.5", "NM_000546.6", "NM_007294.3"]

2. Run the script:

python nucleotide_analyser.py

Or run it cell by cell in a Jupyter notebook.

## Output

A table with one row per accession ID showing:
- Individual nucleotide counts (A, T, G, C)
- A/T ratio
- G/C ratio

## Dependencies

- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
