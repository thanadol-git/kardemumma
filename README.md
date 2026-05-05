# Blood Proteome

A bioinformatics pipeline to compile, query, and export human blood/plasma proteome data from public sources (UniProt, PRIDE).

## Project structure

```
blood_proteome/
├── summon.ipynb                      # Main workflow: aggregate plasma proteome data and export FASTA
├── targeted.ipynb                    # Build FASTA from targeted PRM results
├── src/blood_proteome/
│   ├── __init__.py                   # Re-exports all public symbols
│   └── fn.py                         # UniProt API helpers and FASTA utilities
├── pyproject.toml                    # Package metadata (pip install -e .)
├── config.yml                        # Conda environment definition
├── export/                           # Output FASTA files
│   ├── Blood_proteome_*.fasta
│   ├── Human_proteome_All_*.fasta
│   └── Human_proteome_SwissProt_*.fasta
└── targeted/
    └── Transition Results.csv        # Skyline PRM transition export
```

## Setup

```bash
# 1. Create and activate the conda environment (includes pip install -e .)
conda env create -f config.yml -p ./env
conda activate ./env

# 2. Register the Jupyter kernel
python -m ipykernel install --user --name blood_proteome
```

## Workflows

### 1. Blood / plasma proteome (`summon.ipynb`)

Reads protein group files from Geyer et al. (2019), queries UniProt for sequences and metadata, and writes a FASTA file to `export/`.

### 2. Targeted PRM FASTA (`targeted.ipynb`)

Reads Skyline transition results, extracts UniProt accessions from protein names, fetches sequences in batch, and exports a targeted FASTA for spectral library building.

### 3. Helper library (`src/blood_proteome/fn.py`)

| Function | Description |
|---|---|
| `get_url` | Thin wrapper around `requests.get` with error handling |
| `get_protein_sequences_batch` | Batch-query UniProt for a list of accessions |
| `get_swissprot_sequences_batch` | Same as above, Swiss-Prot reviewed entries only |
| `query_human_proteome` | Stream the full human proteome from UniProt |
| `create_human_proteome_fasta` | Download and write a human proteome FASTA |
| `validate_fasta` | Print summary statistics for a FASTA file |
