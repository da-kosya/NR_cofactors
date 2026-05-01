# NR_cofactors

Structural bioinformatics pipeline for analysing nuclear receptor (NR)
cofactor complexes from the Protein Data Bank, conducted at Nazarbayev
University under Dr. Ferdinand Molnár (BIOL363).

## Research context

Nuclear receptors are a family of 48 transcription factors that regulate
gene expression in response to ligands (hormones, vitamins, metabolites).
This pipeline systematically retrieves and classifies all available PDB
structures for human NRs, filters for ligand-binding domain (LBD)
structures co-crystallised with cofactor peptides, and extracts sequence
and structural features for downstream analysis.

## Pipeline

| Step | Description |
|------|-------------|
| 1 | Fetch UniProt accessions for 48 human NR genes via UniProt REST API |
| 2 | Extract annotated LBD sequences and check for zinc finger domains |
| 3 | Retrieve all associated PDB structure IDs via RCSB Search API |
| 4 | Download PDB XML files and classify as DBD (Zinc present) or LBD |
| 5 | Filter LBD structures for NR–cofactor peptide complexes |
| 6 | Visualise DBD/LBD distribution and chain composition |

## Outputs

- `nuclear_receptor_LBD_analysis.txt` — per-receptor LBD summary table
- `nuclear_receptor_LBDs.fasta` — extracted LBD sequences in FASTA format
- `nr_structural_classification.txt` — DBD vs LBD counts per receptor
- `nr_peptide_complexes.txt` — PDB IDs of valid NR–cofactor complexes
- `nr_peptide_complexes_detailed.txt` — chain-level breakdown of each complex
- `dbd_vs_lbd.png` — stacked bar chart of structure counts
- `chain_distribution.png` — polypeptide chain count distribution

## Usage

1. Set `BASE_DIR` in `nr_cofactors.py` to your local data directory
2. Install dependencies
3. Run the pipeline

```bash
pip install -r requirements.txt
python nr_cofactors.py
```

## Files
