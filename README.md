# pharmacogenomics-signal-enrichment

Cross-referencing FDA FAERS pharmacovigilance signals with PharmGKB 
pharmacogenomic variant-drug-phenotype associations to identify convergent 
evidence of genotype-driven adverse event risk.

## The Scientific Question

FDA FAERS contains millions of spontaneous adverse event reports — no genetic 
data. PharmGKB contains curated evidence that specific gene variants predict 
drug response and toxicity — no real-world reporting volumes.

**This pipeline asks:** Do drugs with known pharmacogenomic risk variants show 
disproportionately stronger FAERS signals for their genotype-predicted adverse 
events compared to their other AEs?

Convergent signal from two independent data sources (spontaneous reporting + 
genomic evidence) strengthens the biological plausibility of an AE without 
requiring individual-level genotype data.

## What This Is (and Isn't)

**Is:** Population-level ecological inference. If Drug X + Variant Z predicts 
AE Y in carriers, and Drug X shows elevated ROR for AE Y in FAERS, that is 
convergent evidence worth investigating.

**Is not:** Individual-level causal inference. FAERS reporters have no genotype 
data. We are identifying statistical enrichment at the population level — the 
same epistemological framework as GWAS.

## Data Sources

- **FAERS:** FDA Adverse Event Reporting System — public quarterly ASCII files
- **PharmGKB:** pharmgkb.org — Creative Commons licensed, free API + bulk download
  - `clinical_annotations.tsv` — curated variant-drug-phenotype associations with evidence grades
  - `relationships.tsv` — gene-drug-disease network edges
  - `drugs.tsv`, `genes.tsv` — entity metadata

## Pipeline

```
faers_signals (ROR/PRR/EBGM per drug-AE pair)
        +
pharmgkb_annotations (variant → drug → phenotype, evidence grade)
        ↓
join on: drug name (normalized) + adverse event term (MedDRA PT ↔ PharmGKB phenotype)
        ↓
enrichment analysis: are pharmacogenomically-annotated AEs over-represented 
                     among strong FAERS signals?
        ↓
visualization + report
```

## Key Concepts

- **ROR (Reporting Odds Ratio):** How much more often is AE Y reported with Drug X 
  than with everything else in the database. >1 = disproportionate reporting.
- **PharmGKB Evidence Level:** 1A (FDA label) > 1B (replicated study) > 2A/2B > 3 > 4
- **Convergence score:** Combined metric weighting ROR signal strength × PharmGKB 
  evidence level — prioritizes drug-AE pairs with strong spontaneous reporting AND 
  strong genomic evidence.

## Requirements

```
pip install pandas numpy scipy requests matplotlib seaborn pyarrow
```

## Usage

```bash
# Fetch latest PharmGKB annotations
python pipeline/fetch_pharmgkb.py

# Run enrichment against FAERS signals
python pipeline/enrichment.py --signals ../faers-pharmacovigilance/data/master/signals_master.csv

# Full analysis notebook
jupyter notebook notebooks/01_pgx_enrichment.ipynb
```
