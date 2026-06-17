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
  - Evidence levels: 1A (FDA label) > 1B (replicated study) > 2A/2B > 3 > 4

## Pipeline

```
faers_signals (ROR/PRR/EBGM per drug-AE pair, 252 drugs, 8 quarters)
        +
pharmgkb_annotations (variant → drug → phenotype, evidence grade)
        ↓
join on: drug name (normalized) + adverse event term (MedDRA PT ↔ PharmGKB phenotype)
         match confidence scored 0.5–1.0 (exact, synonym, token-overlap)
        ↓
convergence score = log₂(ROR) × evidence_weight × match_confidence
        ↓
enrichment test: Mann-Whitney U — are PGx-annotated AEs over-represented
                 among strong FAERS signals? (with BH-FDR + Bonferroni correction)
```

## Results

Enrichment tested across **147 drugs** present in both datasets.

| Multiple-testing correction | Significant drugs |
|---|:-:|
| None (raw p < 0.05) | 30 |
| Benjamini-Hochberg FDR (q < 0.05) | **11** |
| Bonferroni (p < 0.05) | **6** |

19 of 30 nominally-significant hits do not survive correction — the
stricter result is the honest one.

**Bonferroni-surviving drugs and mechanism:**

| Drug | Key gene(s) | Fold enrichment | Biological story |
|---|---|:-:|---|
| Letrozole | CYP19A1 / CYP2A6 | 6.2× | CYP19A1 encodes aromatase — the drug's target enzyme |
| Anastrozole | CYP19A1 class | 4.9× | Same aromatase inhibitor class; gene field unresolved |
| Leflunomide | ESR1 | 4.4× | Estrogen receptor variants modulate RA severity |
| Tocilizumab | FCGR3A | 3.3× | IgG1 mAb Fc receptor — same gene as rituximab |
| Bevacizumab | VEGFA / HSP90AB1 | 2.8× | VEGFA = both antibody target and susceptibility locus |
| Rituximab | FCGR3A / GSTA1 | 2.1× | ADCC efficiency via FcγRIIIA V158F; lymphoma/NMO axis |

**Evidence stratification (NB03):** Tier 1A matched pairs (FDA-label level) show
higher median ROR than tier 3 pairs (candidate-gene level), confirming that
PharmGKB evidence quality is calibrated to FAERS signal strength.

## Notebooks

| Notebook | Content |
|---|---|
| `01_pgx_enrichment.ipynb` | PharmGKB overview, drug-AE join, convergence scores, enrichment test with BH-FDR and Bonferroni correction |
| `02_top_hit_deep_dive.ipynb` | For each Bonferroni survivor: variant → predicted AE → observed FAERS ROR vs all other AEs. Rituximab full section (FCGR3A/ADCC biology). Aromatase inhibitor axis. |
| `03_evidence_stratification.ipynb` | ROR distributions stratified by PharmGKB evidence level (1A → 4). Tests whether annotation quality predicts signal strength — methodology calibration check. |

## Key Concepts

- **ROR (Reporting Odds Ratio):** How much more often is AE Y reported with Drug X
  than with everything else in the database. >1 = disproportionate reporting.
- **PharmGKB Evidence Level:** 1A (FDA label) > 1B (replicated study) > 2A/2B > 3 > 4
- **Convergence score:** log₂(ROR) × evidence_weight × match_confidence — prioritizes
  drug-AE pairs with strong spontaneous reporting AND strong genomic evidence AND
  reliable term matching.
- **Enrichment test:** Mann-Whitney U comparing ROR of PGx-annotated vs unannotated
  AEs per drug. Significant = annotated AEs are over-represented among strong signals.

## Requirements

```
pip install pandas numpy scipy requests matplotlib seaborn pyarrow statsmodels
```

## Usage

```bash
# Fetch latest PharmGKB annotations
python pipeline/fetch_pharmgkb.py

# Regenerate FAERS signal file (if signals_pgx_drugs.csv is missing)
python pgx_expand.py

# Analysis notebooks (run in order)
jupyter notebook notebooks/01_pgx_enrichment.ipynb
jupyter notebook notebooks/02_top_hit_deep_dive.ipynb
jupyter notebook notebooks/03_evidence_stratification.ipynb
```
