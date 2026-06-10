import pandas as pd
from pathlib import Path

pgx = pd.read_parquet('data/pharmgkb/clinical_annotations.parquet')

# Fix: re-explode on ";" instead of "; "
pgx["DRUG"] = pgx["DRUGS_RAW"].fillna("").str.split(";")
pgx = pgx.explode("DRUG")
pgx["DRUG"] = pgx["DRUG"].str.strip().str.upper()
pgx = pgx[pgx["DRUG"] != ""]

print(f"After fix: {len(pgx):,} rows | {pgx['DRUG'].nunique()} unique drugs")
print(f"\nSample drugs: {sorted(pgx['DRUG'].unique())[:20]}")

# Check overlap with FAERS signals
signals = pd.read_csv('../faers-pharmacovigilance/data/master/signals_master.csv')
print(f"\nFAERS signals: {len(signals):,} rows | {signals['DRUG'].nunique()} drugs")
print(f"FAERS drugs: {sorted(signals['DRUG'].unique())}")

faers_drugs = set(signals['DRUG'].str.upper().str.strip())
pgx_drugs   = set(pgx['DRUG'].str.upper().str.strip())
overlap = faers_drugs & pgx_drugs
print(f"\nDirect name overlap: {overlap}")

# Also check toxicity annotations specifically
tox = pgx[pgx['PHENOTYPE_CATEGORY'].str.contains('Toxicity', na=False)]
print(f"\nToxicity/ADR annotations: {len(tox):,}")
print(f"Toxicity drugs: {sorted(tox['DRUG'].unique())[:30]}")
