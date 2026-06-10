import sys
import pandas as pd
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path('../faers-pharmacovigilance')))
from pipeline import SignalDetector

pgx_raw = pd.read_parquet('data/pharmgkb/clinical_annotations_raw.parquet')
pgx_raw.columns = [c.strip().replace(" ","_").replace("/","_").upper() for c in pgx_raw.columns]

drug_col = [c for c in pgx_raw.columns if 'DRUG' in c][0]
pgx_drugs = set(
    pgx_raw[drug_col].fillna("").str.split(";").explode()
    .str.strip().str.upper()
)
pgx_drugs.discard("")
print(f"PharmGKB drugs: {len(pgx_drugs)}")

print("Loading FAERS master...")
drug = pd.read_parquet('../faers-pharmacovigilance/data/master/drug_master.parquet')
reac = pd.read_parquet('../faers-pharmacovigilance/data/master/reac_master.parquet')
drug["_drug"] = drug["PROD_AI"].str.upper().str.strip()
reac["_pt"]   = reac["PT"].str.upper().str.strip()

faers_drug_counts = drug[drug["ROLE_COD"]=="PS"]["_drug"].value_counts()
faers_in_pgx = [d for d in faers_drug_counts.index if d in pgx_drugs]
print(f"FAERS drugs also in PharmGKB: {len(faers_in_pgx)}")
for d in faers_in_pgx[:20]:
    print(f"  {d}: {faers_drug_counts[d]:,} reports")

target_drugs = [d for d in faers_in_pgx if faers_drug_counts.get(d,0) >= 100]
print(f"Running signals for {len(target_drugs)} drugs...")

detector = SignalDetector(drug, reac, drug_col="PROD_AI")
all_signals = []
for i, drug_name in enumerate(target_drugs):
    sigs = detector.run_all(drug_name, min_reports=3, signal_any=True)
    if not sigs.empty:
        sigs["DRUG"] = drug_name
        all_signals.append(sigs)
    if (i+1) % 10 == 0:
        print(f"  {i+1}/{len(target_drugs)} done...")

if all_signals:
    out = pd.concat(all_signals, ignore_index=True)
    out.to_csv("data/signals_pgx_drugs.csv", index=False)
    print(f"Saved {len(out):,} signals for {out['DRUG'].nunique()} drugs")
else:
    print("No signals found")
