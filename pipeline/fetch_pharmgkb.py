"""
PharmGKB Data Fetcher
----------------------
Downloads and parses PharmGKB clinical annotations — the curated database
of variant-drug-phenotype relationships drawn from published clinical studies.

What PharmGKB clinical annotations contain:
    - Variant/haplotype (e.g. CYP2D6 *4, rs3892097)
    - Drug (e.g. Codeine, Warfarin)
    - Phenotype (e.g. "Toxicity/ADR", "Efficacy", "Dosage")
    - Evidence level: 1A (FDA label) > 1B (replicated) > 2A/2B > 3 > 4
    - PMIDs of supporting publications

Why clinical annotations specifically:
    PharmGKB has multiple data tiers. Clinical annotations are the
    highest-confidence tier — each entry has been manually curated from
    primary literature and assigned an evidence grade. This is what
    pharmaceutical companies use for PGx label language decisions.

The core biological logic:
    Gene variants alter drug metabolism or target sensitivity.
    CYP2D6 poor metabolizers can't clear codeine -> morphine accumulates -> 
    respiratory depression. VKORC1 variants change warfarin sensitivity ->
    bleeding risk. SLCO1B1 variants reduce statin transport -> myopathy.
    
    These are mechanistic, not just correlational. The pharmacogenomic
    associations in PharmGKB have biological pathways behind them.
    FAERS signals don't — they're just pattern detection. Convergence
    of both = much stronger evidence.

PharmGKB licensing: Creative Commons Attribution-ShareAlike 4.0
Data: https://www.pharmgkb.org/downloads
"""

import argparse
import io
import logging
import zipfile
from pathlib import Path

import pandas as pd
import requests

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")

PHARMGKB_URLS = {
    "clinical_annotations": "https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip",
    "relationships":         "https://api.pharmgkb.org/v1/download/file/data/relationships.zip",
    "drugs":                 "https://api.pharmgkb.org/v1/download/file/data/drugs.zip",
    "genes":                 "https://api.pharmgkb.org/v1/download/file/data/genes.zip",
}

EVIDENCE_ORDER = {"1A": 6, "1B": 5, "2A": 4, "2B": 3, "3": 2, "4": 1}


def fetch_zip(url: str, label: str) -> dict:
    logger.info(f"Fetching {label}...")
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
    except requests.RequestException as e:
        logger.error(f"Failed: {e}")
        return {}
    tables = {}
    with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
        for name in zf.namelist():
            if name.endswith(".tsv"):
                with zf.open(name) as f:
                    try:
                        df = pd.read_csv(f, sep="\t", low_memory=False, on_bad_lines="warn")
                        tables[Path(name).stem] = df
                        logger.info(f"  {Path(name).name}: {len(df):,} rows")
                    except Exception as e:
                        logger.warning(f"  Could not parse {name}: {e}")
    return tables


def parse_clinical_annotations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean clinical annotations into analysis-ready format.

    The Phenotype Category field is critical for our analysis:
        Toxicity/ADR  -- adverse drug reactions (what we want to match to FAERS AEs)
        Efficacy      -- drug doesn't work as expected
        Dosage        -- variant affects required dose
        Metabolism/PK -- variant affects drug levels (often upstream of toxicity)
        Other

    Evidence levels explained:
        1A: Variant-drug association on FDA/EMA drug label
        1B: Replicated study, expert consensus
        2A: Replicated, but smaller studies
        2B: Single study, good design
        3:  Case reports or weak association
        4:  Preliminary/not replicated

    For enrichment analysis we focus on 1A/1B (strong evidence) and
    Toxicity/ADR category — these are the variants most likely to explain
    FAERS adverse event signals.
    """
    df = df.copy()
    df.columns = [c.strip().replace(" ", "_").replace("/", "_").upper() for c in df.columns]

    # Flexible column mapping (PharmGKB occasionally renames columns)
    col_map = {}
    for col in df.columns:
        if "LEVEL_OF_EVIDENCE" in col or col == "LEVEL":
            col_map[col] = "EVIDENCE_LEVEL"
        elif "DRUG" in col and col not in ("DRUG",):
            col_map[col] = "DRUGS_RAW"
        elif "PHENOTYPE" in col and "CATEGORY" not in col and col not in ("PHENOTYPE",):
            col_map[col] = "PHENOTYPES_RAW"
        elif "VARIANT" in col or "HAPLOTYPE" in col:
            col_map[col] = "VARIANTS"
        elif "PHENOTYPE_CATEGORY" in col:
            col_map[col] = "PHENOTYPE_CATEGORY"
        elif "GENE" == col:
            col_map[col] = "GENE"
    df = df.rename(columns=col_map)

    if "EVIDENCE_LEVEL" in df.columns:
        df["EVIDENCE_SCORE"] = df["EVIDENCE_LEVEL"].map(EVIDENCE_ORDER).fillna(0)

    # Explode multi-drug rows (PharmGKB uses "; " separator)
    drug_col = "DRUGS_RAW" if "DRUGS_RAW" in df.columns else None
    if drug_col:
        df["DRUG"] = df[drug_col].fillna("").str.split("; ")
        df = df.explode("DRUG")
        df["DRUG"] = df["DRUG"].str.strip().str.upper()

    # Normalize phenotype for fuzzy matching against MedDRA terms
    phen_col = "PHENOTYPES_RAW" if "PHENOTYPES_RAW" in df.columns else None
    if phen_col:
        df["PHENOTYPE_NORM"] = (df[phen_col]
                                .fillna("")
                                .str.upper()
                                .str.strip()
                                .str.replace(r"[^\w\s]", " ", regex=True)
                                .str.replace(r"\s+", " ", regex=True)
                                .str.strip())

    n_drugs = df["DRUG"].nunique() if "DRUG" in df.columns else "?"
    logger.info(f"Parsed {len(df):,} annotation rows | {n_drugs} unique drugs")
    return df


def save_pharmgkb(output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {}

    # Clinical annotations
    tables = fetch_zip(PHARMGKB_URLS["clinical_annotations"], "clinical_annotations")
    if tables:
        raw = list(tables.values())[0]
        raw.to_parquet(output_dir / "clinical_annotations_raw.parquet", index=False)
        parsed = parse_clinical_annotations(raw)
        parsed.to_parquet(output_dir / "clinical_annotations.parquet", index=False)
        results["clinical_annotations"] = parsed

    # Relationships
    tables = fetch_zip(PHARMGKB_URLS["relationships"], "relationships")
    if tables:
        raw = list(tables.values())[0]
        raw.to_parquet(output_dir / "relationships.parquet", index=False)
        results["relationships"] = raw

    # Drugs + genes (for name normalization)
    for key in ("drugs", "genes"):
        tables = fetch_zip(PHARMGKB_URLS[key], key)
        if tables:
            raw = list(tables.values())[0]
            raw.to_parquet(output_dir / f"{key}.parquet", index=False)
            results[key] = raw
            logger.info(f"Saved {key}: {len(raw):,} entries")

    logger.info(f"\nAll PharmGKB data saved to {output_dir}/")
    return results


def load_pharmgkb(data_dir: Path) -> dict:
    files = {
        "clinical_annotations": "clinical_annotations.parquet",
        "relationships":         "relationships.parquet",
        "drugs":                 "drugs.parquet",
        "genes":                 "genes.parquet",
    }
    result = {}
    for key, fname in files.items():
        path = data_dir / fname
        if path.exists():
            result[key] = pd.read_parquet(path)
        else:
            logger.warning(f"Missing: {path} - run fetch_pharmgkb.py first")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("data/pharmgkb"))
    parser.add_argument("--load-only", action="store_true")
    args = parser.parse_args()

    if args.load_only:
        data = load_pharmgkb(args.output_dir)
        for k, df in data.items():
            print(f"{k}: {len(df):,} rows | columns: {list(df.columns[:6])}")
    else:
        save_pharmgkb(args.output_dir)
