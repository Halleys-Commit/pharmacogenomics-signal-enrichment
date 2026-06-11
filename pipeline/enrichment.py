"""
Pharmacogenomic Signal Enrichment Analysis
-------------------------------------------
Joins FAERS adverse event signals with PharmGKB pharmacogenomic annotations
to identify drug-AE pairs supported by BOTH spontaneous reporting AND
genomic evidence.

The Matching Problem:
    FAERS uses MedDRA Preferred Terms (standardized): "HAEMORRHAGE", "MYOPATHY"
    PharmGKB uses free-text phenotype descriptions: "Bleeding", "muscle toxicity"
    
    These don't match exactly. We use three matching strategies in order:
    
    1. Exact match (after normalization) -- most reliable
    2. Token overlap -- share >=2 meaningful words
    3. Curated synonym map -- known equivalents (e.g. "bleeding" = "haemorrhage")
    
    Each match gets a confidence score. Exact=1.0, token=0.5-0.9, synonym=0.8.
    We report match confidence alongside the enrichment result so you know
    how much to trust each joined pair.

The Convergence Score:
    For each drug-AE pair that matches across both datasets:
    
    convergence_score = log2(ROR) * evidence_weight * match_confidence
    
    Where evidence_weight = {1A: 1.0, 1B: 0.8, 2A: 0.6, 2B: 0.4, 3: 0.2, 4: 0.1}
    
    High convergence score = strong spontaneous signal + strong genomic evidence
    + reliable term matching. These are the most interesting drug-AE pairs.

The Enrichment Test:
    For each drug in both datasets, we ask:
    "Are the AEs with PharmGKB annotations showing higher ROR than
     the AEs without PharmGKB annotations?"
    
    Statistical test: Mann-Whitney U (non-parametric, ROR is not normally distributed)
    A significant result means: yes, genomically-annotated AEs are
    disproportionately represented among strong signals for this drug.
    This is the key finding — it validates that FAERS captures pharmacogenomic
    signal even without individual genotype data.
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

# ── Known MedDRA PT <-> PharmGKB phenotype synonym mappings ──────────────────
# Manually curated based on common mismatches. Expand as you find more.
# Format: {pharmgkb_term_fragment: [meddra_pt_fragments]}
SYNONYM_MAP = {
    "BLEEDING":         ["HAEMORRHAGE", "HEMORRHAGE", "BLOOD LOSS"],
    "BLEED":            ["HAEMORRHAGE", "HEMORRHAGE"],
    "MUSCLE TOXICITY":  ["MYOPATHY", "MYALGIA", "RHABDOMYOLYSIS"],
    "MYOTOXICITY":      ["MYOPATHY", "MYALGIA", "RHABDOMYOLYSIS", "CREATINE PHOSPHOKINASE INCREASED"],
    "HEPATOTOXICITY":   ["HEPATITIS", "LIVER INJURY", "HEPATIC FAILURE", "ALANINE AMINOTRANSFERASE INCREASED"],
    "LIVER TOXICITY":   ["HEPATITIS", "LIVER INJURY", "HEPATIC FAILURE"],
    "NEPHROTOXICITY":   ["RENAL FAILURE", "ACUTE KIDNEY INJURY", "CREATININE INCREASED"],
    "CARDIOTOXICITY":   ["CARDIAC FAILURE", "ARRHYTHMIA", "QT PROLONGATION"],
    "QT PROLONGATION":  ["QT PROLONGED", "LONG QT SYNDROME"],
    "RESPIRATORY DEPRESSION": ["RESPIRATORY FAILURE", "APNOEA", "HYPOXIA"],
    "NEUROTOXICITY":    ["NEUROPATHY", "PERIPHERAL NEUROPATHY", "ENCEPHALOPATHY"],
    "THROMBOCYTOPENIA": ["PLATELET COUNT DECREASED", "THROMBOCYTOPENIA"],
    "NEUTROPENIA":      ["NEUTROPHIL COUNT DECREASED", "NEUTROPENIA", "FEBRILE NEUTROPENIA"],
    "HYPOGLYCAEMIA":    ["HYPOGLYCAEMIA", "BLOOD GLUCOSE DECREASED", "HYPOGLYCEMIC COMA"],
    "HYPERSENSITIVITY": ["ANAPHYLACTIC REACTION", "ANAPHYLAXIS", "DRUG HYPERSENSITIVITY"],
    "SKIN TOXICITY":    ["RASH", "DERMATITIS", "STEVENS-JOHNSON SYNDROME", "TOXIC EPIDERMAL NECROLYSIS"],
    "PERIPHERAL NEUROPATHY": ["NEUROPATHY PERIPHERAL", "PERIPHERAL SENSORY NEUROPATHY"],
}


class PGxEnrichment:
    """
    Join FAERS signals with PharmGKB annotations and compute enrichment.

    Parameters
    ----------
    signals : pd.DataFrame
        Output of SignalDetector.run_all() — must have columns:
        DRUG, PT, ROR, ROR_CI_lo, ROR_CI_hi, PRR, EBGM, SIGNAL_COUNT, a
    pgx : pd.DataFrame
        Parsed PharmGKB clinical_annotations — from fetch_pharmgkb.py
    min_evidence_score : int
        Minimum PharmGKB evidence score to include (1=4, 2=3, 3=2B, 4=2A, 5=1B, 6=1A)
        Default 2 = grade 3 and above (excludes weakest evidence)
    """

    EVIDENCE_WEIGHTS = {"1A": 1.0, "1B": 0.8, "2A": 0.6, "2B": 0.4, "3": 0.2, "4": 0.1}

    def __init__(
        self,
        signals: pd.DataFrame,
        pgx: pd.DataFrame,
        min_evidence_score: int = 2,
    ):
        self.signals = signals.copy()
        self.pgx = pgx[pgx.get("EVIDENCE_SCORE", 0) >= min_evidence_score].copy()

        # Normalize drug names in both datasets
        if "DRUG" in self.signals.columns:
            self.signals["_drug_norm"] = self.signals["DRUG"].str.upper().str.strip()
        if "DRUG" in self.pgx.columns:
            self.pgx["_drug_norm"] = self.pgx["DRUG"].str.upper().str.strip()

        logger.info(f"Signals: {len(self.signals):,} rows, {self.signals.get('_drug_norm', pd.Series()).nunique()} drugs")
        logger.info(f"PGx annotations (evidence >= {min_evidence_score}): {len(self.pgx):,} rows")

    def _match_terms(self, meddra_pt: str, pgx_phenotype: str) -> float:
        """
        Score how well a MedDRA PT matches a PharmGKB phenotype string.
        Returns confidence score 0.0-1.0.
        
        0.0 = no match
        0.5-0.9 = token overlap (scaled by overlap fraction)  
        0.8 = synonym map match
        1.0 = exact match after normalization
        """
        pt  = str(meddra_pt).upper().strip()
        pgx = str(pgx_phenotype).upper().strip()

        if not pt or not pgx:
            return 0.0

        # Exact match
        if pt == pgx:
            return 1.0
        if pt in pgx or pgx in pt:
            return 0.9

        # Synonym map check
        for pgx_fragment, meddra_fragments in SYNONYM_MAP.items():
            if pgx_fragment in pgx:
                for mf in meddra_fragments:
                    if mf in pt:
                        return 0.8

        # Token overlap (ignore stopwords)
        stopwords = {"AND", "OR", "OF", "TO", "IN", "THE", "A", "AN", "WITH", "FOR"}
        pt_tokens  = set(pt.split()) - stopwords
        pgx_tokens = set(pgx.split()) - stopwords

        if not pt_tokens or not pgx_tokens:
            return 0.0

        overlap = len(pt_tokens & pgx_tokens)
        if overlap == 0:
            return 0.0

        # Scale: need at least 1 shared token, score by Jaccard similarity
        jaccard = overlap / len(pt_tokens | pgx_tokens)
        return round(min(0.85, jaccard * 2), 2)  # cap at 0.85, boost overlap

    def join(self, min_match_confidence: float = 0.5) -> pd.DataFrame:
        """
        Join FAERS signals with PharmGKB annotations by drug name + phenotype matching.

        For each (drug, AE) pair in FAERS signals, find PharmGKB annotations
        for the same drug, then score how well the phenotype descriptions match.

        Returns one row per matched (FAERS signal, PharmGKB annotation) pair,
        with match confidence and convergence score.

        Parameters
        ----------
        min_match_confidence : float
            Minimum phenotype match confidence to include. Default 0.5.
            Raise to 0.8+ for higher-confidence matches only.
        """
        if "_drug_norm" not in self.signals.columns or "_drug_norm" not in self.pgx.columns:
            raise ValueError("Drug normalization failed — check column names")

        # Get drugs present in both datasets
        faers_drugs = set(self.signals["_drug_norm"].dropna())
        pgx_drugs   = set(self.pgx["_drug_norm"].dropna())
        common_drugs = faers_drugs & pgx_drugs

        if not common_drugs:
            logger.warning("No drug name overlap between FAERS signals and PharmGKB. "
                          "Check drug name normalization — FAERS uses active ingredient names "
                          "(PROD_AI), PharmGKB uses brand/generic. May need fuzzy matching.")
            logger.info(f"Sample FAERS drugs: {list(faers_drugs)[:10]}")
            logger.info(f"Sample PGx drugs:   {list(pgx_drugs)[:10]}")
            return pd.DataFrame()

        logger.info(f"Drugs in both datasets: {len(common_drugs)} — {sorted(common_drugs)[:10]}...")

        rows = []
        pgx_phen_col = "PHENOTYPE_NORM" if "PHENOTYPE_NORM" in self.pgx.columns else None

        for drug in common_drugs:
            drug_signals = self.signals[self.signals["_drug_norm"] == drug]
            drug_pgx     = self.pgx[self.pgx["_drug_norm"] == drug]

            for _, sig_row in drug_signals.iterrows():
                pt = str(sig_row.get("PT", ""))
                best_match  = None
                best_conf   = 0.0
                best_pgx_row = None

                for _, pgx_row in drug_pgx.iterrows():
                    pgx_phen = str(pgx_row.get(pgx_phen_col or "PHENOTYPES_RAW", ""))
                    conf = self._match_terms(pt, pgx_phen)
                    if conf > best_conf:
                        best_conf    = conf
                        best_match   = pgx_phen
                        best_pgx_row = pgx_row

                if best_conf >= min_match_confidence and best_pgx_row is not None:
                    evidence_level = str(best_pgx_row.get("EVIDENCE_LEVEL", ""))
                    evidence_weight = self.EVIDENCE_WEIGHTS.get(evidence_level, 0.1)
                    ror = float(sig_row.get("ROR", 1))

                    rows.append({
                        # FAERS signal fields
                        "DRUG":          drug,
                        "PT":            pt,
                        "n_reports":     int(sig_row.get("a", 0)),
                        "ROR":           ror,
                        "ROR_CI_lo":     float(sig_row.get("ROR_CI_lo", np.nan)),
                        "ROR_CI_hi":     float(sig_row.get("ROR_CI_hi", np.nan)),
                        "PRR":           float(sig_row.get("PRR", np.nan)),
                        "EBGM":          float(sig_row.get("EBGM", np.nan)),
                        "SIGNAL_COUNT":  int(sig_row.get("SIGNAL_COUNT", 0)),
                        # PharmGKB fields
                        "PGX_PHENOTYPE":     best_match,
                        "PGX_GENE":          str(best_pgx_row.get("GENE", "")),
                        "PGX_VARIANTS":      str(best_pgx_row.get("VARIANTS", "")),
                        "EVIDENCE_LEVEL":    evidence_level,
                        "EVIDENCE_SCORE":    float(best_pgx_row.get("EVIDENCE_SCORE", 0)),
                        "PHENOTYPE_CATEGORY":str(best_pgx_row.get("PHENOTYPE_CATEGORY", "")),
                        # Match quality
                        "MATCH_CONFIDENCE":  best_conf,
                        # Convergence score — the key metric
                        # High = strong spontaneous signal + strong genomic evidence + good term match
                        "CONVERGENCE_SCORE": round(
                            np.log2(max(ror, 1.01)) * evidence_weight * best_conf, 3
                        ),
                    })

        result = pd.DataFrame(rows)
        if result.empty:
            logger.warning("No matches found — consider lowering min_match_confidence")
            return result

        result = result.sort_values("CONVERGENCE_SCORE", ascending=False).reset_index(drop=True)
        logger.info(f"Matched {len(result):,} drug-AE pairs across {result['DRUG'].nunique()} drugs")
        return result

    def enrichment_test(self, joined: pd.DataFrame) -> pd.DataFrame:
        """
        For each drug, test whether pharmacogenomically-annotated AEs show
        higher ROR than non-annotated AEs.

        Method: Mann-Whitney U test comparing ROR distributions of:
            - AEs WITH PharmGKB annotation (pgx_annotated=True)
            - AEs WITHOUT PharmGKB annotation (pgx_annotated=False)

        A significant result (p < 0.05) means: for this drug, the AEs that
        have known pharmacogenomic variants are reporting at higher rates than
        the AEs without genetic annotations. This is the population-level
        genomic signal we're trying to detect.

        Mann-Whitney U is appropriate here because:
        - ROR is not normally distributed (right-skewed, bounded at 0)
        - Sample sizes per drug are often small
        - We want to test stochastic dominance, not mean difference
        """
        if joined.empty or "DRUG" not in joined.columns:
            return pd.DataFrame()

        results = []
        pgx_drugs = set(joined["DRUG"].unique())

        for drug in pgx_drugs:
            # Annotated AEs: those in the joined table for this drug
            pgx_aes  = set(joined[joined["DRUG"]==drug]["PT"])
            drug_sigs = self.signals[self.signals.get("_drug_norm", pd.Series()) == drug]

            if len(drug_sigs) < 5:
                continue

            annotated     = drug_sigs[drug_sigs["PT"].isin(pgx_aes)]["ROR"].dropna()
            not_annotated = drug_sigs[~drug_sigs["PT"].isin(pgx_aes)]["ROR"].dropna()

            if len(annotated) < 2 or len(not_annotated) < 2:
                continue

            u_stat, p_val = stats.mannwhitneyu(
                annotated, not_annotated, alternative="greater"
            )

            results.append({
                "DRUG":                   drug,
                "n_pgx_annotated_aes":    len(annotated),
                "n_unannotated_aes":      len(not_annotated),
                "median_ror_annotated":   round(annotated.median(), 2),
                "median_ror_unannotated": round(not_annotated.median(), 2),
                "ror_fold_enrichment":    round(annotated.median() / max(not_annotated.median(), 0.01), 2),
                "mannwhitney_u":          round(u_stat, 1),
                "p_value":                round(p_val, 4),
                "significant":            p_val < 0.05,
                "interpretation": (
                    "Genomically-annotated AEs show significantly higher reporting rates"
                    if p_val < 0.05 else
                    "No significant enrichment detected"
                )
            })

        result = pd.DataFrame(results).sort_values("p_value").reset_index(drop=True)

        if not result.empty:
            _, p_fdr_bh, _, _ = multipletests(result["p_value"], method="fdr_bh",    alpha=0.05)
            _, p_bonf,   _, _ = multipletests(result["p_value"], method="bonferroni", alpha=0.05)
            result["p_value_fdr_bh"]        = p_fdr_bh.round(4)
            result["p_value_bonferroni"]     = p_bonf.round(4)
            result["significant_fdr"]        = p_fdr_bh < 0.05
            result["significant_bonferroni"] = p_bonf < 0.05

        n_sig = result["significant"].sum() if not result.empty else 0
        logger.info(f"Enrichment test: {n_sig}/{len(result)} drugs show significant PGx enrichment")
        return result
