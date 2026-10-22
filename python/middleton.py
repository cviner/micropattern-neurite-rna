"""Middleton et al. (2019) Additional file 6 gene lists and the overlap test applied to them."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd
from scipy.stats import hypergeom

from common import EXTERNAL

TABLES = EXTERNAL / "Middleton2019_MOESM6.xlsx"
EXTERNAL_STUDIES = ("Ainsley2014", "Cajigas2012", "Lein2007", "Poon2006",
                    "Taliaferro2016", "Tushev2018", "Zhong2006")


def _symbols(values: pd.Series) -> set[str]:
    return {str(v).strip() for v in values.dropna() if str(v).strip()}


def dendrite_upregulated() -> set[str]:
    """Table S5: the 387 dendrite-upregulated genes. The sheet has no header row."""
    return _symbols(pd.read_excel(TABLES, sheet_name="TableS5-deDend genes", header=None)[0])


def consistently_dendritic() -> set[str]:
    """Table S6: the 1,863 consistently localized genes. The sheet has no header row."""
    return _symbols(pd.read_excel(TABLES, sheet_name="TableS6-consDend genes", header=None)[0])


def published_catalogue(external_only: bool = False) -> set[str]:
    """Table S4. With external_only, keep genes reported by any of the seven external studies,
    which drops genes contributed only by Middleton et al.'s own data ("ThisStudy")."""
    s4 = pd.read_excel(TABLES, sheet_name="TableS4-PublishedDendriticGenes")
    if external_only:
        found = s4[list(EXTERNAL_STUDIES)].apply(lambda c: c.astype(str).str.strip().str.lower() == "yes")
        s4 = s4[found.any(axis=1)]
    return _symbols(s4["GeneName"])


def overlap(query: Iterable[str], reference: Iterable[str], universe: Iterable[str]) -> dict:
    """One-sided hypergeometric (Fisher) over-representation of query in reference within universe."""
    universe = set(universe)
    a, b = set(query) & universe, set(reference) & universe
    n, k, total = len(a), len(b), len(universe)
    hits = a & b
    expected = n * k / total
    return {"universe": total, "query": n, "reference": k, "overlap": len(hits),
            "pct_of_query": 100 * len(hits) / n, "pct_of_reference": 100 * len(hits) / k,
            "expected": expected, "fold": len(hits) / expected,
            "P": float(hypergeom.sf(len(hits) - 1, total, k, n)), "genes": sorted(hits)}
