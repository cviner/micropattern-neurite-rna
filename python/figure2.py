#!/usr/bin/env python3
"""Figure 2 panel data under results/figure2/. R/figure2.R draws the figure from these tables.

a  mean distal versus mean proximal CPM over the 12,396 expressed genes
b  per-set-rate permutation null and the naive symmetric null for the number of genes distal-higher
   in >= 10 of 12 matched sets
c  neurite-mean CPM of the 54 candidates against the expressed universe
d  the 54 by curated category (Supplementary File 4) and log2 mean distal/proximal ratio

Reads results/nulls/per_set_rate_null.tsv, so run proximal_distal_nulls.py first.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import DATA, RESULTS
from middleton import published_catalogue
from neurite import THRESHOLD, load

SUPP4 = DATA / "published" / "SuppFile4(Distal_Proximal).xlsx"
N_SYMMETRIC = 5_000
CATEGORIES = ("Not assigned", "RNA-regulatory", "Cytoskeletal", "Signaling/synaptic",
              "Mitochondrial", "Translation/ribosome")


def panel_data() -> dict[str, pd.DataFrame]:
    nd = load()
    p, d = nd.pair_cpm()
    mean_p, mean_d = p.mean(axis=1), d.mean(axis=1)
    candidates = nd.candidates()
    is54 = nd.universe.isin(candidates)

    scatter = pd.DataFrame({"gene": nd.universe, "prox": np.log10(mean_p + 1),
                            "dist": np.log10(mean_d + 1), "is54": is54})

    # Seed 20260625. Draw 5,000 binomial replicates at rate 0.35 before the symmetric null.
    nontied = (nd.directions() != 0).sum(axis=1)
    rs = np.random.RandomState(20260625)
    for _ in range(N_SYMMETRIC):
        rs.binomial(nontied, 0.35)
    symmetric = np.array([(rs.binomial(nontied, 0.5) >= THRESHOLD).sum() for _ in range(N_SYMMETRIC)])
    per_set = pd.read_csv(RESULTS / "nulls" / "per_set_rate_null.tsv", sep="\t")["distal_higher_ge10of12"]
    nulls = pd.concat([pd.DataFrame({"null": "per_set_rate", "count": per_set}),
                       pd.DataFrame({"null": "symmetric", "count": symmetric})], ignore_index=True)

    abundance = nd.neurite_mean()
    detected = (nd.counts.loc[candidates, nd.dist] >= 1).sum(axis=1)
    abund = pd.concat([pd.DataFrame({"group": "universe", "cpm": abundance[nd.universe]}),
                       pd.DataFrame({"group": "54", "cpm": abundance[candidates],
                                     "detected_distal_sets": detected})], ignore_index=True)
    abund["log10_cpm"] = np.log10(abund.cpm + 1)

    pub = pd.read_excel(SUPP4, sheet_name="54 candidate genes").dropna(subset=["Distal-higher sets (of 12)"])
    category = pub.set_index("Gene symbol")["Curated category"].replace(
        {"Not assigned in this curation": "Not assigned"})
    cat = pd.DataFrame({"gene": nd.universe[is54],
                        "log2fc": np.log2((mean_d[is54] + 1) / (mean_p[is54] + 1))})
    cat["category"] = cat.gene.map(category)
    cat["in_middleton_s4"] = cat.gene.isin(published_catalogue())
    assert cat.category.isin(CATEGORIES).all()
    return {"scatter": scatter, "nulls": nulls, "abundance": abund, "categories": cat}


def main() -> None:
    out = RESULTS / "figure2"
    out.mkdir(parents=True, exist_ok=True)
    for name, frame in panel_data().items():
        frame.to_csv(out / f"{name}.tsv", sep="\t", index=False)
        print(out / f"{name}.tsv")


if __name__ == "__main__":
    main()
