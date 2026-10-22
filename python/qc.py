#!/usr/bin/env python3
"""Sample correlations, outlier screen, and library sizes."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import COMPARTMENTS, cpm, libraries, load_counts, write_results


def whole_control_mean_r(c: pd.DataFrame) -> pd.Series:
    cols = libraries(c, "Whole") + libraries(c, "Control")
    r = c[cols].corr().to_numpy().copy()
    np.fill_diagonal(r, np.nan)
    return pd.Series(np.nanmean(r, axis=1), index=cols)


def compartment_screen(c: pd.DataFrame) -> pd.DataFrame:
    """Within-compartment mean Pearson r on CPM and PC1-PC2 Mahalanobis distance."""
    rows = []
    for comp in COMPARTMENTS:
        cols = libraries(c, comp)
        sub = c[cols]
        r = sub.corr().to_numpy().copy()
        np.fill_diagonal(r, np.nan)
        mean_r = np.nanmean(r, axis=1)
        log = np.log2(sub.to_numpy() + 1)
        _, _, vh = np.linalg.svd(log - log.mean(axis=1, keepdims=True), full_matrices=False)
        pcs = vh[:2].T
        delta = pcs - pcs.mean(axis=0)
        md = np.sqrt(np.einsum("ij,jk,ik->i", delta, np.linalg.inv(np.cov(pcs.T)), delta))
        q1, q3 = np.percentile(mean_r, [25, 75])
        low_r = q1 - 1.5 * (q3 - q1)
        high_md = np.median(md) + 2.5 * np.std(md)
        for i, s in enumerate(cols):
            rows.append(dict(library=s, compartment=comp, mean_r=mean_r[i], mahalanobis=md[i],
                             flag_r=mean_r[i] < low_r, flag_mahalanobis=md[i] > high_md))
    return pd.DataFrame(rows).set_index("library")


def main() -> None:
    counts = load_counts(ercc=True)
    c = cpm(counts)
    wc = whole_control_mean_r(c)
    others = wc.drop(index=["7_Whole", "10_Whole"])
    screen = compartment_screen(c)
    screen.round(4).to_csv("results/qc_compartment_screen.tsv", sep="\t")
    flagged = screen[(screen.flag_r | screen.flag_mahalanobis)]
    sizes = counts.sum(axis=0) / 1e6
    rows = [
        ("mean_r_7_Whole", round(wc["7_Whole"], 2)),
        ("mean_r_10_Whole", round(wc["10_Whole"], 2)),
        ("other_whole_control_min_mean_r", round(others.min(), 3)),
        ("other_whole_control_min_mean_r_above_0.81", int(others.min() > 0.81)),
        ("mean_r_12_Soma", round(screen.loc["12_Soma", "mean_r"], 2)),
        ("flagged_libraries", ",".join(flagged.index)),
        ("outlier_10_Control_r", round(screen.loc["10_Control", "mean_r"], 2)),
        ("outlier_4_Distal_Dendrites_r", round(screen.loc["4_Distal_Dendrites", "mean_r"], 2)),
        ("outlier_1_Proximal_Dendrites_r", round(screen.loc["1_Proximal_Dendrites", "mean_r"], 2)),
        ("outlier_12_Proximal_Dendrites_mahalanobis", round(screen.loc["12_Proximal_Dendrites", "mahalanobis"], 2)),
        ("outlier_7_Soma_r", round(screen.loc["7_Soma", "mean_r"], 2)),
        ("exon_counts_median_millions", round(sizes.median(), 1)),
        ("exon_counts_min_millions", round(sizes.min(), 1)),
        ("exon_counts_max_millions", round(sizes.max(), 1)),
        ("exon_counts_q1_millions", round(sizes.quantile(0.25), 1)),
        ("exon_counts_q3_millions", round(sizes.quantile(0.75), 1)),
    ]
    write_results("qc", rows)
    for k, v in rows:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
