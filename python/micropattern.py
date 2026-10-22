#!/usr/bin/env python3
"""Micropatterned (Whole) versus standard (Control) whole-neuron profiles."""

from __future__ import annotations

import itertools
import math

import numpy as np
import pandas as pd
from scipy import stats

from common import (CULTURES, WHOLE_OUTLIERS, batch, cpm, culture, libraries, load_counts,
                    metadata, pseudo_f, whole_control_pairs, write_results)

STRESS_25 = [
    "Hspa5", "Ddit3", "Atf4", "Xbp1", "Bax", "Casp3", "Hspa1a", "Hspa1b", "Hspa8",
    "Hsp90aa1", "Hsp90ab1", "Cdkn1a", "Trp53", "Bcl2", "Bcl2l1", "Bid", "Hmox1",
    "Sod1", "Sod2", "Cat", "Gpx1", "Nfkb1", "Atm", "Mdm2", "Aifm1",
]


def replay_shuffles(lengths: list[tuple[int, int]]) -> None:
    """Advance numpy's global RNG by the given numbers of shuffles."""
    for n_draws, length in lengths:
        dummy = np.arange(length)
        for _ in range(n_draws):
            np.random.shuffle(dummy)


def nearest_centroid_loo(data: np.ndarray, labels: np.ndarray) -> float:
    correct = 0
    for i in range(len(labels)):
        train = np.arange(len(labels)) != i
        centres = {g: data[train][labels[train] == g].mean(axis=0) for g in np.unique(labels[train])}
        pred = min(centres, key=lambda g: np.sum((data[i] - centres[g]) ** 2))
        correct += int(pred == labels[i])
    return correct / len(labels)


def multivariate_tests(counts: pd.DataFrame) -> list[tuple[str, object]]:
    whole = [c for c in libraries(counts, "Whole") if c not in WHOLE_OUTLIERS]
    ctrl = libraries(counts, "Control")
    sub = counts[whole + ctrl]
    c = cpm(sub)
    keep = c.max(axis=1) >= 1
    x = np.log2(c.loc[keep].to_numpy() + 1).T
    labels = np.array([0] * len(whole) + [1] * len(ctrl))
    f_obs = pseudo_f(x, labels)
    r2 = f_obs / (f_obs + len(labels) - 2)
    acc = nearest_centroid_loo(x, labels)
    chance = max(len(whole), len(ctrl)) / len(labels)

    # Seed 42. The classifier null follows the PERMANOVA shuffles and 2,000 shuffles
    # within each batch, of 24 labels, and of 36 labels.
    np.random.seed(42)
    null_f = []
    for _ in range(2000):
        perm = labels.copy()
        np.random.shuffle(perm)
        null_f.append(pseudo_f(x, perm))
    p_f = (np.sum(np.array(null_f) >= f_obs) + 1) / 2001
    b1 = sum(batch(s) == 1 for s in whole + ctrl)
    b2 = len(labels) - b1
    replay_shuffles([(2000, b1), (2000, b2), (2000, 24), (2000, 36)])
    null_acc = []
    for _ in range(2000):
        perm = labels.copy()
        np.random.shuffle(perm)
        null_acc.append(nearest_centroid_loo(x, perm))
    p_acc = (np.sum(np.array(null_acc) >= acc) + 1) / 2001

    h_dir = abs(2 * math.asin(math.sqrt(acc)) - 2 * math.asin(math.sqrt(chance)))
    return [
        ("permanova_genes", int(keep.sum())),
        ("permanova_pseudo_F", round(f_obs, 4)),
        ("permanova_R2", round(r2, 4)),
        ("permanova_P", round(p_f, 4)),
        ("classifier_accuracy", round(acc, 4)),
        ("classifier_chance", round(chance, 4)),
        ("classifier_P", round(p_acc, 4)),
        ("cohen_h_classifier", round(h_dir, 2)),
    ]


def hotelling(counts: pd.DataFrame) -> list[tuple[str, object]]:
    whole, ctrl = whole_control_pairs(counts)
    c = cpm(counts)
    keep = np.maximum(c[whole].to_numpy(), c[ctrl].to_numpy()).max(axis=1) >= 1
    d = (np.log2(c.loc[keep, whole].to_numpy() + 1) - np.log2(c.loc[keep, ctrl].to_numpy() + 1)).T
    n = d.shape[0]
    t2 = n * np.sum(d.mean(axis=0) ** 2)
    exact = 0
    for signs in itertools.product((-1, 1), repeat=n):
        s = np.array(signs)[:, None]
        exact += int(n * np.sum((d * s).mean(axis=0) ** 2) >= t2 - 1e-9)
    # Seed 42. Consume 5,000 shuffles of 24 and 5,000 sign flips of 12 before the null.
    np.random.seed(42)
    replay_shuffles([(5000, 24)])
    for _ in range(5000):
        np.random.choice([-1, 1], size=12)
    null = []
    for _ in range(5000):
        s = np.random.choice([-1, 1], size=n)[:, None]
        null.append(n * np.sum((d * s).mean(axis=0) ** 2))
    p_mc = (np.sum(np.array(null) >= t2) + 1) / 5001
    return [
        ("hotelling_genes", int(keep.sum())),
        ("hotelling_P", round(p_mc, 4)),
        ("hotelling_P_exact_512", round(exact / 2 ** n, 4)),
    ]


def paired_consistency(counts: pd.DataFrame) -> tuple[list[tuple[str, object]], np.ndarray]:
    whole, ctrl = whole_control_pairs(counts)
    c = cpm(counts)
    expressed = (c[whole].mean(axis=1) >= 0.5) | (c[ctrl].mean(axis=1) >= 0.5)
    diff = (np.log2(c.loc[expressed, whole].to_numpy() + 1)
            - np.log2(c.loc[expressed, ctrl].to_numpy() + 1))
    n = diff.shape[1]

    def consistent(d: np.ndarray) -> tuple[int, int]:
        pos = (d > 0).sum(axis=1) >= 8
        neg = (d < 0).sum(axis=1) >= 8
        return int(pos.sum() + neg.sum()), int(pos.sum())

    total, whole_higher = consistent(diff)
    null = np.array([consistent(diff * np.array(s))[0] for s in itertools.product((-1, 1), repeat=n)])
    p = float((null >= total).sum()) / len(null)
    h = abs(2 * math.asin(math.sqrt(whole_higher / total)) - 2 * math.asin(math.sqrt(0.5)))
    return [
        ("paired_tested_genes", int(expressed.sum())),
        ("consistent_8of9", total),
        ("consistent_whole_higher", whole_higher),
        ("whole_higher_pct", round(100 * whole_higher / total, 1)),
        ("signflip_null_mean", round(float(null.mean()), 1)),
        ("signflip_P", round(p, 4)),
        ("cohen_h_direction", round(h, 2)),
    ], diff


def bootstrap_direction(counts_ercc: pd.DataFrame) -> list[tuple[str, object]]:
    """Paired-culture bootstrap of the whole-higher fraction."""
    whole, ctrl = whole_control_pairs(counts_ercc)
    c = cpm(counts_ercc)
    expressed = (c[whole].mean(axis=1) >= 0.5) | (c[ctrl].mean(axis=1) >= 0.5)
    d1 = (np.log2(c.loc[expressed, whole].to_numpy() + 1)
          - np.log2(c.loc[expressed, ctrl].to_numpy() + 1))
    n1 = d1.shape[1]
    n2 = len(libraries(counts_ercc, "Proximal_Dendrites"))
    rng = np.random.default_rng(42)
    fractions = []
    for _ in range(5000):
        idx = rng.integers(0, n1, size=n1)
        sub = d1[:, idx]
        pos = ((sub > 0).sum(axis=1) >= 8).sum()
        neg = ((sub < 0).sum(axis=1) >= 8).sum()
        fractions.append(pos / max(pos + neg, 1))
        # One unused draw of the 12 neurite pairs per replicate keeps the generator in step.
        rng.integers(0, n2, size=n2)
    return [
        ("bootstrap_ci_low", round(float(np.percentile(fractions, 2.5)), 2)),
        ("bootstrap_ci_high", round(float(np.percentile(fractions, 97.5)), 2)),
    ]


def subsampled_direction(counts_ercc: pd.DataFrame) -> list[tuple[str, object]]:
    """Equalise each pair's library size by multinomial subsampling, then recount direction."""
    whole, ctrl = whole_control_pairs(counts_ercc)
    rng = np.random.default_rng(42)
    c = cpm(counts_ercc)
    expressed = (c[whole].mean(axis=1) >= 0.5) | (c[ctrl].mean(axis=1) >= 0.5)
    # Advance the generator by 3,000 sign flips of the paired libraries before subsampling.
    for _ in range(3000):
        rng.choice([-1, 1], size=len(whole))
    sub = {}
    for w, k in zip(whole, ctrl):
        target = int(min(counts_ercc[w].sum(), counts_ercc[k].sum()))
        for col in (w, k):
            v = counts_ercc[col]
            if v.sum() <= target:
                sub[col] = v.to_numpy()
            else:
                sub[col] = rng.multinomial(target, (v / v.sum()).to_numpy())
    sub = pd.DataFrame(sub, index=counts_ercc.index)
    s = cpm(sub).loc[expressed]
    d = s[whole].to_numpy() - s[ctrl].to_numpy()
    pos = ((d > 0).sum(axis=1) >= 8).sum()
    neg = ((d < 0).sum(axis=1) >= 8).sum()
    return [("subsampled_whole_higher_pct", round(100 * pos / (pos + neg), 1))]


def stress_markers(counts_ercc: pd.DataFrame) -> list[tuple[str, object]]:
    whole, ctrl = whole_control_pairs(counts_ercc)
    c = cpm(counts_ercc)
    expressed = (c[whole].mean(axis=1) >= 0.5) | (c[ctrl].mean(axis=1) >= 0.5)
    lfc = pd.Series(
        (np.log2(c.loc[expressed, whole].to_numpy() + 1)
         - np.log2(c.loc[expressed, ctrl].to_numpy() + 1)).mean(axis=1),
        index=c.index[expressed],
    )
    found = [g for g in STRESS_25 if g in lfc.index]
    markers = lfc.loc[found]
    background = lfc.drop(index=found)
    _, p_w = stats.ranksums(markers.to_numpy(), background.to_numpy())
    rng = np.random.default_rng(42)
    med = float(markers.median())
    null = np.array([np.median(rng.choice(background.to_numpy(), size=len(found), replace=False))
                     for _ in range(5000)])
    return [
        ("stress_markers_detectable", len(found)),
        ("stress_median_log2", round(med, 2)),
        ("stress_wilcoxon_P", round(float(p_w), 3)),
        ("stress_permutation_P", round(float((np.abs(null) >= abs(med)).mean()), 4)),
    ]


def pca_covariates(counts_ercc: pd.DataFrame) -> list[tuple[str, object]]:
    whole, ctrl = whole_control_pairs(counts_ercc)
    cols = [x for pair in zip(whole, ctrl) for x in pair]
    sub = counts_ercc[cols]
    c = cpm(sub)
    pair_max = np.vstack([np.maximum(c[w].to_numpy(), c[k].to_numpy()) for w, k in zip(whole, ctrl)]).max(axis=0)
    x = np.log1p(c.to_numpy()[pair_max >= 0.5])
    x = x - x.mean(axis=1, keepdims=True)
    u, s, _ = np.linalg.svd(x.T, full_matrices=False)
    pcs = u * s
    var = s ** 2 / np.sum(s ** 2) * 100
    meta = metadata()
    lib = sub.sum(axis=0).to_numpy(dtype=float)
    neurons = meta.loc[cols, "Number of neurons collected"].astype(float).to_numpy()
    plating = meta.loc[cols, "Plating concentration (cells/mL)"].astype(float).to_numpy()
    batches = np.array([batch(c_) for c_ in cols])
    plating_p = [stats.pearsonr(pcs[:, k], plating)[1] for k in range(6)]
    r1, p1 = stats.pearsonr(pcs[:, 0], lib)
    r4, p4 = stats.pearsonr(pcs[:, 3], neurons)
    _, p3 = stats.mannwhitneyu(pcs[batches == 1, 2], pcs[batches == 2, 2], alternative="two-sided")
    return [
        ("pca_genes", int((pair_max >= 0.5).sum())),
        ("pc1_variance_pct", round(var[0], 2)),
        ("pc1_library_size_r", round(r1, 2)),
        ("pc1_library_size_P", round(p1, 3)),
        ("pc3_variance_pct", round(var[2], 1)),
        ("pc3_batch_P", round(p3, 3)),
        ("pc4_variance_pct", round(var[3], 1)),
        ("pc4_neurons_r", round(r4, 4)),
        ("pc4_neurons_P", round(p4, 3)),
        ("plating_min_P_top6", round(min(plating_p), 3)),
        ("plating_min_P_above_0.10", int(min(plating_p) > 0.10)),
    ]


def control_batch_genes(counts: pd.DataFrame) -> list[tuple[str, object]]:
    ctrl = libraries(counts, "Control")
    c = cpm(counts[ctrl])
    c = c.loc[c.mean(axis=1) >= 0.5]
    b1 = [x for x in ctrl if batch(x) == 1]
    b2 = [x for x in ctrl if batch(x) == 2]
    p = np.array([stats.mannwhitneyu(r[b1], r[b2], alternative="two-sided").pvalue for _, r in c.iterrows()])
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    return [("control_batch_tested", len(p)), ("control_batch_de_genes", int((q < 0.05).sum()))]


def main() -> None:
    counts = load_counts()
    counts_ercc = load_counts(ercc=True)
    rows = []
    rows += multivariate_tests(counts)
    rows += hotelling(counts)
    paired, _ = paired_consistency(counts)
    rows += paired
    rows += bootstrap_direction(counts_ercc)
    rows += subsampled_direction(counts_ercc)
    rows += stress_markers(counts_ercc)
    rows += pca_covariates(counts_ercc)
    rows += control_batch_genes(counts)
    path = write_results("micropattern", rows)
    for k, v in rows:
        print(f"{k}\t{v}")
    print(path)


if __name__ == "__main__":
    main()
