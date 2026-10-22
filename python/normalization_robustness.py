#!/usr/bin/env python3
"""Proximal-versus-distal direction under alternative library normalizations.

Each normalization is fitted on all 59 libraries; the gene universe stays the CPM-defined 12,396 so
only the normalization varies. A gene is proximal-higher if it is higher proximally in more matched
sets than distally; equal counts are tied and belong to neither direction.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import RESULTS, write_results
from neurite import THRESHOLD, load


def cpm(m: np.ndarray) -> np.ndarray:
    return m / m.sum(axis=0) * 1e6


def median_of_ratios(m: np.ndarray) -> np.ndarray:
    logs = np.log(m[(m > 0).all(axis=1)])
    return m / np.exp(np.median(logs - logs.mean(axis=1, keepdims=True), axis=0))


def upper_quartile(m: np.ndarray) -> np.ndarray:
    sf = np.array([np.percentile(col[col > 0], 75) for col in m.T])
    return m / (sf / np.exp(np.log(sf).mean()))


def tmm(m: np.ndarray, logratio_trim: float = 0.3, sum_trim: float = 0.05) -> np.ndarray:
    """edgeR TMM with the reference library whose upper quartile is closest to the mean."""
    lib = m.sum(axis=0)
    f75 = np.array([np.percentile(col[col > 0] / n, 75) for col, n in zip(m.T, lib)])
    ref = int(np.argmin(np.abs(f75 - f75.mean())))
    sf = np.ones(m.shape[1])
    for j in range(m.shape[1]):
        if j == ref:
            continue
        o, r = m[:, j], m[:, ref]
        keep = (o > 0) & (r > 0)
        oj, rj = o[keep] / lib[j], r[keep] / lib[ref]
        mv, av = np.log2(oj / rj), 0.5 * np.log2(oj * rj)
        var = (1 - oj) / o[keep] + (1 - rj) / r[keep]
        ok = np.isfinite(mv) & np.isfinite(av) & np.isfinite(var) & (var > 0)
        mv, av, var = mv[ok], av[ok], var[ok]
        if len(mv) < 10:
            continue
        lo_m, hi_m = np.quantile(mv, [logratio_trim, 1 - logratio_trim])
        lo_a, hi_a = np.quantile(av, [sum_trim, 1 - sum_trim])
        k = (mv >= lo_m) & (mv <= hi_m) & (av >= lo_a) & (av <= hi_a)
        if k.any():
            sf[j] = 2 ** (np.sum(mv[k] / var[k]) / np.sum(1 / var[k]))
    sf = sf * lib
    return m / (sf / np.exp(np.log(sf).mean()))


def ranks(m: np.ndarray) -> np.ndarray:
    return pd.DataFrame(m).rank(method="average").to_numpy()


NORMALIZATIONS = {"cpm": cpm, "median_of_ratios": median_of_ratios, "tmm": tmm,
                  "upper_quartile": upper_quartile, "rank": ranks}


def main() -> None:
    nd = load()
    raw = nd.counts.to_numpy(dtype=np.float64)
    rows_u = np.flatnonzero(nd.counts.index.isin(nd.universe))
    pi = [nd.counts.columns.get_loc(c) for c in nd.prox]
    di = [nd.counts.columns.get_loc(c) for c in nd.dist]
    candidates = set(nd.candidates())

    table, members = [], {}
    for name, fn in NORMALIZATIONS.items():
        norm = fn(raw)
        p, d = norm[np.ix_(rows_u, pi)], norm[np.ix_(rows_u, di)]
        n_p, n_d = (p > d).sum(axis=1), (d > p).sum(axis=1)
        members[name] = set(nd.counts.index[rows_u[n_d >= THRESHOLD]])
        table.append({
            "normalization": name,
            "proximal_higher_pct": 100 * (n_p > n_d).mean(),
            "distal_higher_pct": 100 * (n_d > n_p).mean(),
            "tied_pct": 100 * (n_p == n_d).mean(),
            "distal_higher_ge10of12": len(members[name]),
            "overlap_with_54": len(members[name] & candidates),
        })
    table = pd.DataFrame(table)
    table.to_csv(RESULTS / "normalization_robustness_table.tsv", sep="\t", index=False)
    tiers = pd.DataFrame([(g, sum(g in s for s in members.values()),
                           ",".join(n for n, s in members.items() if g in s)) for g in sorted(candidates)],
                         columns=["gene", "n_normalizations", "retained_in"])
    tiers.to_csv(RESULTS / "candidates_by_normalization.tsv", sep="\t", index=False)

    t = table.set_index("normalization")
    scaling = t.drop(index="rank")
    rows = [(f"proximal_higher_pct_{n}", t.loc[n, "proximal_higher_pct"]) for n in t.index]
    rows += [("tied_pct_scaling_min", scaling["tied_pct"].min()),
             ("tied_pct_scaling_max", scaling["tied_pct"].max()),
             ("candidates_retained_all_normalizations", int((tiers.n_normalizations == len(members)).sum()))]
    print(write_results("normalization_robustness", rows))
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
