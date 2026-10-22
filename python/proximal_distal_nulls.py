#!/usr/bin/env python3
"""Calibration nulls for the size of the distal-higher candidate set (10 of 12 matched sets).

global rate     each gene's non-tied directions redrawn at the single pooled distal-higher rate
abundance       the same, at each neurite-abundance decile's own rate
per-set rate    within each matched set, observed outcomes permuted among genes inside strata of
                abundance decile x distal detection; the null reported in the Results

Seeds: 0 for the global-rate null, 1 for the abundance null, and 2026 for the per-set-rate null.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import RESULTS, write_results
from neurite import SEED_ORDER, THRESHOLD, load

SEED_GLOBAL, SEED_ABUNDANCE, SEED_SET_RATE = 0, 1, 2026
N_GLOBAL = N_ABUNDANCE = 5_000
N_SET_RATE = 20_000
DETECTION_EDGES = (1, 5, 9, 12)


def binomial_null(nontied: np.ndarray, rate: float, rng: np.random.Generator, n: int) -> np.ndarray:
    return np.array([(rng.binomial(nontied, rate) >= THRESHOLD).sum() for _ in range(n)])


def set_rate_null(obs: np.ndarray, strata: np.ndarray, rng: np.random.Generator, n: int) -> np.ndarray:
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    by_set = [[obs[g, j] for g in groups] for j in range(obs.shape[1])]
    out = np.empty(n, dtype=np.int64)
    for b in range(n):
        hits = np.zeros(len(obs), dtype=np.int64)
        for values in by_set:
            hits += np.concatenate([rng.permutation(v) for v in values]) == 1
        out[b] = (hits >= THRESHOLD).sum()
    return out


def summary(null: np.ndarray, observed: int) -> dict:
    exceed = int((null >= observed).sum())
    return {"mean": null.mean(), "sd": null.std(), "max": int(null.max()),
            "exceedances": exceed, "P": (exceed + 1) / (len(null) + 1)}


def main() -> None:
    nd = load(SEED_ORDER)
    obs = nd.directions()
    observed = int(((obs == 1).sum(axis=1) >= THRESHOLD).sum())
    nontied = (obs != 0).sum(axis=1)
    distal_higher = obs == 1
    rate = distal_higher.sum() / (obs != 0).sum()
    per_set = distal_higher.sum(axis=0) / (obs != 0).sum(axis=0)

    global_null = binomial_null(nontied, rate, np.random.default_rng(SEED_GLOBAL), N_GLOBAL)

    abundance = nd.neurite_mean().loc[nd.universe].to_numpy()
    deciles = np.array_split(np.argsort(abundance), 10)
    decile_rates = [distal_higher[d].sum() / (obs[d] != 0).sum() for d in deciles]
    rng = np.random.default_rng(SEED_ABUNDANCE)
    abundance_null = np.zeros(N_ABUNDANCE, dtype=np.int64)
    for d, r in zip(deciles, decile_rates):
        abundance_null += binomial_null(nontied[d], r, rng, N_ABUNDANCE)

    decile = pd.qcut(abundance, 10, labels=False, duplicates="drop").astype(int)
    detected = (nd.cpm.loc[nd.universe, nd.dist].to_numpy() > 0).sum(axis=1)
    strata = decile * 10 + np.digitize(detected, DETECTION_EDGES)
    set_null = set_rate_null(obs, strata, np.random.default_rng(SEED_SET_RATE), N_SET_RATE)

    out = RESULTS / "nulls"
    out.mkdir(parents=True, exist_ok=True)
    for name, null in (("global_rate", global_null), ("abundance_decile", abundance_null),
                       ("per_set_rate", set_null)):
        pd.Series(null, name="distal_higher_ge10of12").to_csv(out / f"{name}_null.tsv", sep="\t", index_label="draw")
    pd.DataFrame({"culture": SEED_ORDER, "distal_higher_rate": per_set}).sort_values("culture").to_csv(
        out / "per_set_distal_higher_rate.tsv", sep="\t", index=False)
    pd.DataFrame({"decile": range(1, 11), "distal_higher_rate": decile_rates}).to_csv(
        out / "abundance_decile_rates.tsv", sep="\t", index=False)

    g, a, s = summary(global_null, observed), summary(abundance_null, observed), summary(set_null, observed)
    rows = [
        ("observed_distal_higher", observed),
        ("global_rate", rate),
        ("global_rate_null_mean", g["mean"]),
        ("global_rate_null_max", g["max"]),
        ("global_rate_null_exceedances", g["exceedances"]),
        ("global_rate_null_P", g["P"]),
        ("global_rate_null_P_below_0.0002", int(g["P"] < 0.0002)),
        ("abundance_null_mean", a["mean"]),
        ("abundance_null_exceedances", a["exceedances"]),
        ("per_set_rate_min", per_set.min()),
        ("per_set_rate_max", per_set.max()),
        ("set_rate_null_mean", s["mean"]),
        ("set_rate_null_sd", s["sd"]),
        ("set_rate_null_max", s["max"]),
        ("set_rate_null_exceedances", s["exceedances"]),
        ("set_rate_null_P", s["P"]),
    ]
    print(write_results("proximal_distal_nulls", rows))
    for k, v in rows:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
