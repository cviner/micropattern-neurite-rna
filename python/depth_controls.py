#!/usr/bin/env python3
"""Read-depth controls for the 54 distal-higher candidates.

Depth matching   thin the deeper library of each matched set to the pair's smaller depth, then count
                 how many of the 54 stay distal-higher in at least 10 of 12 sets.
Pseudo-distal    replace each distal library by a thinned copy of its proximal partner, so that no
                 gene is truly distal-higher, then count calls and their overlap with the 54:
                   paired      thinned to the partner's distal depth (changes only the pairs in which
                               proximal is deeper, so 10 of 12 is unreachable by construction)
                   global_min  thinned to the shallowest distal depth
                   half        thinned to half of its own depth

Thinning models: binomial per-gene, exact hypergeometric draw of reads, multinomial with replacement.
Seeds are default_rng(200 + seed) for depth matching and default_rng(100 + seed) for each
pseudo-distal variant. Matched sets are in sorted culture-string order. pseudo_distal_max_overlap
is the maximum overlap_with_54 over seeds 0-4 of the global_min and half variants under binomial,
hypergeometric, and multinomial thinning. The paired variant is excluded.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import RESULTS, write_results
from neurite import SEED_ORDER, THRESHOLD, load

MODELS = ("binomial", "hypergeometric", "multinomial")
SEEDS = range(5)


def thin(counts: np.ndarray, target: int, rng: np.random.Generator, model: str) -> np.ndarray:
    depth = int(counts.sum())
    if target >= depth:
        return counts
    if model == "binomial":
        return rng.binomial(counts, target / depth)
    if model == "hypergeometric":
        return rng.multivariate_hypergeometric(counts, int(target))
    if model == "multinomial":
        return rng.multinomial(int(target), counts / depth)
    raise ValueError(model)


def cpm(m: np.ndarray) -> np.ndarray:
    return m / m.sum(axis=0) * 1e6


def depth_matched(prox: np.ndarray, dist: np.ndarray, rows: np.ndarray, seed: int, model: str) -> int:
    rng = np.random.default_rng(200 + seed)
    p, d = prox.copy(), dist.copy()
    for j in range(p.shape[1]):
        target = min(p[:, j].sum(), d[:, j].sum())
        p[:, j] = thin(prox[:, j], target, rng, model)
        d[:, j] = thin(dist[:, j], target, rng, model)
    higher = (cpm(d)[rows] > cpm(p)[rows]).sum(axis=1)
    return int((higher >= THRESHOLD).sum())


def pseudo_distal(prox: np.ndarray, targets: np.ndarray, universe: np.ndarray, rows54: np.ndarray,
                  seed: int, model: str) -> tuple[int, int]:
    rng = np.random.default_rng(100 + seed)
    fake = np.column_stack([thin(prox[:, j], int(targets[j]), rng, model) for j in range(prox.shape[1])])
    higher = cpm(fake) > cpm(prox)
    calls = higher[universe].sum(axis=1) >= THRESHOLD
    return int(calls.sum()), int((higher[rows54].sum(axis=1) >= THRESHOLD).sum())


def main() -> None:
    nd = load(SEED_ORDER)
    genes = nd.counts.index
    prox = nd.counts[nd.prox].to_numpy(dtype=np.int64)
    dist = nd.counts[nd.dist].to_numpy(dtype=np.int64)
    universe = np.flatnonzero(genes.isin(nd.universe))
    rows54 = np.flatnonzero(genes.isin(nd.candidates()))
    assert len(universe) == 12_396 and len(rows54) == 54

    matched = pd.DataFrame([(m, s, depth_matched(prox, dist, rows54, s, m)) for m in MODELS for s in SEEDS],
                           columns=["model", "seed", "retained_of_54"])

    pdepth, ddepth = prox.sum(axis=0), dist.sum(axis=0)
    targets = {"paired": ddepth, "global_min": np.full(12, ddepth.min()), "half": pdepth // 2}
    pseudo = pd.DataFrame(
        [(v, m, s, *pseudo_distal(prox, t, universe, rows54, s, m))
         for m in MODELS for v, t in targets.items() for s in SEEDS],
        columns=["variant", "model", "seed", "calls", "overlap_with_54"])

    matched.to_csv(RESULTS / "depth_matching_runs.tsv", sep="\t", index=False)
    pseudo.to_csv(RESULTS / "pseudo_distal_runs.tsv", sep="\t", index=False)

    retained = matched["retained_of_54"]
    thinned_all = pseudo[pseudo.variant != "paired"]
    rows = [
        ("pairs_distal_deeper", int((ddepth > pdepth).sum())),
        ("depth_matched_retained_min", int(retained.min())),
        ("depth_matched_retained_max", int(retained.max())),
        ("depth_matched_retained_median", float(retained.median())),
        *[(f"depth_matched_retained_{m}", ",".join(map(str, g["retained_of_54"])))
          for m, g in matched.groupby("model", sort=False)],
        ("pseudo_distal_paired_max_calls", int(pseudo.loc[pseudo.variant == "paired", "calls"].max())),
        ("pseudo_distal_max_overlap", int(thinned_all["overlap_with_54"].max())),
        *[(f"pseudo_distal_max_overlap_{m}", int(g["overlap_with_54"].max()))
          for m, g in thinned_all.groupby("model", sort=False)],
        *[(f"pseudo_distal_{v}_{m}_calls", ",".join(map(str, g["calls"])))
          for (v, m), g in pseudo.groupby(["variant", "model"], sort=False)],
        *[(f"pseudo_distal_{v}_{m}_overlap", ",".join(map(str, g["overlap_with_54"])))
          for (v, m), g in pseudo.groupby(["variant", "model"], sort=False)],
    ]
    print(write_results("depth_controls", rows))
    for k, v in rows:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
