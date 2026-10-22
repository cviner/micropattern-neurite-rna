#!/usr/bin/env python3
"""Proximal-versus-distal neurite comparison over the 12 matched sets.

Compartment separation (exact restricted PERMANOVA and leave-one-pair-out), the 10-of-12 consistency
sets and their exact sign-flip null, the Middleton Table S4 overlap of the 54 candidates, marker
panels, and descriptive checks (abundance, detection, read concentration, dispersion, pairing).

Writes results/proximal_distal.tsv together with the keys from proximal_distal_nulls.py,
depth_controls.py and normalization_robustness.py, which must be run first.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from common import CULTURES, DATA, RESULTS, pseudo_f, write_results
from middleton import overlap, published_catalogue
from neurite import THRESHOLD, Neurite, load

SECTIONS = ("proximal_distal_nulls", "depth_controls", "normalization_robustness")
SUPP4 = DATA / "published" / "SuppFile4(Distal_Proximal).xlsx"
PERMANOVA_MIN_CPM = 1.0
POLARITY = {"Kif5a": 1, "Ppfia3": 1, "Camk2a": -1, "Bdnf": -1, "Gria1": -1, "Map2": -1}
REFERENCE = ("Gapdh", "Actb", "Hprt", "B2m", "Ywhaz", "Tbp", "Sdha", "Pgk1")
AXON = ("Mapt", "Snap25", "Syt1", "Dlg4", "Shank3", "Map2", "Gap43", "Nefl")
GFAP_ELEVATED_CPM = 500
UNRESTRICTED_PERMUTATIONS = 1_000
# Relative tolerance for ties: complementary label sets give the same statistic up to rounding.
TIE = 1e-9


def pseudo_f_many(d2: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Two-group PERMANOVA pseudo-F for every row of a 0/1 label matrix, from squared distances."""
    n = d2.shape[0]
    within = 0.0
    for g in (0, 1):
        u = (labels == g).astype(float)
        within = within + np.einsum("pi,ij,pj->p", u, d2, u) / (2 * u.sum(axis=1))
    return (d2.sum() / (2 * n) - within) / (within / (n - 2))


def pair_swaps(n_pairs: int) -> np.ndarray:
    """All 2**n_pairs within-pair label exchanges for samples ordered proximal 1..n, distal 1..n.
    Row 0 is the observed labelling; label 1 marks the distal group."""
    bits = (np.arange(2 ** n_pairs)[:, None] >> np.arange(n_pairs)) & 1
    return np.hstack([bits, 1 - bits])


def drop_pair(n_pairs: int, k: int) -> np.ndarray:
    keep = [i for i in range(n_pairs) if i != k]
    return np.array(keep + [n_pairs + i for i in keep])


def cauchy_combination(p: np.ndarray) -> float:
    return float(0.5 - np.arctan(np.mean(np.tan((0.5 - p) * np.pi))) / np.pi)


def unrestricted_lopo(d2: np.ndarray) -> np.ndarray:
    """Unrestricted label shuffles per omitted set, P = (b + 1) / 1001.
    RandomState(42) first draws 3,000 permutations of 24 and 1,000 of 21."""
    rs = np.random.RandomState(42)
    for _ in range(3000):
        rs.permutation(24)
    for _ in range(1000):
        rs.permutation(21)
    labels = np.repeat([0, 1], 11)
    p = []
    for k in range(12):
        sub = d2[np.ix_(drop_pair(12, k), drop_pair(12, k))]
        null = pseudo_f_many(sub, np.array([rs.permutation(labels) for _ in range(UNRESTRICTED_PERMUTATIONS)]))
        observed = pseudo_f_many(sub, labels[None])[0]
        p.append(((null >= observed * (1 - TIE)).sum() + 1) / (UNRESTRICTED_PERMUTATIONS + 1))
    return np.array(p)


def separation(nd: Neurite) -> list[tuple[str, object]]:
    cols = nd.prox + nd.dist
    keep = nd.cpm[cols].max(axis=1) >= PERMANOVA_MIN_CPM
    x = np.log2(nd.cpm.loc[keep, cols].to_numpy().T + 1)
    labels = np.repeat([0, 1], 12)
    f = pseudo_f(x, labels)
    d2 = squareform(pdist(x, "sqeuclidean"))
    swaps = pair_swaps(12)
    null = pseudo_f_many(d2, swaps)
    p = float(np.mean(null >= null[0] * (1 - TIE)))

    lopo = []
    swaps11 = pair_swaps(11)
    for k in range(12):
        rows = drop_pair(12, k)
        f_k = pseudo_f(x[rows], np.repeat([0, 1], 11))
        null_k = pseudo_f_many(d2[np.ix_(rows, rows)], swaps11)
        lopo.append((CULTURES[k], f_k, f_k / (f_k + 20), float(np.mean(null_k >= null_k[0] * (1 - TIE)))))
    lopo = pd.DataFrame(lopo, columns=["omitted_culture", "pseudo_F", "R2", "exact_P"])
    lopo["unrestricted_P"] = unrestricted_lopo(d2)
    lopo.to_csv(RESULTS / "proximal_distal_lopo.tsv", sep="\t", index=False)
    omit12 = lopo.set_index("omitted_culture").loc[12]
    return [
        ("permanova_genes", int(keep.sum())),
        ("permanova_pseudo_F", f),
        ("permanova_R2", f / (f + 22)),
        ("permanova_P", p),
        ("lopo_P_min", lopo.exact_P.min()),
        ("lopo_P_max", lopo.exact_P.max()),
        ("lopo_R2_min", lopo.R2.min()),
        ("lopo_R2_max", lopo.R2.max()),
        ("exact_omit12_P", omit12.exact_P),
        ("unrestricted_omit12_P", omit12.unrestricted_P),
        ("unrestricted_cauchy_P", cauchy_combination(lopo.unrestricted_P.to_numpy())),
    ]


def symmetric_null(dirs: np.ndarray) -> pd.DataFrame:
    """Exact within-pair sign-flip null: every one of the 4,096 flips of the 12 matched sets."""
    flips = 1 - 2 * ((np.arange(4096)[:, None] >> np.arange(12)) & 1).astype(np.int8)
    out = np.empty((4096, 2), dtype=np.int64)
    for i, f in enumerate(flips):
        s = dirs * f
        out[i] = ((s == 1).sum(axis=1) >= THRESHOLD).sum(), ((s == -1).sum(axis=1) >= THRESHOLD).sum()
    null = pd.DataFrame(out, columns=["distal_higher", "proximal_higher"])
    null["total"] = null.sum(axis=1)
    return null


def candidate_table(nd: Neurite, dirs: np.ndarray, s4: set[str]) -> pd.DataFrame:
    p, d = nd.pair_cpm()
    table = pd.DataFrame({
        "gene": nd.universe,
        "distal_higher_sets": (dirs == 1).sum(axis=1),
        "mean_distal_cpm": d.mean(axis=1),
        "mean_proximal_cpm": p.mean(axis=1),
        "detected_distal_sets": (nd.counts.loc[nd.universe, nd.dist] >= 1).sum(axis=1).to_numpy(),
    })
    table["log2_mean_ratio"] = np.log2((table.mean_distal_cpm + 1) / (table.mean_proximal_cpm + 1))
    table["middleton_s4"] = table.gene.isin(s4)
    return table[table.distal_higher_sets >= THRESHOLD].reset_index(drop=True)


def check_against_supp4(table: pd.DataFrame) -> None:
    pub = pd.read_excel(SUPP4, sheet_name="54 candidate genes").dropna(subset=["Distal-higher sets (of 12)"])
    pub = pub.set_index("Gene symbol")
    assert set(pub.index) == set(table.gene), "candidate set differs from Supplementary File 4"
    ours = table.set_index("gene").loc[pub.index]
    assert (ours.distal_higher_sets == pub["Distal-higher sets (of 12)"]).all()
    assert (ours.detected_distal_sets == pub["Detected distal sets (of 12)"]).all()
    assert (ours.middleton_s4 == (pub["Middleton et al. 2019 Table S4"] == "Yes")).all()
    assert np.allclose(ours.mean_distal_cpm.round(1), pub["Mean distal CPM"])


def consistency(nd: Neurite) -> list[tuple[str, object]]:
    dirs = nd.directions()
    n_d, n_p = (dirs == 1).sum(axis=1), (dirs == -1).sum(axis=1)
    distal, proximal = nd.universe[n_d >= THRESHOLD], nd.universe[n_p >= THRESHOLD]

    null = symmetric_null(dirs)
    null.to_csv(RESULTS / "nulls" / "symmetric_signflip_null.tsv", sep="\t", index_label="flip")
    observed = int(null.total[0])

    s4, s4_external = published_catalogue(), published_catalogue(external_only=True)
    table = candidate_table(nd, dirs, s4)
    check_against_supp4(table)
    table.sort_values(["distal_higher_sets", "gene"], ascending=[False, True]).to_csv(
        RESULTS / "candidates_distal_54.tsv", sep="\t", index=False)
    pd.Series(sorted(proximal), name="gene").to_csv(RESULTS / "proximal_higher_1154.tsv", sep="\t", index=False)
    cat = overlap(distal, s4, nd.universe)
    ext = overlap(distal, s4_external, nd.universe)

    omit = dirs[:, [i for i, n in enumerate(CULTURES) if n != 12]]
    omit_d, omit_p = (omit == 1).sum(axis=1), (omit == -1).sum(axis=1)
    abundance = nd.neurite_mean()
    return [
        ("expressed_universe", len(nd.universe)),
        ("consistent_10of12_total", len(distal) + len(proximal)),
        ("consistent_10of12_proximal", len(proximal)),
        ("consistent_10of12_distal", len(distal)),
        ("symmetric_null_mean", null.total.mean()),
        ("symmetric_null_distal_mean", null.distal_higher.mean()),
        ("symmetric_null_P", float((null.total >= observed).mean())),
        ("symmetric_null_fold", observed / null.total.mean()),
        ("distal_higher_9of12", int((n_d >= 9).sum())),
        ("distal_higher_11of12", int((n_d >= 11).sum())),
        ("distal_higher_12of12", int((n_d >= 12).sum())),
        ("omit12_consistent_9of11", int(((omit_d >= 9) | (omit_p >= 9)).sum())),
        ("omit12_consistent_10of11", int(((omit_d >= 10) | (omit_p >= 10)).sum())),
        ("middleton_s4_reference_in_universe", cat["reference"]),
        ("middleton_s4_overlap", cat["overlap"]),
        ("middleton_s4_expected", cat["expected"]),
        ("middleton_s4_fold", cat["fold"]),
        ("middleton_s4_P", cat["P"]),
        ("middleton_s4_external_reference_in_universe", ext["reference"]),
        ("middleton_s4_external_overlap", ext["overlap"]),
        ("middleton_s4_external_fold", ext["fold"]),
        ("middleton_s4_external_P", ext["P"]),
        ("candidates_median_cpm", abundance[distal].median()),
        ("universe_median_cpm", abundance[nd.universe].median()),
        ("candidates_detected_11of12_distal", int((table.detected_distal_sets >= 11).sum())),
    ]


def log2_ratios(nd: Neurite, gene: str) -> np.ndarray:
    """Per matched set log2[(distal CPM + 1) / (proximal CPM + 1)]."""
    return np.log2((nd.cpm.loc[gene, nd.dist].to_numpy() + 1) / (nd.cpm.loc[gene, nd.prox].to_numpy() + 1))


def markers(nd: Neurite) -> list[tuple[str, object]]:
    panel = []
    for gene, sign in POLARITY.items():
        d, p = nd.cpm.loc[gene, nd.dist].to_numpy(), nd.cpm.loc[gene, nd.prox].to_numpy()
        panel.append((gene, "distal" if sign > 0 else "proximal", np.median(log2_ratios(nd, gene)),
                      100 * np.mean(np.sign(d - p) == sign)))
    panel += [(gene, "reference", np.median(log2_ratios(nd, gene)), np.nan) for gene in REFERENCE]
    panel = pd.DataFrame(panel, columns=["gene", "expected_higher", "median_log2_distal_over_proximal",
                                         "pct_sets_expected_direction"])
    axon = pd.DataFrame({"gene": AXON, "median_distal_over_proximal": [2 ** np.median(log2_ratios(nd, g))
                                                                       for g in AXON]})
    panel.to_csv(RESULTS / "polarity_markers.tsv", sep="\t", index=False)
    axon.to_csv(RESULTS / "axon_markers.tsv", sep="\t", index=False)
    polar = panel[panel.expected_higher != "reference"]
    ref = panel[panel.expected_higher == "reference"]
    gfap = nd.cpm.loc["Gfap", nd.dist]
    return [
        ("polarity_median_pct_expected", polar.pct_sets_expected_direction.median()),
        ("polarity_markers_median_abs_log2", polar.median_log2_distal_over_proximal.abs().median()),
        ("polarity_reference_median_abs_log2", ref.median_log2_distal_over_proximal.abs().median()),
        *[(f"axon_{g}", v) for g, v in zip(axon.gene, axon.median_distal_over_proximal)],
        ("gfap_elevated_distal", int((gfap > GFAP_ELEVATED_CPM).sum())),
        ("gfap_elevated_distal_cultures", ",".join(c.split("_")[0] for c in gfap.index[gfap > GFAP_ELEVATED_CPM])),
    ]


def top100_share(counts: pd.DataFrame, libs: list[str]) -> float:
    return float(np.mean([np.sort(counts[c].to_numpy())[::-1][:100].sum() / counts[c].sum() for c in libs]))


def dispersion(d2: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Mean Euclidean distance to the group centroid, per labelling (rows) and group (columns)."""
    out = []
    for g in (0, 1):
        u = (labels == g).astype(float)
        m = u.sum(axis=1, keepdims=True)
        spread = np.einsum("pi,ij,pj->p", u, d2, u)[:, None] / (2 * m ** 2)
        dist = np.sqrt(np.clip(u @ d2 / m - spread, 0, None))
        out.append((dist * u).sum(axis=1) / m[:, 0])
    return np.column_stack(out)


def within_minus_across(c: pd.DataFrame, a: str, b: str) -> float:
    r = np.corrcoef(c[[f"{n}_{a}" for n in CULTURES]].T, c[[f"{n}_{b}" for n in CULTURES]].T)[:12, 12:]
    return float(np.diag(r).mean() - r[~np.eye(12, dtype=bool)].mean())


def descriptive(nd: Neurite) -> list[tuple[str, object]]:
    x = np.log2(nd.cpm.loc[nd.universe, nd.prox + nd.dist].to_numpy().T + 1)
    disp = dispersion(squareform(pdist(x, "sqeuclidean")), pair_swaps(12))
    diff = disp[:, 1] - disp[:, 0]
    pd.DataFrame(disp, columns=["proximal", "distal"]).to_csv(
        RESULTS / "nulls" / "permdisp_null.tsv", sep="\t", index_label="swap")
    return [
        ("top100_distal_pct", 100 * top100_share(nd.counts, nd.dist)),
        ("top100_proximal_pct", 100 * top100_share(nd.counts, nd.prox)),
        ("permdisp_proximal", disp[0, 0]),
        ("permdisp_distal", disp[0, 1]),
        ("permdisp_P", float(np.mean(np.abs(diff) >= abs(diff[0]) * (1 - TIE)))),
        ("within_minus_across_r_prox_dist", within_minus_across(nd.cpm, "Proximal_Dendrites", "Distal_Dendrites")),
        ("within_minus_across_r_prox_soma", within_minus_across(nd.cpm, "Proximal_Dendrites", "Soma")),
        ("within_minus_across_r_dist_soma", within_minus_across(nd.cpm, "Distal_Dendrites", "Soma")),
    ]


def section_rows() -> list[tuple[str, object]]:
    rows = []
    for name in SECTIONS:
        path = RESULTS / f"{name}.tsv"
        if not path.exists():
            raise SystemExit(f"{path} is missing; run python/{name}.py first")
        rows += list(pd.read_csv(path, sep="\t", dtype=str).itertuples(index=False, name=None))
    return rows


def main() -> None:
    (RESULTS / "nulls").mkdir(parents=True, exist_ok=True)
    nd = load()
    rows = separation(nd) + consistency(nd) + markers(nd) + descriptive(nd)
    for k, v in rows:
        print(f"{k}\t{v}")
    rows += section_rows()
    keys = [k for k, _ in rows]
    assert len(keys) == len(set(keys)), "duplicate result keys"
    print(write_results("proximal_distal", rows))


if __name__ == "__main__":
    main()
