#!/usr/bin/env python3
"""Neurite-versus-soma comparisons against Middleton et al. (2019) and the stress-marker order check.

Table S5 overlap   DESeq2 neurite-higher genes (Supplementary File 2) versus Middleton's 387
                   dendrite-upregulated genes, within the 12,396 expressed genes
Table S6 recovery  Middleton's consistent-detection rule (>= 1 raw count in >= 11 of 12 libraries)
                   applied to our proximal and distal libraries, versus their 1,863 genes
Stress order       canonical stress markers whose mean CPM rises distal < proximal < soma, among the
                   markers passing the micropattern detectability filter
"""

from __future__ import annotations

import pandas as pd

from common import DATA, RESULTS, cpm, expressed_universe, libraries, load_counts, whole_control_pairs, write_results
from middleton import consistently_dendritic, dendrite_upregulated, overlap

SUPP2 = DATA / "published" / "SuppFile2(Dendrites_Soma).xlsx"
SUPP3 = DATA / "published" / "SuppFile3(ConsDistal_ConsProximal).xlsx"
MIN_LIBRARIES = 11
STRESS_MARKERS = ("Hspa5", "Ddit3", "Atf4", "Xbp1", "Bax", "Casp3", "Hspa1a", "Hspa1b", "Hspa8",
                  "Hsp90aa1", "Hsp90ab1", "Cdkn1a", "Trp53", "Bcl2", "Bcl2l1", "Bid", "Hmox1",
                  "Sod1", "Sod2", "Cat", "Gpx1", "Nfkb1", "Atm", "Mdm2", "Aifm1")
ORDER = ("Distal_Dendrites", "Proximal_Dendrites", "Soma")


def consistently_detected(counts: pd.DataFrame, libs: list[str]) -> set[str]:
    return set(counts.index[(counts[libs] >= 1).sum(axis=1) >= MIN_LIBRARIES])


def published_symbols(sheet: str) -> set[str]:
    return set(pd.read_excel(SUPP3, sheet_name=sheet)["Gene Symbol"].dropna().astype(str))


def stress_order(counts: pd.DataFrame, c: pd.DataFrame) -> pd.DataFrame:
    whole, ctrl = whole_control_pairs(counts)
    markers = list(STRESS_MARKERS)
    table = pd.DataFrame({comp: c.loc[markers, libraries(counts, comp)].mean(axis=1) for comp in ORDER})
    table["detectable"] = (c.loc[markers, whole].mean(axis=1) >= 0.5) | (c.loc[markers, ctrl].mean(axis=1) >= 0.5)
    table["rising_toward_soma"] = (table[ORDER[0]] < table[ORDER[1]]) & (table[ORDER[1]] < table[ORDER[2]])
    return table


def main() -> None:
    counts = load_counts()
    c = cpm(counts)
    universe = expressed_universe(counts)

    neurite_higher = set(pd.read_excel(SUPP2, sheet_name="Neurite DE")["Gene Symbol"].dropna().astype(str))
    s5 = overlap(neurite_higher, dendrite_upregulated(), universe)

    prox = consistently_detected(counts, libraries(counts, "Proximal_Dendrites"))
    dist = consistently_detected(counts, libraries(counts, "Distal_Dendrites"))
    for name, ours in (("ConsProximal", prox), ("ConsDistal", dist), ("ConsDistal Specific", dist - prox)):
        assert ours == published_symbols(name), f"{name} differs from Supplementary File 3"
    union = prox | dist
    s6 = consistently_dendritic()
    recovered = union & s6

    stress = stress_order(counts, c)
    stress.to_csv(RESULTS / "stress_marker_compartment_order.tsv", sep="\t", index_label="gene")
    for name, genes in (("consistent_proximal", prox), ("consistent_distal", dist),
                        ("consistent_distal_only", dist - prox), ("middleton_s6_recovered", recovered),
                        ("neurite_higher_x_middleton_s5", s5["genes"])):
        pd.Series(sorted(genes), name="gene").to_csv(RESULTS / f"{name}.tsv", sep="\t", index=False)

    detectable = stress[stress.detectable]
    rows = [
        ("neurite_higher_genes", len(neurite_higher)),
        ("neurite_higher_in_universe", s5["query"]),
        ("middleton_s5_genes", len(dendrite_upregulated())),
        ("middleton_s5_overlap", s5["overlap"]),
        ("middleton_s5_testable", s5["reference"]),
        ("middleton_s5_pct_of_ours", s5["pct_of_query"]),
        ("middleton_s5_pct_of_testable", s5["pct_of_reference"]),
        ("middleton_s5_expected", s5["expected"]),
        ("middleton_s5_fold", s5["fold"]),
        ("middleton_s5_hypergeometric_P", s5["P"]),
        ("consistent_proximal", len(prox)),
        ("consistent_distal", len(dist)),
        ("consistent_both", len(prox & dist)),
        ("consistent_union", len(union)),
        ("consistent_distal_only", len(dist - prox)),
        ("middleton_s6_recovered", len(recovered)),
        ("middleton_s6_total", len(s6)),
        ("middleton_s6_recovered_pct", 100 * len(recovered) / len(s6)),
        ("stress_markers_detectable", int(stress.detectable.sum())),
        ("stress_markers_not_detectable", ",".join(stress.index[~stress.detectable])),
        ("stress_markers_rising_toward_soma", int(detectable.rising_toward_soma.sum())),
        ("stress_markers_rising_toward_soma_all25", int(stress.rising_toward_soma.sum())),
    ]
    print(write_results("neurite_soma", rows))
    for k, v in rows:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
