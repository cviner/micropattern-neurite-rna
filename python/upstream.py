#!/usr/bin/env python3
"""Read-level deduplication: share removed, and per-gene agreement before and after."""

import pandas as pd
from scipy import stats

from common import DATA, load_counts, write_results


def main() -> None:
    dedup = pd.read_csv(DATA / "upstream" / "deduplication_stats.tsv", sep="\t")
    removed = 100 * (1 - dedup.unique_read_pairs / dedup.initial_read_pairs)
    raw = pd.read_csv(DATA / "upstream" / "nondeduplicated_exon_counts_5_libraries.tsv", sep="\t", index_col=0)
    counts = load_counts(ercc=True)
    pearson, spearman, lost = [], [], []
    for lib in raw.columns:
        a, b = raw[lib], counts.loc[raw.index, lib]
        pearson.append(stats.pearsonr(a, b)[0])
        spearman.append(stats.spearmanr(a, b)[0])
        lost.append(int(((a > 0) & (b == 0)).sum()))
    rows = [
        ("libraries", len(dedup)),
        ("dedup_removed_median_pct", round(removed.median(), 0)),
        ("dedup_removed_min_pct", round(removed.min(), 0)),
        ("dedup_removed_max_pct", round(removed.max(), 0)),
        ("dedup_pearson_min", round(min(pearson), 3)),
        ("dedup_pearson_max", round(max(pearson), 3)),
        ("dedup_spearman_min_above_0.999", int(min(spearman) > 0.999)),
        ("dedup_genes_lost_max", max(lost)),
    ]
    write_results("upstream", rows)
    for k, v in rows:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
