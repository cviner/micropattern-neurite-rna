#!/usr/bin/env python3
"""Compare every number in claims.tsv with the value the scripts wrote to results/."""

import sys

import pandas as pd

from common import RESULTS, ROOT


def main() -> int:
    claims = pd.read_csv(ROOT / "claims.tsv", sep="\t", dtype=str)
    cache: dict[str, dict[str, str]] = {}
    rows = []
    for c in claims.itertuples(index=False):
        if c.result_file not in cache:
            path = RESULTS / f"{c.result_file}.tsv"
            cache[c.result_file] = (
                dict(pd.read_csv(path, sep="\t", dtype=str).values) if path.exists() else {}
            )
        got = cache[c.result_file].get(c.key)
        note = c.note if isinstance(c.note, str) else ""
        if got is None:
            status = "missing"
        elif abs(float(got) - float(c.expected)) <= float(c.tolerance) + 1e-12:
            status = "match"
        else:
            status = "documented" if note else "differs"
        rows.append((c.claim, c.expected, got, status, note))
    out = pd.DataFrame(rows, columns=["claim", "paper", "computed", "status", "note"])
    out.to_csv(RESULTS / "claims_check.tsv", sep="\t", index=False)
    summary = out["status"].value_counts().to_dict()
    print(summary)
    for r in out[out.status != "match"].itertuples(index=False):
        print(f"{r.status}\t{r.claim}\tpaper {r.paper}\tcomputed {r.computed}\t{r.note}")
    return 0 if set(summary) <= {"match", "documented"} else 1


if __name__ == "__main__":
    sys.exit(main())
