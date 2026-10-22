"""Shared inputs for every analysis script."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
COUNTS = DATA / "Patterned_Neuron_rmdup_exon_counts.tsv"
METADATA = DATA / "published" / "SuppTable1_collection_metadata.xlsx"
EXTERNAL = DATA / "external"

COMPARTMENTS = ("Distal_Dendrites", "Proximal_Dendrites", "Soma", "Whole", "Control")
CULTURES = tuple(range(1, 13))
WHOLE_OUTLIERS = ("7_Whole", "10_Whole")


def load_counts(ercc: bool = False) -> pd.DataFrame:
    """Exon counts, genes by libraries. ERCC reference rows are dropped unless asked for."""
    df = pd.read_csv(COUNTS, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if not ercc:
        df = df.loc[~df.index.str.startswith("spikeIn.ERCC")]
    return df


def culture(library: str) -> int:
    return int(library.split("_", 1)[0])


def compartment(library: str) -> str:
    return library.split("_", 1)[1]


def batch(library: str) -> int:
    return 1 if culture(library) <= 8 else 2


def libraries(df: pd.DataFrame, comp: str) -> list[str]:
    cols = [c for c in df.columns if compartment(c) == comp]
    return sorted(cols, key=culture)


def cpm(df: pd.DataFrame) -> pd.DataFrame:
    return df / df.sum(axis=0) * 1e6


def proximal_distal_pairs(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    prox = libraries(df, "Proximal_Dendrites")
    dist = libraries(df, "Distal_Dendrites")
    assert [culture(c) for c in prox] == [culture(c) for c in dist] == list(CULTURES)
    return prox, dist


def whole_control_pairs(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    """The nine matched Whole and Control libraries used by the paired micropattern tests."""
    whole, ctrl = [], []
    for n in CULTURES:
        w, c = f"{n}_Whole", f"{n}_Control"
        if w in df.columns and c in df.columns and w not in WHOLE_OUTLIERS:
            whole.append(w)
            ctrl.append(c)
    return whole, ctrl


def expressed_universe(df: pd.DataFrame) -> pd.Index:
    """The 12,396 genes with mean CPM >= 0.5 in the proximal or distal libraries."""
    prox, dist = proximal_distal_pairs(df)
    c = cpm(df)
    keep = (c[prox].mean(axis=1) >= 0.5) | (c[dist].mean(axis=1) >= 0.5)
    return df.index[keep]


def metadata() -> pd.DataFrame:
    meta = pd.read_excel(METADATA)
    kind = {"DistDend": "Distal_Dendrites", "ProxDend": "Proximal_Dendrites",
            "Soma": "Soma", "Whole": "Whole", "Control": "Control"}
    meta["library"] = meta["Culture/Plate"].astype(int).astype(str) + "_" + meta["Sample type"].map(kind)
    return meta.set_index("library")


def pseudo_f(data: np.ndarray, labels: np.ndarray) -> float:
    """PERMANOVA pseudo-F on Euclidean distances; rows are samples."""
    overall = data.mean(axis=0)
    groups = np.unique(labels)
    between = within = 0.0
    for g in groups:
        x = data[labels == g]
        centre = x.mean(axis=0)
        between += len(x) * np.sum((centre - overall) ** 2)
        within += np.sum((x - centre) ** 2)
    k, n = len(groups), len(labels)
    return (between / (k - 1)) / (within / (n - k))


def write_results(name: str, rows: list[tuple[str, object]]) -> Path:
    RESULTS.mkdir(exist_ok=True)
    path = RESULTS / f"{name}.tsv"
    with path.open("w") as fh:
        fh.write("key\tvalue\n")
        for key, value in rows:
            fh.write(f"{key}\t{value}\n")
    return path
