"""Matched proximal and distal neurite libraries shared by the proximal-versus-distal scripts."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pandas as pd

from common import CULTURES, cpm, culture, expressed_universe, load_counts

# Seeded draws list matched sets in sorted culture-string order (1, 10, 11, 12, 2, ..., 9).
SEED_ORDER = tuple(sorted(CULTURES, key=str))
THRESHOLD = 10


class Neurite(NamedTuple):
    counts: pd.DataFrame
    cpm: pd.DataFrame
    universe: pd.Index
    prox: list[str]
    dist: list[str]

    def pair_cpm(self) -> tuple[np.ndarray, np.ndarray]:
        """Proximal and distal CPM over the expressed universe, genes by matched sets."""
        return (self.cpm.loc[self.universe, self.prox].to_numpy(),
                self.cpm.loc[self.universe, self.dist].to_numpy())

    def directions(self) -> np.ndarray:
        """+1 where distal CPM is higher, -1 where proximal CPM is higher, 0 for ties."""
        p, d = self.pair_cpm()
        return np.sign(d - p).astype(np.int8)

    def candidates(self) -> pd.Index:
        """The distal-higher candidates: distal > proximal in at least 10 of 12 matched sets."""
        return self.universe[(self.directions() == 1).sum(axis=1) >= THRESHOLD]

    def neurite_mean(self) -> pd.Series:
        """Mean CPM over all 24 neurite libraries, proximal then distal."""
        cols = sorted(self.prox, key=culture) + sorted(self.dist, key=culture)
        return pd.Series(self.cpm[cols].to_numpy().mean(axis=1), index=self.cpm.index)


def load(order: tuple[int, ...] = CULTURES) -> Neurite:
    counts = load_counts()
    return Neurite(counts, cpm(counts), expressed_universe(counts),
                   [f"{n}_Proximal_Dendrites" for n in order],
                   [f"{n}_Distal_Dendrites" for n in order])
