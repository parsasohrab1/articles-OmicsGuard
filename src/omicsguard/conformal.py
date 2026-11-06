from __future__ import annotations

import numpy as np


def _q(scores: np.ndarray, alpha: float) -> float:
    n = len(scores)
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return float(np.sort(scores)[k - 1]) if k <= n else float("inf")


def classwise_thresholds(P: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    th = np.full(P.shape[1], np.inf)
    for c in range(P.shape[1]):
        m = y == c
        if m.sum():
            th[c] = _q(1 - P[m, c], alpha)
    return th


def prediction_sets(P: np.ndarray, th: np.ndarray) -> np.ndarray:
    return (1 - P) <= th


def refer(S: np.ndarray) -> np.ndarray:
    """True when the set is not a singleton (empty or several labels)."""
    return S.sum(axis=1) != 1


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """P(score_pos > score_neg) with ties counted half (Mann-Whitney)."""
    allv = np.concatenate([pos, neg])
    ranks = allv.argsort().argsort().astype(float) + 1
    # average ranks for ties
    order = np.argsort(allv, kind="stable")
    sv = allv[order]
    i = 0
    r = np.empty(len(allv))
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        r[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    del ranks
    u = r[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2
    return float(u / (len(pos) * len(neg)))
