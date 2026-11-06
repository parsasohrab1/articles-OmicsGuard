import numpy as np

from omicsguard.conformal import auroc, classwise_thresholds, prediction_sets, refer


def test_auroc_basic_cases():
    assert auroc(np.array([3.0, 4.0]), np.array([1.0, 2.0])) == 1.0
    assert auroc(np.array([1.0, 2.0]), np.array([3.0, 4.0])) == 0.0
    assert auroc(np.array([1.0, 1.0]), np.array([1.0, 1.0])) == 0.5


def test_classwise_coverage_and_refer():
    rng = np.random.default_rng(0)

    def draw(n):
        y = rng.integers(0, 3, n)
        logit = rng.normal(0, 1, (n, 3))
        logit[np.arange(n), y] += 1.2
        P = np.exp(logit) / np.exp(logit).sum(1, keepdims=True)
        return P, y

    Pc, yc = draw(6000)
    Pt, yt = draw(9000)
    S = prediction_sets(Pt, classwise_thresholds(Pc, yc, 0.1))
    for c in range(3):
        assert abs(S[yt == c, c].mean() - 0.9) < 0.03
    assert list(refer(np.array([[1, 0], [1, 1], [0, 0]], dtype=bool))) == [False, True, True]
