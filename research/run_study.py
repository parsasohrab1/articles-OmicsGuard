"""TCGA PanCancer RNA-Seq (UCI 401; 801 tumours, 5 types, 20 531 genes): conformal sets and unseen-type referral.

usage: python research/run_study.py <dir with data.csv and labels.csv>
(1) iid: 20 random splits 50/25/25 train/calibration/test; class-wise conformal at alpha 0.05.
(2) open-set: each tumour type in turn is removed from training and calibration and appears only in the test set;
    referral = prediction set is empty or has several labels. Compared with a max-probability rule with threshold at the
    5 % quantile of calibration max-probabilities.
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from omicsguard.conformal import auroc, classwise_thresholds, prediction_sets, refer

warnings.filterwarnings("ignore")
ALPHA = 0.05
REPS = 10
TOP_GENES = 2000


def top_var_idx(X, k):
    v = X.var(axis=0)
    return np.argsort(-v)[:k]


def fit(Xtr, ytr, genes):
    m = make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=2000))
    return m.fit(Xtr[:, genes], ytr)


def main(folder: str):
    X = pd.read_csv(f"{folder}/data.csv", index_col=0).to_numpy(np.float32)
    lab = pd.read_csv(f"{folder}/labels.csv", index_col=0)["Class"]
    classes = sorted(lab.unique())
    y = lab.map({c: i for i, c in enumerate(classes)}).to_numpy()
    out = {
        "note": "real TCGA PanCancer RNA-Seq (UCI 401); top-2000 variance genes chosen on the training part; "
        "logistic regression C=0.05; alpha=0.05",
        "classes": classes,
        "iid": [],
        "open_set": {},
    }
    # (1) iid
    accs, covs, covc, refs = [], [], [], []
    for r in range(20):
        idx = np.arange(len(y))
        tr, rest = train_test_split(idx, test_size=0.5, stratify=y, random_state=r)
        ca, te = train_test_split(rest, test_size=0.5, stratify=y[rest], random_state=r)
        g = top_var_idx(X[tr], TOP_GENES)
        m = fit(X[tr], y[tr], g)
        Pc, Pt = m.predict_proba(X[ca][:, g]), m.predict_proba(X[te][:, g])
        S = prediction_sets(Pt, classwise_thresholds(Pc, y[ca], ALPHA))
        accs.append(float((Pt.argmax(1) == y[te]).mean()))
        covs.append(float(S[np.arange(len(te)), y[te]].mean()))
        covc.append([float(S[y[te] == c, c].mean()) for c in range(len(classes))])
        refs.append(float(refer(S).mean()))
    out["iid"] = {
        "accuracy": float(np.mean(accs)),
        "coverage": float(np.mean(covs)),
        "coverage_by_class": dict(zip(classes, np.mean(covc, axis=0).tolist(), strict=True)),
        "referral_rate": float(np.mean(refs)),
    }
    print(
        "iid",
        {k: v for k, v in out["iid"].items() if k != "coverage_by_class"},
        out["iid"]["coverage_by_class"],
    )
    # (2) open set
    for c_out, cname in enumerate(classes):
        res = {
            "auroc_maxprob": [],
            "unseen_referred_conformal": [],
            "seen_referred_conformal": [],
            "unseen_referred_maxprob": [],
            "seen_referred_maxprob": [],
        }
        known = [c for c in range(len(classes)) if c != c_out]
        remap = {c: i for i, c in enumerate(known)}
        for r in range(REPS):
            idx_k = np.where(y != c_out)[0]
            tr, rest = train_test_split(idx_k, test_size=0.5, stratify=y[idx_k], random_state=r)
            ca, te = train_test_split(rest, test_size=0.5, stratify=y[rest], random_state=r)
            un = np.where(y == c_out)[0]
            yk = np.array([remap[c] for c in y[idx_k]])
            ymap = dict(zip(idx_k, yk, strict=True))
            g = top_var_idx(X[tr], TOP_GENES)
            m = fit(X[tr], np.array([ymap[i] for i in tr]), g)
            Pc = m.predict_proba(X[ca][:, g])
            yc = np.array([ymap[i] for i in ca])
            P_seen, P_un = m.predict_proba(X[te][:, g]), m.predict_proba(X[un][:, g])
            th = classwise_thresholds(Pc, yc, ALPHA)
            res["unseen_referred_conformal"].append(float(refer(prediction_sets(P_un, th)).mean()))
            res["seen_referred_conformal"].append(float(refer(prediction_sets(P_seen, th)).mean()))
            cut = np.quantile(Pc.max(1), 0.05)
            res["unseen_referred_maxprob"].append(float((P_un.max(1) < cut).mean()))
            res["seen_referred_maxprob"].append(float((P_seen.max(1) < cut).mean()))
            res["auroc_maxprob"].append(auroc(-P_un.max(1), -P_seen.max(1)))
        out["open_set"][cname] = {k: float(np.mean(v)) for k, v in res.items()} | {
            "n_unseen": len(un)
        }
        print(cname, {k: round(v, 3) for k, v in out["open_set"][cname].items()})
    Path(__file__).with_name("results.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
