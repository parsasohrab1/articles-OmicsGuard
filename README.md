# articles-OmicsGuard

Class-wise conformal tumour-type classification from RNA-Seq with referral of samples from an **unseen** tumour type. Real data:
TCGA PanCancer RNA-Seq as distributed by UCI (801 tumours, 5 types, 20 531 genes). TRL 4 (public real data, laboratory validation only).

```
pip install -e ".[dev]" && python -m pytest
# data: https://archive.ics.uci.edu/dataset/401 (TCGA-PANCAN-HiSeq-801x20531; not committed)
python research/run_study.py <dir with data.csv and labels.csv>
```
