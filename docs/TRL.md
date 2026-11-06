# TRL assessment – OmicsGuard

**Claimed level: TRL 4** (laboratory validation on public real data). Not TRL 5: one cohort-level data release, no external platform, no clinical setting.

Data: TCGA PanCancer RNA-Seq via UCI 401 (BRCA 300, KIRC 146, LUAD 141, PRAD 136, COAD 78).

## Closed-set (20 random splits, alpha = 0.05)
Accuracy 0.998; coverage 0.964; referral rate 0.036; class-wise coverage: BRCA 0.957, COAD 0.945, KIRC 0.972, LUAD 0.975, PRAD 0.969.

## Open-set: one tumour type removed from training and calibration
Referral = conformal set empty or several labels. Max-prob rule = max probability below the 5 % calibration quantile.

| Unseen type | n unseen | AUROC (max-prob) | Unseen referred (conformal) | Seen referred (conformal) | Unseen referred (max-prob) | Seen referred (max-prob) |
|---|---|---|---|---|---|---|
| BRCA | 300 | 0.998 | 0.880 | 0.029 | 1.000 | 0.051 |
| COAD | 78 | 0.991 | 0.500 | 0.029 | 1.000 | 0.049 |
| KIRC | 146 | 0.999 | 0.569 | 0.031 | 1.000 | 0.045 |
| LUAD | 141 | 0.998 | 0.999 | 0.034 | 1.000 | 0.048 |
| PRAD | 136 | 0.993 | 0.987 | 0.047 | 1.000 | 0.067 |

## Findings and caveats
- The closed-set task is nearly solved (accuracy 0.998); class-wise coverage is 0.94–0.98 for a 0.95 target (small calibration sets, 10–40 samples per class, so ±0.03 is expected).
- Conformal referral of an unseen tumour type is **uneven**: LUAD and PRAD are referred almost always (0.999 and 0.987) but COAD only 0.50 and KIRC 0.57. Half of the unseen COAD samples and 43 % of unseen KIRC samples receive a confident single wrong label.
- The max-probability rule with a 5 % false-referral budget refers every unseen sample in all five experiments (1.000) at 4.5–6.7 % false referrals and AUROC 0.99–1.00: on this data it is both simpler and better than conformal referral for novelty detection. Conformal sets guarantee coverage of known classes; they do not guarantee detection of unknown ones.
- Cancer types here are very distinct in expression; the result may not hold for closely related subtypes. The 2000-gene filter and C are fixed, not tuned.
- Samples are TCGA tumours without information on batch or centre; test and calibration may share batches with training data, so closed-set accuracy may be optimistic.
