# Software Requirements Specification – OmicsGuard

- SR-1 Features: the 2000 highest-variance genes on the training part; standardisation; multinomial logistic regression (C = 0.05).
- SR-2 `conformal.classwise_thresholds`, `prediction_sets`: class-wise split-conformal sets at alpha 0.05; `refer`: set empty or several labels; `auroc`: Mann–Whitney AUROC with tie handling.
- SR-3 `research/run_study.py` (1): 20 random splits 50/25/25, coverage per class and referral rate; (2) leave-one-tumour-type-out: the type is absent from training and calibration and appears only in the test set; 10 repeats per type; comparison with a max-probability rule (threshold = 5 % quantile of calibration max-probabilities).
- SR-4 Two unit tests.
