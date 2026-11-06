# Business Requirements Specification – OmicsGuard

| ID | Requirement |
|----|-------------|
| BR-1 | Classify tumour type from RNA-Seq with a per-class coverage guarantee, and report referral to an expert. |
| BR-2 | Detect samples from a tumour type that the model never saw, instead of forcing a label. |
| BR-3 | Compare conformal referral with a simple maximum-probability rejection rule. |
| BR-4 | Report real-data results only; no diagnostic claim. |

Out of scope: clinical deployment, survival or treatment prediction, other cancer types, batch correction across platforms.
