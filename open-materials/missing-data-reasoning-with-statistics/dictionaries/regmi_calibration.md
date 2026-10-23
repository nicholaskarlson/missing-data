# Regression-MI calibration case

Fictional numerical calibration records, not Tonmick participants or measured
Likert scores. row_id: P001..P032; x: completely observed fixed continuous design
variable; y: unbounded continuous response, blank for eight unobserved values.
The 24 observed response numbers are hand-authored teaching inputs. There is no
complete-y file or known generating coefficient. The fixed missingness pattern
does not establish MAR. Model-based inference assumes ignorability given x.
The target is the conditional model slope, not the finite mean of these 32 rows.
See model-specs/regmi.json for all model, variance and failure conventions.
