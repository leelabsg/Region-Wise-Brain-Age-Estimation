#!/usr/bin/env python
"""Stage 07 - proteomic-age gap definitions, computed once on the full out-of-fold clock cohort.
  protage_gap_raw      = protage_oof - age_target
  protage_gap_residual = residual of OLS(protage_gap_raw ~ 1 + age_target)
  protage_gap_int      = rank-based inverse-normal transform of protage_gap_residual
                         (rankdata 'average', norm.ppf((rank - 0.5) / n); regionbae/postprocess.py)
Genetic eligibility is applied afterwards, in 08_build_gwas_pheno.py.
Writes results/protage204_predictions.parquet and results/calibration.json.
"""
import os, json, datetime as dt
import numpy as np, pandas as pd
from scipy.stats import rankdata, norm
from _paths import WORK as RUN

assert os.path.exists(f'{RUN}/logs/s05.done')

def inverse_normal_transformation(values):
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    valid_mask = np.isfinite(values)
    valid_values = values[valid_mask]
    if len(valid_values) == 0:
        return out
    ranks = rankdata(valid_values, method='average')
    transformed = norm.ppf((ranks - 0.5) / len(valid_values))
    out[valid_mask] = transformed
    return out

O = pd.read_parquet(f'{RUN}/results/oof_predictions_stage05.parquet')
O['protage_gap_raw'] = O.protage_oof - O.age_target
A = np.column_stack([np.ones(len(O)), O.age_target.to_numpy()])
coef, *_ = np.linalg.lstsq(A, O.protage_gap_raw.to_numpy(), rcond=None)
O['protage_gap_residual'] = O.protage_gap_raw - A @ coef
O['protage_gap_int'] = inverse_normal_transformation(O.protage_gap_residual)
assert np.allclose(O.protage_gap_raw, O.protage_oof - O.age_target)
assert abs(np.corrcoef(O.protage_gap_residual, O.age_target)[0, 1]) < 1e-8 and abs(O.protage_gap_residual.mean()) < 1e-8
cols = ['eid', 'age_granular_unrounded', 'age_target', 'outer_fold', 'protage_oof', 'protage_gap_raw',
        'protage_gap_residual', 'protage_gap_int']
O[cols].to_parquet(f'{RUN}/results/protage204_predictions.parquet', index=False)
json.dump(dict(model='protage_gap_raw ~ 1 + age_target', intercept=float(coef[0]), slope=float(coef[1]),
               reference_n=len(O), int_convention="rankdata(method='average'); norm.ppf((rank-0.5)/n)"),
          open(f'{RUN}/results/calibration.json', 'w'), indent=1)
open(f'{RUN}/logs/s06.done', 'w').write(dt.datetime.utcnow().isoformat() + 'Z\n')
print('clock cohort N:', len(O), flush=True)
