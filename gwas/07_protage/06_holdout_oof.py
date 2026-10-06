#!/usr/bin/env python
"""Stage 06 - (1) fit on the 70% development set with the selected hyperparameters, predict the 30% holdout once;
(2) full clock cohort 5-fold OOF with the same fixed hyperparameters (lgbm_cv_regression): each participant gets one
prediction from the fold model that excluded them. Models saved as LightGBM text; reload predictions compared.
Metrics use the unadjusted predicted age.
"""
import os, json, datetime as dt
import numpy as np, pandas as pd, lightgbm as lgb
from scipy.stats import pearsonr
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from concurrent.futures import ProcessPoolExecutor

from _paths import WORK as RUN
SEED = 3456
THREADS_PER_FIT = int(os.environ.get('THREADS_PER_FIT', '16'))
assert os.path.exists(f'{RUN}/logs/s04.done')
D = pd.read_parquet(f'{RUN}/preprocessing/clock_cohort.parquet').set_index('eid')
FEAT = json.load(open(f'{RUN}/preprocessing/feature_order.json'))
sel = json.load(open(f'{RUN}/optuna/selected_hyperparameters.json'))
P = dict(sel['best_params']); P.update(sel['base_params']); P.update(n_jobs=THREADS_PER_FIT)
assert P['boosting_type'] == 'dart' and P['n_estimators'] == 5000

def metrics(y, p):
    return dict(N=len(y), pearson_r=float(pearsonr(y, p)[0]), R2=float(r2_score(y, p)),
                MAE=float(mean_absolute_error(y, p)), RMSE=float(np.sqrt(mean_squared_error(y, p))))

def fit_and_save(args):
    name, tr_ids, te_ids = args
    m = lgb.LGBMRegressor(**P).fit(D.loc[tr_ids, FEAT], D.loc[tr_ids, 'age_target'])
    path = f'{RUN}/models/{name}.txt'; m.booster_.save_model(path)
    pred = m.predict(D.loc[te_ids, FEAT])
    re = lgb.Booster(model_file=path).predict(D.loc[te_ids, FEAT].to_numpy())
    assert np.allclose(pred, re, rtol=0, atol=1e-9), f'reload mismatch {name}'
    assert list(m.booster_.feature_name()) == FEAT, 'feature order'
    return name, te_ids, pred, m.booster_.num_trees()

hs = pd.read_csv(f'{RUN}/splits/holdout_split.tsv', sep='\t')
of = pd.read_csv(f'{RUN}/splits/oof_folds.tsv', sep='\t')
jobs = [('holdout_model_dev70', hs.eid[hs.set == 'development_70'].to_numpy(), hs.eid[hs.set == 'holdout_30'].to_numpy())]
for k in range(5):
    jobs.append((f'oof_fold{k}', of.eid[of.outer_fold != k].to_numpy(), of.eid[of.outer_fold == k].to_numpy()))
with ProcessPoolExecutor(len(jobs)) as ex:
    out = list(ex.map(fit_and_save, jobs))

rows, oof_parts = [], []
for name, te, pred, ntree in out:
    y = D.loc[te, 'age_target'].to_numpy()
    rows.append(dict(set=name, n_trees=ntree, **metrics(y, pred)))
    if name.startswith('oof_fold'):
        oof_parts.append(pd.DataFrame(dict(eid=te, outer_fold=int(name[-1]), protage_oof=pred)))
O = pd.concat(oof_parts).sort_values('eid')
# coverage / uniqueness / train-test disjointness
assert O.eid.is_unique and len(O) == len(D) and set(O.eid) == set(D.index), 'OOF coverage'
assert np.isfinite(O.protage_oof).all()
for k in range(5):
    assert not set(of.eid[of.outer_fold != k]) & set(O.eid[O.outer_fold == k]), 'train/test overlap'
O = O.merge(D[['age_granular_unrounded', 'age_target']], left_on='eid', right_index=True)
rows.append(dict(set='pooled_oof', n_trees=np.nan, **metrics(O.age_target.to_numpy(), O.protage_oof.to_numpy())))
pd.DataFrame(rows).to_csv(f'{RUN}/results/prediction_metrics.tsv', sep='\t', index=False)
O.to_parquet(f'{RUN}/results/oof_predictions_stage05.parquet', index=False)
json.dump(dict(params=P, feature_order=FEAT), open(f'{RUN}/models/model_config.json', 'w'), indent=1)
open(f'{RUN}/logs/s05.done', 'w').write(dt.datetime.utcnow().isoformat() + 'Z\n')
print(pd.DataFrame(rows).round(4).to_string(index=False), flush=True)
