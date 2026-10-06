#!/usr/bin/env python
"""Stage 04 - scaling, clock cohort and splits (UKB-proteomic-age-model-both-sexes.py order):
per-protein MinMaxScaler + median centering fitted on the imputed preprocessing pool -> batch NA/0/7 removed ->
post-imputation >20% missing check on the 204 features -> 70:30 split (seed 3456) and 5-fold OOF / tuning folds.
"""
import os, json, datetime as dt
import numpy as np, pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split, KFold

from _paths import WORK as RUN, PRE, QC, SPL
SEED = 3456
assert os.path.exists(f'{RUN}/logs/s02.done'), 'imputation not complete'

Xi = pd.read_parquet(f'{PRE}/npx_pool_imputed.parquet').set_index('eid').sort_index()
prot = list(Xi.columns)
# ---------------- scaling on the full imputed pool ----------------
rows, S = [], {}
for p in prot:
    v = Xi[p].to_numpy().reshape(-1, 1).copy()
    sc = MinMaxScaler().fit(v); u = sc.transform(v).ravel(); med = np.nanmedian(u)
    S[p] = u - med
    rows.append(dict(protein=p, min=float(sc.data_min_[0]), max=float(sc.data_max_[0]), median_after_minmax=float(med),
                     constant=bool(sc.data_max_[0] == sc.data_min_[0])))
Xs = pd.DataFrame(S, index=Xi.index)
par = pd.DataFrame(rows); par['fit_sample'] = 'preprocessing pool (imputed), n=%d' % len(Xi)
par.to_csv(f'{PRE}/scaler_parameters.tsv', sep='\t', index=False)
assert not par.constant.any(), 'constant protein column'

# ---------------- clock cohort ----------------
flow = pd.read_csv(f'{QC}/sample_flow_stage01.tsv', sep='\t').to_dict('records')
def step(name, n, note=''): flow.append(dict(step=name, n=int(n), note=note)); print(f'  {name}: {n:,} {note}', flush=True)
b = pd.read_parquet(f'{PRE}/batch.parquet').set_index('eid').reindex(Xs.index)
a = pd.read_parquet(f'{PRE}/age_granular.parquet').set_index('eid').reindex(Xs.index)
step('excluded: olink batch missing', b.olink_batch.isna().sum())
step('excluded: batch 0', (b.olink_batch == 0).sum())
step('excluded: batch 7', (b.olink_batch == 7).sum())
keep = b.olink_batch.notna() & ~b.olink_batch.isin([0, 7])
pm = pd.read_csv(f'{QC}/protein_mapping.tsv', sep='\t')
feat = pm.source_column.tolist(); assert len(feat) == 204 and len(set(feat)) == 204
C = Xs.loc[keep, feat]
rm = C.isna().sum(axis=1) > 0.20 * len(feat); pmiss = C.isna().mean() > 0.20
step('excluded: >20% of 204 features missing after imputation', rm.sum()); assert not pmiss.any()
C = C[~rm]
av = a.loc[C.index, 'age_valid']; step('excluded: no valid granular age', (~av).sum()); C = C[av]
C = C.sort_index()
step('clock cohort', len(C))
cohort = pd.concat([a.loc[C.index, ['age_granular_unrounded', 'age_target']], C], axis=1)
cohort.reset_index().to_parquet(f'{PRE}/clock_cohort.parquet', index=False)
json.dump(feat, open(f'{PRE}/feature_order.json', 'w'))
chk = C.describe().T[['min', '50%', 'max']]
print('  scaled features: median of medians %.2e; min %.3f; max %.3f' % (chk['50%'].median(), chk['min'].min(), chk['max'].max()), flush=True)

# ---------------- splits ----------------
ids = C.index.to_numpy()
tr, te = train_test_split(ids, train_size=0.7, random_state=SEED, shuffle=True)
pd.DataFrame(dict(eid=ids, set=np.where(np.isin(ids, tr), 'development_70', 'holdout_30'))).to_csv(f'{SPL}/holdout_split.tsv', sep='\t', index=False)
step('development (70%)', len(tr)); step('holdout (30%)', len(te))
tr_sorted = np.array(tr)   # keep train_test_split order, as the reference passes X_train to KFold in this order
tf = np.empty(len(tr_sorted), int)
for k, (_, v) in enumerate(KFold(5, shuffle=True, random_state=SEED).split(tr_sorted)): tf[v] = k
pd.DataFrame(dict(eid=tr_sorted, tuning_fold=tf)).to_csv(f'{SPL}/tuning_folds.tsv', sep='\t', index=False)
of = np.empty(len(ids), int)
for k, (_, v) in enumerate(KFold(5, shuffle=True, random_state=SEED).split(ids)): of[v] = k
pd.DataFrame(dict(eid=ids, outer_fold=of)).to_csv(f'{SPL}/oof_folds.tsv', sep='\t', index=False)
pd.DataFrame(flow).to_csv(f'{QC}/sample_flow.tsv', sep='\t', index=False)
open(f'{RUN}/logs/s03.done', 'w').write(dt.datetime.utcnow().isoformat() + 'Z\n')
print('stage 03 done', flush=True)
