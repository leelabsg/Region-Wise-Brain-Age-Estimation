#!/usr/bin/env python
"""Stage 02 - protein-only missingness QC, applied in order:
 1) drop protein assays with missingness >= 50% in the preprocessing pool
 2) participant missingness over the QC-passed proteins only; drop participants with >= 54%
 3) re-check assay missingness after participant removal; drop assays with >= 50%
Only protein columns enter the missingness denominators. All 204 panel proteins must survive.
Writes preprocessing/knn_input.parquet and qc/{protein_qc,sample_flow_stage01}.tsv, qc/qc_summary.json.
"""
import json, datetime as dt
import pandas as pd
from _paths import PRE, QC, LOGS

X = pd.read_parquet(f'{PRE}/npx_pool_raw.parquet').set_index('eid').sort_index()
flow = pd.read_csv(f'{QC}/sample_flow_stage01.tsv', sep='\t').to_dict('records')
def step(name, n, note=''): flow.append(dict(step=name, n=int(n), note=note)); print(f'  {name}: {n:,} {note}', flush=True)

m1 = X.isna().mean()
drop1 = sorted(m1.index[m1 >= 0.50]); X1 = X.drop(columns=drop1)
step('protein assays dropped (missing >= 50%, pool)', len(drop1), str(drop1))
rm = X1.isna().mean(axis=1)
drop_p = rm.index[rm >= 0.54]; X2 = X1.drop(index=drop_p)
step('participants dropped (missing >= 54% of QC-passed proteins)', len(drop_p))
m3 = X2.isna().mean()
drop3 = sorted(m3.index[m3 >= 0.50]); X3 = X2.drop(columns=drop3)
step('protein assays dropped at re-check (missing >= 50%)', len(drop3), str(drop3))
step('KNN input participants', X3.shape[0]); step('KNN input proteins', X3.shape[1])

pm = pd.read_csv(f'{QC}/protein_mapping.tsv', sep='\t')
missing_panel = [c for c in pm.source_column if c not in X3.columns]
assert not missing_panel, f'panel proteins removed by QC: {missing_panel}'

X3.reset_index().to_parquet(f'{PRE}/knn_input.parquet', index=False)
pd.DataFrame(dict(protein=X.columns, missing_pool=m1.values, dropped_step1=X.columns.isin(drop1),
                  dropped_step3=X.columns.isin(drop3))).to_csv(f'{QC}/protein_qc.tsv', sep='\t', index=False)
pd.DataFrame(flow).to_csv(f'{QC}/sample_flow_stage01.tsv', sep='\t', index=False)
json.dump(dict(thresholds=dict(assay_ge=0.50, participant_ge=0.54, assay_recheck_ge=0.50),
               n_participants=X3.shape[0], n_proteins=X3.shape[1], proteins_dropped=drop1 + drop3,
               n_participants_dropped=len(drop_p)), open(f'{QC}/qc_summary.json', 'w'), indent=1)
open(f'{LOGS}/k01.done', 'w').write(dt.datetime.utcnow().isoformat() + 'Z\n')
