#!/usr/bin/env python
"""Stage 01 - inputs: schema checks, data-use exclusion, ProtAge-204 panel mapping, Olink batch, granular age.

Inputs (environment):
  PROTAGE_NPX            UK Biobank Olink NPX, instance 0, wide (eid + one column per assay)
  PROTAGE_PLATES         field 30901 (plate id), instance 0: eid, plate
  PROTAGE_BATCHMAP       UK Biobank olink_batch_number.dat (PlateID, Batch)
  PROTAGE_PARTICIPANTS   current participant list with eid and Date_inst_0 (withdrawn participants absent)
  PROTAGE_PHENO          fields 34 (birth year), 52 (birth month), 53 (assessment date) as f.34.0.0, f.52.0.0, f.53.0.0
  PROTAGE_REF_DIR        clone of the published ProtAge model repository, providing
                         files/ProtAge_proteins_2023-12-29.csv and files/olink_names_oct_30_2023.csv
Writes preprocessing/{npx_pool_raw,batch,age_granular}.parquet and qc/{protein_mapping,sample_flow_stage01}.tsv.
Granular age = (assessment date - 1st of birth month) / 365.25; the model target is this age rounded to 2 decimals.
"""
import json, os, datetime as dt
import numpy as np, pandas as pd
from _paths import env, PRE, QC, LOGS, WORK

SRC, PLATES, BATCHMAP = env('PROTAGE_NPX'), env('PROTAGE_PLATES'), env('PROTAGE_BATCHMAP')
PARTICIPANTS, PHENO, REF = env('PROTAGE_PARTICIPANTS'), env('PROTAGE_PHENO'), env('PROTAGE_REF_DIR')
PANEL = f'{REF}/files/ProtAge_proteins_2023-12-29.csv'
OLINK_NAMES = f'{REF}/files/olink_names_oct_30_2023.csv'

flow = []
def step(name, n, note=''):
    flow.append(dict(step=name, n=int(n), note=note)); print(f'  {name}: {n:,} {note}', flush=True)

npx = pd.read_csv(SRC, low_memory=False)
assert npx.columns[0] == 'eid', npx.columns[0]
prot_cols = list(npx.columns[1:])
assert all(pd.api.types.is_numeric_dtype(npx[c]) for c in prot_cols), 'non-numeric protein columns'
assert npx.eid.is_unique, 'duplicate eid in wide matrix'
lower = [c.lower() for c in prot_cols]
assert len(set(lower)) == len(lower), 'case-folding collision among protein columns'
step('source rows (wide, one row per eid)', len(npx))
step('source protein columns', len(prot_cols))

# data-use exclusion before any preprocessing
cur = pd.read_csv(PARTICIPANTS, usecols=['eid', 'Date_inst_0'])
keep = npx.eid.isin(set(cur.eid))
step('excluded: eid absent from current participant list', (~keep).sum(), '(data-use exclusion, applied before imputation)')
npx = npx[keep].reset_index(drop=True)
step('preprocessing pool (imputation + scaling reference)', len(npx))

# ProtAge-204 panel -> source columns (case-insensitive, one-to-one)
official = pd.read_csv(OLINK_NAMES, header=None)[0].tolist()
panel = pd.read_csv(PANEL)
pm = []
for _, r in panel.iterrows():
    g = r['Gene_name']; hits = [c for c in prot_cols if c.lower() == g.lower()]
    pm.append(dict(gene_name=g, uniprot=r['UniProt ID'], source_column=hits[0] if len(hits) == 1 else '',
                   n_matches=len(hits), in_official_list=g in official))
pm = pd.DataFrame(pm)
assert pm.gene_name.is_unique and len(pm) == 204
pm['observed_n'] = [int(npx[c].notna().sum()) if c else 0 for c in pm.source_column]
assert (pm.n_matches == 1).all() and (pm.observed_n > 0).all(), 'panel mapping incomplete'
pm.to_csv(f'{QC}/protein_mapping.tsv', sep='\t', index=False)

# Olink batch: field 30901 plate -> PlateID/Batch map
pl = pd.read_csv(PLATES); pl.columns = ['eid', 'plate']
assert pl.eid.is_unique
bm = pd.read_csv(BATCHMAP, sep=r'\s+')
assert list(bm.columns) == ['PlateID', 'Batch'] and bm.PlateID.is_unique
b = npx[['eid']].merge(pl, on='eid', how='left').merge(bm.rename(columns={'PlateID': 'plate', 'Batch': 'olink_batch'}),
                                                        on='plate', how='left')
b.to_parquet(f'{PRE}/batch.parquet', index=False)

# granular age
ph = pd.read_csv(PHENO, usecols=['f.eid', 'f.34.0.0', 'f.52.0.0', 'f.53.0.0'])
ph.columns = ['eid', 'birth_year', 'birth_month', 'visit_date_i0']
ph = ph[ph.eid.isin(set(npx.eid))]
a = npx[['eid']].merge(ph, on='eid', how='left')
a['visit_date_i0'] = pd.to_datetime(a.visit_date_i0, errors='coerce')
ok = a.birth_year.notna() & a.birth_month.between(1, 12) & a.visit_date_i0.notna()
est_birth = pd.to_datetime(dict(year=a.birth_year.where(ok, 1900).astype(int), month=a.birth_month.where(ok, 1).astype(int), day=1))
a['age_granular_unrounded'] = np.where(ok, (a.visit_date_i0 - est_birth).dt.days / 365.25, np.nan)
a['age_target'] = np.round(a.age_granular_unrounded, 2)
a['age_valid'] = ok
a[['eid', 'age_granular_unrounded', 'age_target', 'age_valid']].to_parquet(f'{PRE}/age_granular.parquet', index=False)

npx.to_parquet(f'{PRE}/npx_pool_raw.parquet', index=False)
pd.DataFrame(flow).to_csv(f'{QC}/sample_flow_stage01.tsv', sep='\t', index=False)
open(f'{LOGS}/s01.done', 'w').write(dt.datetime.utcnow().isoformat() + 'Z\n')
print('stage 01 done', flush=True)
