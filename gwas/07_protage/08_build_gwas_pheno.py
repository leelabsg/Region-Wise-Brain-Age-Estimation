#!/usr/bin/env python
"""Stage 08 - GWAS sample and phenotype/covariate file for SAIGE.

Sample = out-of-fold clock cohort  n  self-reported British (field 21000, instance 0)
         n  the SAIGE step-1 PLINK set (UKB_GRM_PLINK.fam)  n  the imputed BGEN sample (UKB_SAMPLE_FILE).
The same sample is used for step 1 and step 2.
Columns: IID, Age (age_target: granular age at protein sampling, 2 dp - not the predicted age), Sex (0 female,
1 male), predicted_protein_age, raw_protage204_gap, protage204_gap_residual, protage204_gap_int (the GWAS
phenotype), PC1-PC10.

Environment: PROTAGE_WORK_DIR, UKB_GRM_PLINK, UKB_SAMPLE_FILE,
  PROTAGE_ETHNICITY   CSV with 'Participant ID' and 'Ethnic background | Instance 0'
  PROTAGE_COVARIATES  CSV with eid, Sex (Female/Male) and PC1-PC10 (UK Biobank genetic principal components)
  PROTAGE_GWAS_PHENO  output CSV
  PROTAGE_EXPECTED_N  optional assertion on the final N (38806 in the published run)
"""
import os
import numpy as np, pandas as pd
from _paths import env, WORK

PRED = f'{WORK}/results/protage204_predictions.parquet'
FAM, SAMPLE = env('UKB_GRM_PLINK') + '.fam', env('UKB_SAMPLE_FILE')
ETH, MAT, OUT = env('PROTAGE_ETHNICITY'), env('PROTAGE_COVARIATES'), env('PROTAGE_GWAS_PHENO')
pcs = [f'PC{i}' for i in range(1, 11)]
p = pd.read_parquet(PRED).rename(columns={'eid': 'IID'})
fam = set(pd.read_csv(FAM, sep=r'\s+', header=None, usecols=[1])[1])
bg = set(pd.read_csv(SAMPLE, sep=r'\s+', skiprows=2, header=None, usecols=[0])[0])
eth = pd.read_csv(ETH, usecols=['Participant ID', 'Ethnic background | Instance 0'])
brit = set(eth.loc[eth.iloc[:, 1].eq('British'), 'Participant ID'])
flow = [('OOF clock cohort', len(p))]
p = p[p.IID.isin(brit)]; flow.append(('self-reported British (21000 i0)', len(p)))
p = p[p.IID.isin(fam)]; flow.append(('in step-1 PLINK .fam', len(p)))
p = p[p.IID.isin(bg)]; flow.append(('in step-2 imputed BGEN sample', len(p)))
m = pd.read_csv(MAT, usecols=['eid', 'Sex'] + pcs).rename(columns={'eid': 'IID'}); m['Sex'] = m.Sex.map({'Female': 0, 'Male': 1})
d = p.merge(m, on='IID', how='left', validate='1:1')
out = pd.DataFrame({'IID': d.IID, 'Age': d.age_target, 'Sex': d.Sex, 'predicted_protein_age': d.protage_oof,
                    'raw_protage204_gap': d.protage_gap_raw, 'protage204_gap_residual': d.protage_gap_residual,
                    'protage204_gap_int': d.protage_gap_int})
for c in pcs: out[c] = d[c]
need = ['Age', 'Sex', 'protage204_gap_int'] + pcs
assert out.IID.is_unique
assert out.IID.isin(fam).all() and out.IID.isin(bg).all() and out.IID.isin(brit).all()
assert np.isfinite(out[need].to_numpy(float)).all(), 'missing phenotype/covariate'
assert set(out.Sex.unique()) <= {0, 1}
if os.environ.get('PROTAGE_EXPECTED_N'):
    assert len(out) == int(os.environ['PROTAGE_EXPECTED_N']), len(out)
out.to_csv(OUT, index=False)
pd.DataFrame(flow, columns=['step', 'n']).to_csv(f'{WORK}/qc/gwas_sample_flow.tsv', sep='\t', index=False)
print(flow)
