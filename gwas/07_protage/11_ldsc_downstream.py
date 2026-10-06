#!/usr/bin/env python
"""Stage 11 - LDSC for the ProtAge-204 GWAS and the external-ageing correlation family.

  preprocess (04_ldsc/01_preprocess.py: 0.01 < AF_Allele2 < 0.99, per chromosome)
  -> munge_sumstats.py per chromosome (--a1 Allele2 --a2 Allele1, --merge-alleles own SNP list, --N = the N column)
  -> 04_ldsc/03_merge_munged.py -> ldsc --h2
  -> for each BAG: 04_ldsc/align_flip.py (BAG -> ProtAge orientation), ldsc --rg ProtAge,BAG
  -> telomere length PC1: align_flip.py, ldsc --rg ProtAge,telomere
  -> BH over the 21 external-ageing tests: telomere x 10 BAG (from the existing telomere x BAG --rg logs),
     ProtAge x 10 BAG, ProtAge x telomere; P is the two-sided normal P of the printed Z
  -> phenotypic partial correlations ProtAge x BAG on the raw (pre-INT) gaps, adjusted for sex, centred age at
     protein sampling and the centred protein-MRI interval; BH across the ten BAGs
Cross-trait intercepts are unconstrained.

Environment: LDSC_DIR, LDSC_PYTHON, LDSC_REF_LD_CHR, LDSC_W_LD_CHR, LDSC_TELOMERE_SUMSTATS,
  PROTAGE_SUMSTATS            merged SAIGE step-2 output (14 columns, Allele2 = effect allele)
  PROTAGE_LDSC_DIR            working/output directory
  LDSC_BAG_SUMSTATS_TMPL      munged BAG sumstats, {region} substituted
  LDSC_TELOMERE_BAG_LOG_TMPL  existing ldsc --rg log for telomere x BAG, {region} substituted
  PROTAGE_GWAS_PHENO          output of 08_build_gwas_pheno.py
  BAG_PHENO_FILE              BAG table with IID, Age (at MRI), Sex, {region}_corrected_delta_age
"""
import os, re, sys, json, math, subprocess
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import env

HERE = env('PROTAGE_LDSC_DIR'); MERGED = env('PROTAGE_SUMSTATS'); NAME = 'protage204'
LDSC_DIR, LDSC_PY = env('LDSC_DIR'), env('LDSC_PYTHON', 'python')
REF_LD, W_LD = env('LDSC_REF_LD_CHR'), env('LDSC_W_LD_CHR')
TEL, PADS, TEL_LOGS = env('LDSC_TELOMERE_SUMSTATS'), env('LDSC_BAG_SUMSTATS_TMPL'), env('LDSC_TELOMERE_BAG_LOG_TMPL')
L4 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '04_ldsc')
REGIONS = ['global', 'caudate', 'cerebellum', 'frontal_lobe', 'insula', 'occipital_lobe', 'parietal_lobe', 'putamen', 'temporal_lobe', 'thalamus']
DISPLAY = {'global': 'Whole Brain', 'caudate': 'Caudate', 'cerebellum': 'Cerebellum', 'frontal_lobe': 'Frontal Lobe', 'insula': 'Insula',
           'occipital_lobe': 'Occipital Lobe', 'parietal_lobe': 'Parietal Lobe', 'putamen': 'Putamen', 'temporal_lobe': 'Temporal Lobe', 'thalamus': 'Thalamus'}
TELDISP, PROTDISP = 'Telomere Length PC1', 'ProtAge-204'

def run(cmd, out):
    if subprocess.run(cmd, stdout=open(out, 'w'), stderr=subprocess.STDOUT).returncode != 0:
        sys.exit(f'FAIL: {" ".join(cmd[:3])} ... (see {out})')
for d in ('files', 'snplist_A1effect', 'munge_output', 'sumstats/aligned', 'h2', 'rg', 'tables', 'phenotypic', 'logs'):
    os.makedirs(f'{HERE}/{d}', exist_ok=True)
SUM = f'{HERE}/sumstats/{NAME}.sumstats.gz'

# 1. preprocess
n = pd.read_csv(MERGED, sep=r'\s+', usecols=['N']).N
if n.min() != n.max(): sys.exit('FAIL: N column not constant')
N = int(n.iloc[0]); json.dump(dict(N=N, rows=int(len(n))), open(f'{HERE}/logs/merged_N.json', 'w'))
run([sys.executable, f'{L4}/01_preprocess.py', '--save_path', f'{HERE}/files', '--save_file_name', NAME, '--gwas_results_path', MERGED], f'{HERE}/logs/1_preprocess.out')

# 2. munge (A1 = Allele2 = effect allele) + merge
procs = []
for c in range(1, 23):
    inp = f'{HERE}/files/chr{c}_{NAME}_GWAS_SummaryStatistics.txt'
    with open(inp) as fi, open(f'{HERE}/snplist_A1effect/chr{c}_snplist.txt', 'w') as fo:
        h = fi.readline().rstrip('\n').split('\t'); iS, i1, i2 = h.index('MarkerID'), h.index('Allele1'), h.index('Allele2')
        fo.write('SNP\tA1\tA2\n')
        for line in fi:
            v = line.rstrip('\n').split('\t'); fo.write(f'{v[iS]}\t{v[i2]}\t{v[i1]}\n')
    cmd = [LDSC_PY, f'{LDSC_DIR}/munge_sumstats.py', '--sumstats', inp, '--N', str(N), '--snp', 'MarkerID', '--a1', 'Allele2', '--a2', 'Allele1',
           '--out', f'{HERE}/munge_output/{NAME}_chr{c}_munge', '--merge-alleles', f'{HERE}/snplist_A1effect/chr{c}_snplist.txt']
    procs.append((c, subprocess.Popen(cmd, stdout=open(f'{HERE}/munge_output/{NAME}_chr{c}.out', 'w'), stderr=subprocess.STDOUT)))
for c, p in procs:
    if p.wait() != 0: sys.exit(f'FAIL: munge chr{c}')
run([sys.executable, f'{L4}/03_merge_munged.py', '--input_dir', f'{HERE}/munge_output', '--output_file', SUM, '--munge_file_name', NAME], f'{HERE}/logs/2_merge.out')

def ldsc(args, prefix):
    run([LDSC_PY, f'{LDSC_DIR}/ldsc.py'] + args + ['--ref-ld-chr', REF_LD, '--w-ld-chr', W_LD, '--out', prefix], prefix + '.stdout')

# 3. h2
ldsc(['--h2', SUM], f'{HERE}/h2/heritability_{NAME}')

# 4. rg
def align(src, out):
    run([sys.executable, f'{L4}/align_flip.py', '--source_file', src, '--reference_file', SUM, '--output_file', out], out + '.align.out')
for r in REGIONS:
    al = f'{HERE}/sumstats/aligned/{NAME}_{r}_aligned.sumstats.gz'; align(PADS.replace('{region}', r), al)
    ldsc(['--rg', f'{SUM},{al}'], f'{HERE}/rg/gene_corr_protein_{r}')
al = f'{HERE}/sumstats/aligned/{NAME}_telomere_aligned.sumstats.gz'; align(TEL, al)
ldsc(['--rg', f'{SUM},{al}'], f'{HERE}/rg/gene_corr_{NAME}_telomere')

# 5. 21-test BH
flt = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
RE = {'rg': re.compile(r"Genetic Correlation.*?Genetic Correlation:\s*(" + flt + r")\s*\((" + flt + r")\)", re.S),
      'z': re.compile(r"Genetic Correlation.*?Z-score:\s*(" + flt + r")", re.S),
      'p': re.compile(r"Genetic Correlation.*?P:\s*(" + flt + r")", re.S),
      'gencov_int': re.compile(r"Genetic Covariance.*?Intercept:\s*(" + flt + r")\s*\((" + flt + r")\)", re.S)}
def parse_rg(path):
    t = open(path).read()
    if 'Analysis finished' not in t: sys.exit(f'FAIL: {path} incomplete')
    m = RE['rg'].search(t); gi = RE['gencov_int'].search(t)
    return dict(rg=float(m.group(1)), rg_SE=float(m.group(2)), Z=float(RE['z'].search(t).group(1)), P_ldsc=float(RE['p'].search(t).group(1)),
                gencov_intercept=float(gi.group(1)), gencov_intercept_SE=float(gi.group(2)))
def bh(pvals):
    n = len(pvals); order = sorted(range(n), key=lambda i: pvals[i]); adj = [0.0] * n; prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]; prev = min(prev, pvals[i] * n / (rank + 1)); adj[i] = prev
    return adj
rows = [dict(parse_rg(TEL_LOGS.replace('{region}', r)), Trait1=TELDISP, Trait2=DISPLAY[r]) for r in REGIONS]
rows += [dict(parse_rg(f'{HERE}/rg/gene_corr_protein_{r}.log'), Trait1=PROTDISP, Trait2=DISPLAY[r]) for r in REGIONS]
rows.append(dict(parse_rg(f'{HERE}/rg/gene_corr_{NAME}_telomere.log'), Trait1=PROTDISP, Trait2=TELDISP))
for o in rows: o['P'] = math.erfc(abs(o['Z']) / math.sqrt(2.0))
for o, q in zip(rows, bh([o['P'] for o in rows])): o['FDR_BH_21'] = q
RG = pd.DataFrame(rows); RG.to_csv(f'{HERE}/tables/external_ageing_rg_21_BH.tsv', sep='\t', index=False)

# 6. phenotypic partial correlations (raw gaps; sex, centred age at protein sampling, centred protein-MRI interval)
def residualize(C, v):
    beta, *_ = np.linalg.lstsq(C, v, rcond=None); return v - C @ beta
def pearson(x, y):
    xc, yc = x - x.mean(), y - y.mean(); return float(np.sum(xc * yc) / np.sqrt(np.sum(xc * xc) * np.sum(yc * yc)))
def stats_from_r(r, n, k):
    df = n - k - 2; rc = min(max(r, -1.0 + 1e-12), 1.0 - 1e-12); t = rc * np.sqrt(df / (1.0 - rc * rc))
    p = float(2.0 * stats.t.sf(abs(t), df)); z = np.arctanh(rc); se = 1.0 / np.sqrt(n - k - 3)
    return p, float(np.tanh(z - 1.959963984540054 * se)), float(np.tanh(z + 1.959963984540054 * se))
prot = pd.read_csv(env('PROTAGE_GWAS_PHENO'), dtype={'IID': str})
reg = pd.read_csv(env('BAG_PHENO_FILE'), dtype={'IID': str})
m = prot[['IID', 'Age', 'Sex', 'protage204_gap_residual']].merge(reg[['IID', 'Age', 'Sex'] + [f'{t}_corrected_delta_age' for t in REGIONS]],
                                                              on='IID', suffixes=('_protein', '_brain'))
assert m.IID.is_unique and (m.Sex_protein.astype(int).to_numpy() == m.Sex_brain.astype(int).to_numpy()).all()
ageP = m.Age_protein.to_numpy(float); interval = m.Age_brain.to_numpy(float) - ageP
C = np.column_stack([np.ones(len(m)), m.Sex_protein.to_numpy(float), ageP - ageP.mean(), interval - interval.mean()])
X = m.protage204_gap_residual.to_numpy(float); ph = []
for t in REGIONS:
    r = pearson(residualize(C, X), residualize(C, m[f'{t}_corrected_delta_age'].to_numpy(float)))
    p, lo, hi = stats_from_r(r, len(m), 3)
    ph.append(dict(Trait1=PROTDISP, Trait2=DISPLAY[t], phenotypic_r=r, p=p, n=len(m), ci_low=lo, ci_high=hi))
PH = pd.DataFrame(ph); PH['q_bh_10'] = bh(PH.p.tolist())
PH.to_csv(f'{HERE}/phenotypic/protage204_bag_phenotypic_corr.tsv', sep='\t', index=False)
json.dump(dict(n_overlap=len(m), interval_years=dict(mean=float(interval.mean()), sd=float(interval.std(ddof=1)))),
          open(f'{HERE}/phenotypic/overlap_summary.json', 'w'), indent=1)
print(RG[['Trait1', 'Trait2', 'rg', 'rg_SE', 'P', 'FDR_BH_21']].to_string(index=False)); print(PH.to_string(index=False))
