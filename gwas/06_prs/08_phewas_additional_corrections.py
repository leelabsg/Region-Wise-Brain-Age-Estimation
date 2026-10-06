#!/usr/bin/env python
"""Additional multiple-testing summaries for the PRS PheWAS, computed only from existing raw P values.
No regression is refitted. Each sample is its own family; samples are never pooled.

  1) pair-level global BH:  q_global_pairs_BH = BH over all valid PRS x phecode raw P within a sample
  2) phecode-level omnibus: p_any_PRS_bonferroni10 = min(1, 10 * min raw P over the PRSs tested for that phecode)
                            q_any_PRS_BH = BH over all evaluable phecodes within a sample
     n_PRS_tested records the number of valid tests; a phecode with no valid P gets NA (not 0 or 1).
The per-PRS q (existing 'fdr' column) and beta/SE/OR/P/N are carried over unchanged.

Supplementary Table 18 reports (1) as "Sample-wide BH-adjusted P value" and (2) as "Phecode-level BH-adjusted
P value". These are sensitivity summaries; the primary results remain the per-PRS BH of 05_analyze_ST18_phewas.py.

Usage: python gwas/06_prs/08_phewas_additional_corrections.py OUTDIR name=path [name=path ...]
Each path is a full (all valid tests) PheWAS result table written by 05_analyze_ST18_phewas.py for one sample;
required columns: Region, phecode, p_value, coef, and fdr (per-PRS q).
"""
import sys, os, json
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests

OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
SAMPLES = dict(a.split('=', 1) for a in sys.argv[2:])
ALPHA = 0.05

def bh(p):
    p = np.asarray(p, float); q = np.full(p.shape, np.nan); ok = np.isfinite(p)
    if ok.any(): q[ok] = multipletests(p[ok], method='fdr_bh')[1]
    return q

pairs, phe, summ, checks = [], [], [], []
for name, path in SAMPLES.items():
    d = pd.read_csv(path, sep='\t', dtype={'phecode': str})      # CSV parser: Phenotype may contain newlines
    d['sample'] = name
    d['p_value'] = pd.to_numeric(d.p_value, errors='coerce')
    dup = int(d.duplicated(['Region', 'phecode']).sum())
    if dup: raise SystemExit(f'{name}: {dup} duplicated PRS x phecode rows (validation failure)')
    valid = np.isfinite(d.p_value) & (d.p_value >= 0) & (d.p_value <= 1)
    # sanity: the existing per-PRS q must be reproducible from the same raw P
    reproduced = d.groupby('Region', group_keys=False).apply(lambda g: pd.Series(bh(g.p_value), index=g.index))
    maxdiff = float(np.nanmax(np.abs(reproduced - pd.to_numeric(d.fdr, errors='coerce'))))
    d['q_per_PRS'] = pd.to_numeric(d.fdr, errors='coerce')
    d['q_global_pairs_BH'] = bh(d.p_value.where(valid))
    pairs.append(d)
    # phecode-level omnibus
    g = d[valid].groupby('phecode').agg(n_PRS_tested=('p_value', 'size'), min_raw_p=('p_value', 'min'),
                                        min_q_per_PRS=('q_per_PRS', 'min'))
    allp = pd.DataFrame(index=pd.Index(sorted(d.phecode.unique()), name='phecode')).join(g)
    allp['p_any_PRS_bonferroni10'] = np.minimum(1.0, 10 * allp.min_raw_p)
    allp['q_any_PRS_BH'] = bh(allp.p_any_PRS_bonferroni10)
    allp['sample'] = name
    phe.append(allp.reset_index())
    sig_pp = d.q_per_PRS < ALPHA; sig_g = d.q_global_pairs_BH < ALPHA
    phe_pp = set(d.phecode[sig_pp]); phe_any = set(allp.index[allp.q_any_PRS_BH < ALPHA])
    summ.append(dict(sample=name, source=path, n_tests=len(d), n_valid_p=int(valid.sum()), n_PRS=d.Region.nunique(),
                     n_phecodes=len(allp), n_phecodes_evaluable=int(allp.min_raw_p.notna().sum()),
                     pairs_sig_per_PRS=int(sig_pp.sum()), pairs_sig_global=int(sig_g.sum()),
                     pairs_kept=int((sig_pp & sig_g).sum()), pairs_added=int((~sig_pp & sig_g).sum()),
                     pairs_dropped=int((sig_pp & ~sig_g).sum()),
                     phecodes_any_sig_per_PRS=len(phe_pp), phecodes_sig_any_PRS_omnibus=len(phe_any),
                     phecodes_kept=len(phe_pp & phe_any), phecodes_added=len(phe_any - phe_pp),
                     phecodes_dropped=len(phe_pp - phe_any)))
    checks.append(dict(sample=name, per_PRS_q_reproduced_max_abs_diff=maxdiff, duplicates=dup))

P = pd.concat(pairs, ignore_index=True)
keep = [c for c in ['sample', 'Region', 'phecode', 'Phenotype', 'Category', 'model', 'sex_rule', 'n_case', 'n_total',
                    'coef', 'se', 'odds_ratio', 'ci95_lower', 'ci95_upper', 'p_value', 'q_per_PRS', 'q_global_pairs_BH']
        if c in P.columns]
P[keep].to_csv(f'{OUT}/pairs_with_global_q.tsv', sep='\t', index=False)
H = pd.concat(phe, ignore_index=True)
H[['sample', 'phecode', 'n_PRS_tested', 'min_raw_p', 'min_q_per_PRS', 'p_any_PRS_bonferroni10', 'q_any_PRS_BH']].to_csv(
    f'{OUT}/phecode_any_PRS_summary.tsv', sep='\t', index=False)
S = pd.DataFrame(summ); S.to_csv(f'{OUT}/sample_summary.tsv', sep='\t', index=False)
pd.DataFrame(checks).to_csv(f'{OUT}/validation.tsv', sep='\t', index=False)

# cross-cohort: phecodes significant (omnibus) in more than one sample, and direction of identical PRS-phecode pairs
names = list(SAMPLES)
cross = []
for i, a in enumerate(names):
    for b in names[i + 1:]:
        A = H[(H['sample'] == a) & (H.q_any_PRS_BH < ALPHA)].phecode; B = H[(H['sample'] == b) & (H.q_any_PRS_BH < ALPHA)].phecode
        pa, pb = P[P['sample'] == a], P[P['sample'] == b]
        m = pa.merge(pb, on=['Region', 'phecode'], suffixes=('_a', '_b'))
        for lab, sa, sb in (('per_PRS', 'q_per_PRS_a', 'q_per_PRS_b'), ('global_pairs', 'q_global_pairs_BH_a', 'q_global_pairs_BH_b')):
            both = m[(m[sa] < ALPHA) & (m[sb] < ALPHA)]
            cross.append(dict(sample_a=a, sample_b=b, correction=lab, pairs_sig_in_both=len(both),
                              same_direction=int((np.sign(both.coef_a) == np.sign(both.coef_b)).sum())))
        cross.append(dict(sample_a=a, sample_b=b, correction='any_PRS_omnibus_phecodes',
                          pairs_sig_in_both=len(set(A) & set(B)), same_direction=np.nan))
pd.DataFrame(cross).to_csv(f'{OUT}/cross_sample_overlap.tsv', sep='\t', index=False)
json.dump(dict(samples=SAMPLES, alpha=ALPHA, families='each sample separately; never pooled across samples',
               omnibus='min(1, 10*min raw P) per phecode, BH across evaluable phecodes within sample'),
          open(f'{OUT}/run_config.json', 'w'), indent=1)
print(S.T.to_string()); print(pd.DataFrame(checks).to_string(index=False)); print(pd.DataFrame(cross).to_string(index=False))
