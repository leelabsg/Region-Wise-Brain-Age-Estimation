#!/usr/bin/env python
"""Relatedness pruning for the PRS-PheWAS sensitivity analysis (ST18 continued). Weights are not re-estimated.

Related pair: KING kinship > 0.0884 in the UK Biobank-provided relatedness file (no recomputation).
 1) remove target participants related to anyone in the PRS-derivation sample (n = 30,800)
 2) within the remaining target, greedily remove the highest-degree participant of the relatedness graph until
    no related pair remains; degree ties are broken by a fixed random priority (seed 20260930). No phenotype
    information is used.
The keep list is then passed to 05_analyze_ST18_phewas.py as PRS_PHEWAS_KEEP, with the same models and BH
correction as the primary run. Published: White British non-imaging 389,045 -> 359,334; held-out imaging
10,267 -> 9,997.

    python gwas/06_prs/09_relatedness_prune.py <ukb_rel.dat> <derivation_ids.txt> <target_ids.txt> <out_prefix>

ID files: one participant id per line. Writes <out_prefix>_keep_ids.txt and <out_prefix>_summary.json.
"""
import sys, json, heapq
import numpy as np, pandas as pd

THR, SEED = 0.0884, 20260930
REL, DER, TGT, OUT = sys.argv[1:5]
rel = pd.read_csv(REL, sep=r'\s+'); r = rel[rel.Kinship > THR][['ID1', 'ID2', 'Kinship']]
der = set(pd.read_csv(DER, header=None)[0].astype(int))
T = set(pd.read_csv(TGT, header=None)[0].astype(int)); n0 = len(T)
assert not (T & der), 'participant overlap with derivation'
ext = r[(r.ID1.isin(T) & r.ID2.isin(der)) | (r.ID2.isin(T) & r.ID1.isin(der))]
ext_rm = (set(ext.ID1) | set(ext.ID2)) & T
T1 = T - ext_rm
e = r[r.ID1.isin(T1) & r.ID2.isin(T1)]
adj = {}
for a, b in zip(e.ID1, e.ID2):
    if a == b: continue
    adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
nodes = sorted(adj)                                              # deterministic node order
prio = dict(zip(nodes, np.random.default_rng(SEED).permutation(len(nodes))))
heap = [(-len(adj[v]), prio[v], v) for v in nodes]; heapq.heapify(heap)
removed = set()
while heap:
    d, p, v = heapq.heappop(heap)
    if v in removed or -d != len(adj[v]): continue               # stale entry
    if len(adj[v]) == 0: continue
    removed.add(v)
    for u in adj[v]:
        adj[u].discard(v); heapq.heappush(heap, (-len(adj[u]), prio[u], u))
    adj[v] = set()
T2 = T1 - removed
left = r[r.ID1.isin(T2) & r.ID2.isin(T2) & (r.ID1 != r.ID2)]
assert len(left) == 0, f'{len(left)} related pairs remain'
pd.Series(sorted(T2)).to_csv(f'{OUT}_keep_ids.txt', index=False, header=False)
s = dict(original_N=n0, removed_related_to_derivation=len(ext_rm), removed_internal_pruning=len(removed),
         final_N=len(T2), final_internal_related_pairs=0, threshold=f'Kinship > {THR}', tie_break_seed=SEED)
json.dump(s, open(f'{OUT}_summary.json', 'w'), indent=1); print(s)
