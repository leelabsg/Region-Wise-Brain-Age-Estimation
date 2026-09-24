#!/usr/bin/env python3
"""
Manuscript restructure -- STEP A: mapping-INDEPENDENT audit + drafts only.
Preflight+audit, retained-set reconstruction, Step3 exact pair-set audit, region
signal-flow, ST10 draft, legacy 16-group inventory, empty 11-group template, manifest.
No GWAS/MAGMA/LDSC/lambda rerun; no modify/overwrite confirmed inputs; no 16->11 merge;
no harmonization map / Fig4b / mapping-dependent ST11 / prose.
Analysis key = (region, Dataset, Cell_type). Negation-free retained logic.
"""
import os, re, sys, json, ast, hashlib, subprocess
from datetime import datetime, timezone
import numpy as np
import pandas as pd

def _env(name):
    import os as _os, sys as _sys
    v = _os.environ.get(name)
    if not v:
        _sys.exit(f"ERROR: {name} is not set (see gwas/config/paths.env.example)")
    return v


# ------------------------------------------------------------------ constants
REPO = _env("MAGMA_ROOT")
OUT  = os.path.join(REPO, "outputs/manuscript_tables")
RUNS = os.path.join(OUT, "runs")
A_PATH   = os.path.join(OUT, "tableA_step1_marginal_all_regions.tsv")
SUM_PATH = os.path.join(OUT, "supp_step1_2_summary_all.tsv")
S3_PATH  = os.path.join(OUT, "supp_step3_cross_dataset_all.tsv")
B_PATH   = os.path.join(OUT, "tableB_independent_significant_set.tsv")
LEGACY_MAP         = os.path.join(REPO, "outputs/combined/cellgroup_mapping.tsv")
LEGACY_RULE_SOURCE = os.path.join(REPO, "src/regional_pipeline/cellgroup_map.py")
CONFIG_PATH        = os.path.join(OUT, "harmonization_config_v2_APPROVED.tsv")
APPROVED_MAP_PATH  = os.path.join(OUT, "harmonization_map_v2_APPROVED.tsv")

# Derived from the data in main(), not pinned. A value that happens to match one run is
# still a hard-coded constant the next.
M_BONF = None
RAW_THRESHOLD = None
EXPECTED_REGIONS = ["global","caudate","frontal_lobe","insula","occipital_lobe",
                    "parietal_lobe","putamen","temporal_lobe","thalamus","cerebellum"]
B_STATE_ALLOW      = {"dropped_at_step2","entered_step3","step2_retained","step1_single"}
B_RETAINED_STATES  = {"entered_step3","step2_retained","step1_single"}
SINGLETON_ALLOWED  = {"step1_single","entered_step3"}
MULTI_ALLOWED      = {"dropped_at_step2","entered_step3","step2_retained"}
COND_STATE_ALLOW   = {"single","joint","partial-joint","colinear",
                      "joint-drop","colinear-drop","partial-joint-drop"}
REQ_COLS = {
    "A":  ["region","Dataset","Cell_type","NGENES","BETA","BETA_STD","SE","P","P.adj.pds","P.adj"],
    "SUM":["region","Dataset","Cell_type","cond_state","cond_cell_type","step3"],
    "S3": ["region","MODEL","Dataset","Cell_type"],
    "B":  ["region","Dataset","Cell_type","survived_at"],
    "LEGACY_MAP":["label","group","rule_id"],
}
# fixed schemas so 0-row conflict TSVs still carry headers
OFFENDING_COLUMNS   = ["invariant","reason","region","Dataset","Cell_type","MODEL","observed","expected"]
MISSING_STEP2_COLS  = ["region","Dataset","Cell_type","in_SUM","in_B"]
STATE_CONFLICT_COLS = ["region","Dataset","Cell_type","B_state","SUM_step3"]
STEP3_CONFLICT_COLS = ["relation_type","region","MODEL","Dataset","Cell_type","pair","rows"]
# provenance / dataset-specific review tokens (single source for risk_flag AND dset_specific)
PROV_TOKENS = ["_as_","_per_","_peric_","_h_bg_","_sox6_","_matrix_","imolgs","_rb_","_hba"]
# option-(a): source-reported NA handling. STAT_COLS are the six statistic columns MAGMA emits.
STAT_COLS = ["BETA","BETA_STD","SE","P","P.adj.pds","P.adj"]
SOURCE_NA_COLS = ["region","Dataset","Cell_type","NGENES","BETA_raw","BETA_STD_raw","SE_raw",
                  "P_raw","P_adj_pds_raw","P_adj_raw","step1_test_status"]

# ------------------------------------------------------------------ helpers
def sha256(path):
    if not path or not os.path.exists(path): return None
    h = hashlib.sha256()
    with open(path,"rb") as f:
        for c in iter(lambda:f.read(1<<20), b""): h.update(c)
    return h.hexdigest()

def read_str(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

def empty_key_mask(df, cols):
    m = pd.Series(False, index=df.index)
    for c in cols:
        m = m | df[c].astype(str).str.strip().eq("")
    return m

def now_utc():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

def git_hash():
    try:
        return subprocess.run(["git","-C",REPO,"rev-parse","HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or "unavailable"
    except Exception:
        return "unavailable"

def write_manifest(run, inputs, status, halt_reason="", extra=None):
    man = dict(utc_timestamp=now_utc(), status=status, run_dir=run,
               script_path=os.path.abspath(__file__), script_sha256=sha256(os.path.abspath(__file__)),
               git_hash=git_hash(), python=sys.version.split()[0], pandas=pd.__version__, numpy=np.__version__,
               inputs={k:{"path":v,"sha256":sha256(v),"exists":os.path.exists(v)} for k,v in inputs.items()},
               outputs={fn:sha256(os.path.join(run,fn)) for fn in sorted(os.listdir(run))},
               M_BONF=M_BONF, raw_threshold=RAW_THRESHOLD, expected_regions=EXPECTED_REGIONS,
               B_state_allowlist=sorted(B_STATE_ALLOW), cond_state_allowlist=sorted(COND_STATE_ALLOW),
               state_compatibility=dict(singleton=sorted(SINGLETON_ALLOWED), multi=sorted(MULTI_ALLOWED)),
               config_path_present=os.path.exists(CONFIG_PATH), approved_map_present=os.path.exists(APPROVED_MAP_PATH),
               user_11group_config_present=os.path.exists(CONFIG_PATH),
               mapping_dependent_outputs="HALTED", halt_reason=halt_reason)
    if extra: man.update(extra)
    with open(os.path.join(run,"run_manifest.json"),"w") as fh:
        json.dump(man, fh, indent=2, default=str)

def match_type(pat):
    if pat.startswith("^_"): return "ANCHORED_PREFIX"
    toks = re.findall(r"_[a-z0-9]+_", pat)
    phrase = bool(re.search(r"[a-z]{5,}", re.sub(r"_[a-z0-9]{1,4}_","",pat)))
    return "TOKEN" if (toks and not phrase) else ("REGEX_MIXED" if (toks and phrase) else "SUBSTRING")

def risk_flag(pat):
    return ";".join(f"PROV_OR_SHORT:{t}" for t in PROV_TOKENS if t in pat)

def dset_specific(pat):
    return any(t in pat for t in PROV_TOKENS)

# ================================================================== main
def main():
    os.makedirs(RUNS, exist_ok=True)
    base = f"restructure_{now_utc()}"
    run = os.path.join(RUNS, base); i = 1
    while os.path.exists(run):
        run = os.path.join(RUNS, f"{base}_{i}"); i += 1
    os.makedirs(run)
    inputs = {"A":A_PATH,"SUM":SUM_PATH,"S3":S3_PATH,"B":B_PATH,
              "LEGACY_MAP":LEGACY_MAP,"LEGACY_RULE_SOURCE":LEGACY_RULE_SOURCE}

    # -------------------- preflight --------------------
    pre_rows = []; pre_ok = True
    def preflight_one(name, path, required):
        nonlocal pre_ok
        exists = os.path.exists(path); n_rows=""; cols=""; missing=""
        if not exists:
            pre_ok = False
        else:
            try:
                df = read_str(path); n_rows=len(df); cols=",".join(df.columns)
                miss = [c for c in required if c not in df.columns]; missing=",".join(miss)
                if miss: pre_ok = False
            except Exception as e:
                missing = f"READ_ERROR:{e}"; pre_ok = False
        pre_rows.append(dict(input=name, path=path, exists=exists, n_rows=n_rows, columns=cols,
                             missing_required_columns=missing, sha256=sha256(path),
                             status=("OK" if (exists and not missing) else "FAIL")))
    for nm,pth in [("A",A_PATH),("SUM",SUM_PATH),("S3",S3_PATH),("B",B_PATH),("LEGACY_MAP",LEGACY_MAP)]:
        preflight_one(nm, pth, REQ_COLS[nm])
    rs_exists = os.path.exists(LEGACY_RULE_SOURCE)
    if not rs_exists: pre_ok = False
    pre_rows.append(dict(input="LEGACY_RULE_SOURCE", path=LEGACY_RULE_SOURCE, exists=rs_exists, n_rows="",
                         columns="", missing_required_columns=("" if rs_exists else "FILE_MISSING"),
                         sha256=sha256(LEGACY_RULE_SOURCE), status=("OK" if rs_exists else "FAIL")))
    pd.DataFrame(pre_rows).to_csv(os.path.join(run,"preflight_report.tsv"), sep="\t", index=False)
    if not pre_ok:
        write_manifest(run, inputs, "PRECHECK_FAILED", "missing input file or required column")
        print("PREFLIGHT FAILED ->", os.path.join(run,"preflight_report.tsv"))
        print(pd.DataFrame(pre_rows)[["input","exists","missing_required_columns","status"]].to_string(index=False))
        sys.exit(1)

    # -------------------- load --------------------
    A_s = read_str(A_PATH); SUM_s = read_str(SUM_PATH); S3_s = read_str(S3_PATH); B_s = read_str(B_PATH)
    A = A_s.copy()
    # ---- derive the pooled test count from THIS run -------------------------------
    global M_BONF, RAW_THRESHOLD
    _npr = A.groupby("region").size()
    if _npr.nunique() != 1:
        print("FAIL: regions do not have equal step-1 row counts:", _npr.to_dict()); sys.exit(1)
    M_BONF = int(_npr.iloc[0])
    RAW_THRESHOLD = 0.05 / M_BONF
    print(f"[derived] rows per region={M_BONF}  regions={_npr.size}  total={len(A)}  "
          f"raw Bonferroni threshold={RAW_THRESHOLD:.6e}")
    for c in ["P","P.adj","P.adj.pds"]:
        A[c+"_num"] = pd.to_numeric(A[c], errors="coerce")
    A["NGENES_num"] = pd.to_numeric(A["NGENES"], errors="coerce")
    SUM = SUM_s.copy(); SUM["step3_num"] = pd.to_numeric(SUM["step3"], errors="coerce")

    # ---- separate raw source-missing tokens ("" / case-insensitive "NA") from malformed non-numeric ----
    def _raw_missing(col):
        s = A_s[col].astype(str).str.strip()
        return s.eq("") | s.str.upper().eq("NA")
    raw_missing   = {c: _raw_missing(c) for c in STAT_COLS}
    stat_num      = {c: pd.to_numeric(A_s[c], errors="coerce") for c in STAT_COLS}
    malformed     = {c: (stat_num[c].isna() & ~raw_missing[c]) for c in STAT_COLS}  # NaN but NOT '' / 'NA'
    malformed_any = pd.concat([malformed[c] for c in STAT_COLS], axis=1).any(axis=1)
    n_missing_per_row     = pd.concat([raw_missing[c] for c in STAT_COLS], axis=1).sum(axis=1)
    all_stats_missing     = n_missing_per_row.eq(len(STAT_COLS))
    any_stats_missing     = n_missing_per_row.ge(1)
    partial_stats_missing = any_stats_missing & ~all_stats_missing
    all_stats_numeric     = pd.concat([stat_num[c].notna() for c in STAT_COLS], axis=1).all(axis=1)
    ngenes_ok             = A["NGENES_num"].notna()
    numeric_stat_pattern  = all_stats_numeric & n_missing_per_row.eq(0)
    tested_ngenes_bad     = numeric_stat_pattern & ~ngenes_ok   # numeric stats but NGENES missing/malformed
    tested_row            = numeric_stat_pattern & ngenes_ok
    source_reported_na    = all_stats_missing & ngenes_ok
    tested_by_region = tested_row.groupby(A["region"]).sum().reindex(EXPECTED_REGIONS).fillna(0).astype(int).to_dict()
    sna_by_region    = source_reported_na.groupby(A["region"]).sum().reindex(EXPECTED_REGIONS).fillna(0).astype(int).to_dict()

    inv=[]; offend=[]; missing_step2=[]; state_conf=[]; step3_conf=[]; blocking={"v":False}
    def check(name, passed, detail=""):
        inv.append(dict(invariant=name, result="PASS" if passed else "FAIL", detail=str(detail)))
        if not passed: blocking["v"] = True
    def add_offend(invariant, reason, region="", Dataset="", Cell_type="", MODEL="", observed="", expected=""):
        offend.append(dict(invariant=invariant, reason=reason, region=region, Dataset=Dataset,
                           Cell_type=Cell_type, MODEL=MODEL, observed=str(observed), expected=str(expected)))

    # -------------------- Step1 audit --------------------
    check("regions_exactly_expected10", set(A.region)==set(EXPECTED_REGIONS) and A.region.nunique()==10, f"got={sorted(set(A.region))}")
    nt = A.groupby("region").size()
    # equal across regions, whatever the value: catches every region being wrong by the
    # same amount, which pinning a single number does not.
    check("each_region_equal_row_count", nt.nunique()==1, nt.to_dict())
    check("total_rows_eq_regions_x_per_region", len(A)==len(EXPECTED_REGIONS)*M_BONF,
          f"rows={len(A)} regions={len(EXPECTED_REGIONS)} per_region={M_BONF}")
    for nm,df in [("A",A_s),("SUM",SUM_s),("B",B_s)]:
        em = empty_key_mask(df,["region","Dataset","Cell_type"]); ne=int(em.sum())
        check(f"{nm}_key_no_empty", ne==0, f"{ne} empty-key rows")
        for _,r in df[em].iterrows(): add_offend(f"{nm}_key_no_empty","empty_key",r.region,r.Dataset,r.Cell_type)
        dk = df.duplicated(subset=["region","Dataset","Cell_type"], keep=False)
        nd = int(df.duplicated(subset=["region","Dataset","Cell_type"]).sum())
        check(f"{nm}_key_unique", nd==0, f"{nd} duplicate keys")
        for _,r in df[dk].iterrows(): add_offend(f"{nm}_key_unique","duplicate_key",r.region,r.Dataset,r.Cell_type)
    # ---- malformed STAT_COLS (non-numeric AND not a recognized '' / 'NA' token) -> blocking ----
    for c in STAT_COLS:
        nbad = int(malformed[c].sum()); cc = c.replace(" ","")
        check(f"A_{cc}_malformed_non_numeric", nbad==0, f"{nbad} malformed (not '' / 'NA')")
        for _,r in A[malformed[c]].iterrows():
            add_offend(f"A_{cc}_malformed_non_numeric","malformed_non_numeric",r.region,r.Dataset,r.Cell_type,"",r[c],"numeric or recognized-NA")
    # ---- partial (all-or-none) STAT_COLS missing within a row -> blocking ----
    check("A_partial_stat_missing", int(partial_stats_missing.sum())==0, f"{int(partial_stats_missing.sum())} partial-missing rows")
    for _,r in A[partial_stats_missing].iterrows():
        add_offend("A_partial_stat_missing","partial_stat_missing",r.region,r.Dataset,r.Cell_type,"",
                   ";".join(f"{c}={r[c]}" for c in STAT_COLS),"all-or-none")
    # ---- NGENES must be present+numeric for BOTH statuses ----
    check("A_tested_ngenes_present", int(tested_ngenes_bad.sum())==0, f"{int(tested_ngenes_bad.sum())} numeric-stat rows lacking numeric NGENES")
    for _,r in A[tested_ngenes_bad].iterrows():
        add_offend("A_tested_ngenes_present","tested_row_ngenes_missing_or_non_numeric",r.region,r.Dataset,r.Cell_type,"",r["NGENES"],"numeric NGENES")
    allmiss_ngenes_bad = all_stats_missing & ~ngenes_ok
    check("A_source_na_ngenes_present", int(allmiss_ngenes_bad.sum())==0, f"{int(allmiss_ngenes_bad.sum())} all-stats-missing rows lacking numeric NGENES")
    for _,r in A[allmiss_ngenes_bad].iterrows():
        add_offend("A_source_na_ngenes_present","ngenes_missing",r.region,r.Dataset,r.Cell_type,"",r["NGENES"],"numeric NGENES")
    # ---- every row is EXACTLY one of TESTED / SOURCE_REPORTED_NA (exhaustive + mutually exclusive) ----
    unclassified_row = ~(tested_row | source_reported_na)
    overlap_row      = tested_row & source_reported_na
    check("A_step1_test_status_exhaustive", int(unclassified_row.sum())==0, f"{int(unclassified_row.sum())} rows not classified as TESTED or SOURCE_REPORTED_NA")
    check("A_step1_test_status_mutually_exclusive", int(overlap_row.sum())==0, f"{int(overlap_row.sum())} rows assigned to both statuses")
    for _,r in A[unclassified_row].iterrows():
        add_offend("A_step1_test_status_exhaustive","unclassified_step1_row",r.region,r.Dataset,r.Cell_type,"",
                   ";".join(f"{c}={r[c]}" for c in ["NGENES"]+STAT_COLS),"exactly one of TESTED or SOURCE_REPORTED_NA")
    # ---- source-reported NA: all six STAT_COLS missing AND NGENES present+numeric ----
    # The invariant is STRUCTURAL: the same (Dataset, Cell_type) keys are untestable in every
    # region, so rows == keys x regions and every region carries exactly `keys` of them. Tying the
    # three counts to each other is stronger than pinning them to a previous run's values.
    sna = A[source_reported_na]
    sna_full_keys = set(zip(sna.region, sna.Dataset, sna.Cell_type))
    sna_dc_keys   = set(zip(sna.Dataset, sna.Cell_type))
    regions_per_key = {k:set() for k in sna_dc_keys}
    for rg,d,c in sna_full_keys: regions_per_key[(d,c)].add(rg)
    n_keys             = len(sna_dc_keys)
    per_region_ok      = all(int((sna.region==rg).sum())==n_keys for rg in EXPECTED_REGIONS)
    each_key_10        = all(regions_per_key[(d,c)]==set(EXPECTED_REGIONS) for (d,c) in sna_dc_keys)
    keyset_by_region   = {rg:set(zip(sna[sna.region==rg].Dataset, sna[sna.region==rg].Cell_type)) for rg in EXPECTED_REGIONS}
    cross_region_equal = all(keyset_by_region[rg]==sna_dc_keys for rg in EXPECTED_REGIONS)
    pattern_ok = (int(source_reported_na.sum())==n_keys*len(EXPECTED_REGIONS)
                  and set(sna.region)==set(EXPECTED_REGIONS) and per_region_ok and each_key_10 and cross_region_equal)
    check("A_source_reported_na_pattern", pattern_ok,
          f"rows={int(source_reported_na.sum())}(=keys x regions? {int(source_reported_na.sum())==n_keys*len(EXPECTED_REGIONS)}) keys={n_keys} "
          f"regions={sna.region.nunique()} per_region_ok={per_region_ok} each_key_10={each_key_10} cross_region_equal={cross_region_equal}")
    if not pattern_ok:
        add_offend("A_source_reported_na_pattern","source_na_pattern_summary_mismatch",observed=(
                   f"rows={int(source_reported_na.sum())};keys={len(sna_dc_keys)};regions={sna.region.nunique()};"
                   f"per_region_ok={per_region_ok};each_key_10={each_key_10};cross_region_equal={cross_region_equal}"),
                   expected=f"rows=keys*regions={n_keys*len(EXPECTED_REGIONS)};per_region=keys={n_keys}")
        for _,r in sna.iterrows():
            add_offend("A_source_reported_na_pattern","source_na_row_context",r.region,r.Dataset,r.Cell_type,"",
                       ";".join(f"{c}={r[c]}" for c in ["NGENES"]+STAT_COLS),"fixed 10-keys x 10-regions pattern")
    # ---- numeric range + significance/discrepancy on TESTED rows only (no impute of source-NA) ----
    for c in ["P","P.adj","P.adj.pds"]:
        rng_mask = tested_row & ((A[c+"_num"]<0)|(A[c+"_num"]>1)); nrng=int(rng_mask.sum())
        check(f"A_{c}_in_0_1_tested", nrng==0, f"{nrng} tested rows out of [0,1]")
        for _,r in A[rng_mask].iterrows(): add_offend(f"A_{c}_in_0_1_tested","out_of_range",r.region,r.Dataset,r.Cell_type,"",r[c],"[0,1]")
    sig_padj = tested_row & (A["P.adj_num"]<0.05)
    sig_raw  = tested_row & (A["P_num"]<RAW_THRESHOLD)
    mm = sig_padj.ne(sig_raw)
    check("significance_padj_eq_rawthreshold", int(mm.sum())==0, f"{int(mm.sum())} classification mismatches (tested rows)")
    for _,r in A[mm].iterrows():
        add_offend("significance_padj_eq_rawthreshold","class_mismatch",r.region,r.Dataset,r.Cell_type,"",
                   f"P={r.P};Padj={r['P.adj']}", f"(P<{RAW_THRESHOLD:.6e})={bool(r['P_num']<RAW_THRESHOLD)}")
    disc_t = (A.loc[tested_row,"P.adj_num"]-np.minimum(A.loc[tested_row,"P_num"]*M_BONF,1.0)).abs()
    disc_stats = dict(max=float(disc_t.max()), median=float(disc_t.median()),
                      p95=float(np.nanpercentile(disc_t,95)), n_class_mismatch=int(mm.sum()), n_computed_step1_rows=int(tested_row.sum()))
    top100 = (A.loc[tested_row].assign(abs_discrepancy=disc_t.values)
                .sort_values("abs_discrepancy", ascending=False)
                [["region","Dataset","Cell_type","P","P.adj","abs_discrepancy"]].head(100))
    top100.to_csv(os.path.join(run,"padj_formula_discrepancy_top100.tsv"), sep="\t", index=False)
    inv.append(dict(invariant="padj_formula_discrepancy", result="DIAGNOSTIC",
                    detail=f"max={disc_stats['max']:.3e} median={disc_stats['median']:.3e} p95={disc_stats['p95']:.3e} n_computed_step1_rows={disc_stats['n_computed_step1_rows']} "
                           f"(descriptive only; classification equality is the blocking criterion)"))
    sigA = A[sig_padj].copy()
    sigc = sigA.groupby("region").size().reindex(EXPECTED_REGIONS).fillna(0).astype(int).to_dict()
    # Observed, not compared to a previous run's counts.
    inv.append(dict(invariant="sig_counts_observed", result="DIAGNOSTIC",
                    detail="; ".join(f"{r}={sigc[r]}" for r in EXPECTED_REGIONS) + f"; total={sum(sigc.values())}"))
    check("sig_counts_sum_consistent", sum(sigc.values())==int(sig_padj.sum()),
          f"per_region_sum={sum(sigc.values())} mask_sum={int(sig_padj.sum())}")

    # -------------------- retained set --------------------
    sig_keys = set(zip(sigA.region,sigA.Dataset,sigA.Cell_type))
    SUM_keys = set(zip(SUM.region,SUM.Dataset,SUM.Cell_type)); B_keys = set(zip(B_s.region,B_s.Dataset,B_s.Cell_type))
    check("B_keys_eq_sig_keys", B_keys==sig_keys, f"B\\sig={len(B_keys-sig_keys)} sig\\B={len(sig_keys-B_keys)}")
    for k in sorted(B_keys-sig_keys): add_offend("B_keys_eq_sig_keys","B_not_in_sig",k[0],k[1],k[2])
    for k in sorted(sig_keys-B_keys): add_offend("B_keys_eq_sig_keys","sig_not_in_B",k[0],k[1],k[2])
    check("SUM_subset_sig", SUM_keys.issubset(sig_keys), f"{len(SUM_keys-sig_keys)} SUM not in sig")
    for k in sorted(SUM_keys-sig_keys): add_offend("SUM_subset_sig","SUM_nonsignificant_extra",k[0],k[1],k[2])
    check("B_subset_sig", B_keys.issubset(sig_keys), f"{len(B_keys-sig_keys)} B not in sig")

    per_rd = sigA.groupby(["region","Dataset"]).size().rename("nsig").reset_index()
    sigA = sigA.merge(per_rd, on=["region","Dataset"])
    singleton_keys = set(map(tuple, sigA.loc[sigA.nsig==1,["region","Dataset","Cell_type"]].values))
    multi_keys     = set(map(tuple, sigA.loc[sigA.nsig>=2,["region","Dataset","Cell_type"]].values))
    for k in sorted(multi_keys):
        inS = k in SUM_keys; inB = k in B_keys
        if not (inS and inB):
            missing_step2.append(dict(region=k[0],Dataset=k[1],Cell_type=k[2],in_SUM=inS,in_B=inB))
            add_offend("multi_keys_have_SUM_and_B","missing_SUM_or_B",k[0],k[1],k[2],"",f"in_SUM={inS};in_B={inB}","present")
    check("multi_keys_have_SUM_and_B", len(missing_step2)==0, f"{len(missing_step2)} missing")

    Bstate = {(r,d,c):s for r,d,c,s in zip(B_s.region,B_s.Dataset,B_s.Cell_type,B_s.survived_at)}
    bad_state = [(k,v) for k,v in Bstate.items() if v not in B_STATE_ALLOW or str(v).strip()==""]
    check("B_states_in_allowlist", len(bad_state)==0, f"{len(bad_state)} unknown/blank states")
    for k,v in bad_state: add_offend("B_states_in_allowlist","unknown_or_blank_state",k[0],k[1],k[2],"",v,"|".join(sorted(B_STATE_ALLOW)))
    incompat = []
    for k in singleton_keys:
        s = Bstate.get(k)
        if s is not None and s not in SINGLETON_ALLOWED:
            incompat.append(k); add_offend("state_compatibility","singleton_bad_state",k[0],k[1],k[2],"",s,"/".join(sorted(SINGLETON_ALLOWED)))
    for k in multi_keys:
        s = Bstate.get(k)
        if s is not None and s not in MULTI_ALLOWED:
            incompat.append(k); add_offend("state_compatibility","multi_bad_state",k[0],k[1],k[2],"",s,"/".join(sorted(MULTI_ALLOWED)))
    check("state_compatibility_singleton_multi", len(incompat)==0, f"{len(incompat)} incompatible")
    explicit_step2_retained_keys = {k for k in multi_keys if Bstate.get(k) in {"entered_step3","step2_retained"}}
    check("singleton_step2retained_disjoint", len(singleton_keys & explicit_step2_retained_keys)==0,
          f"overlap={len(singleton_keys & explicit_step2_retained_keys)}")
    retained_signal_set = singleton_keys | explicit_step2_retained_keys
    check("retained_subset_sig", retained_signal_set.issubset(sig_keys), f"{len(retained_signal_set-sig_keys)} not sig")
    B_retained_positive = {k for k,s in Bstate.items() if s in B_RETAINED_STATES}
    check("retained_eq_B_positive_states", retained_signal_set==B_retained_positive, f"symdiff={len(retained_signal_set ^ B_retained_positive)}")

    # -------------------- SUM cross-check --------------------
    dupSUM = int(SUM.duplicated(subset=["region","Dataset","Cell_type"]).sum())
    check("SUM_key_unique_precheck", dupSUM==0, f"{dupSUM} SUM dups")
    bad_step3_mask = SUM["step3_num"].isna() | ~SUM["step3_num"].isin([0,1])
    check("SUM_step3_binary", int(bad_step3_mask.sum())==0,
          f"{int(bad_step3_mask.sum())} missing/non-binary values; raw_values={sorted(SUM['step3'].astype(str).unique())}")
    for _,row in SUM[bad_step3_mask].iterrows():
        add_offend("SUM_step3_binary","missing_non_numeric_or_non_binary_step3",
                   row["region"],row["Dataset"],row["Cell_type"],"",row["step3"],"0 or 1")
    cond_counts = SUM["cond_state"].value_counts(dropna=False).to_dict()
    bad_cond = []
    for r,d,c,cs in zip(SUM.region,SUM.Dataset,SUM.Cell_type,SUM.cond_state):
        raw = str(cs)
        if raw.strip()=="":
            bad_cond.append((r,d,c,raw,"BLANK")); continue
        if raw=="NA":
            bad_cond.append((r,d,c,raw,"LITERAL_NA")); continue
        if any(tok.strip() not in COND_STATE_ALLOW for tok in raw.split(";")):
            bad_cond.append((r,d,c,raw,"UNKNOWN_TOKEN"))
    check("cond_state_in_allowlist", len(bad_cond)==0, f"{len(bad_cond)} unknown/blank/NA cond_state")
    for r,d,c,raw,why in bad_cond: add_offend("cond_state_in_allowlist",why,r,d,c,"",raw,"|".join(sorted(COND_STATE_ALLOW)))
    SUMstep3 = {(r,d,c):v for r,d,c,v in zip(SUM.region,SUM.Dataset,SUM.Cell_type,SUM["step3_num"])}
    for k in multi_keys:
        bs = Bstate.get(k); s3 = SUMstep3.get(k)
        if bs is None: continue
        if s3 is None or pd.isna(s3):
            state_conf.append(dict(region=k[0],Dataset=k[1],Cell_type=k[2],B_state=bs,SUM_step3=""))
            add_offend("B_SUM_retained_agree","missing_or_non_numeric_SUM_step3",k[0],k[1],k[2],"",
                       f"B={bs};SUM.step3={s3}","SUM.step3 in {0,1}")
            continue
        b_retained = bs in {"entered_step3","step2_retained"}; sum_retained = s3==1
        if b_retained != sum_retained:
            state_conf.append(dict(region=k[0],Dataset=k[1],Cell_type=k[2],B_state=bs,SUM_step3=int(s3)))
            add_offend("B_SUM_retained_agree","B_vs_SUM_mismatch",k[0],k[1],k[2],"",f"B={bs}",f"SUM.step3={int(s3)}")
    check("B_SUM_retained_agree", len(state_conf)==0, f"{len(state_conf)} conflicts")
    sum_single_keys = {(r,d,c) for r,d,c,cs in zip(SUM.region,SUM.Dataset,SUM.Cell_type,SUM.cond_state) if cs=="single"}
    singleton_in_sum = singleton_keys & SUM_keys
    check("cond_single_eq_singleton_in_SUM", sum_single_keys==singleton_in_sum, f"symdiff={len(sum_single_keys ^ singleton_in_sum)}")
    for k in sorted(sum_single_keys-singleton_in_sum): add_offend("cond_single_eq_singleton_in_SUM","single_not_singleton",k[0],k[1],k[2],"","cond=single","Table-A singleton")
    for k in sorted(singleton_in_sum-sum_single_keys): add_offend("cond_single_eq_singleton_in_SUM","singleton_not_cond_single",k[0],k[1],k[2],"","Table-A singleton","cond=single")
    multi_single = multi_keys & sum_single_keys
    check("multi_not_cond_single", len(multi_single)==0, f"{len(multi_single)} multi are cond=single")
    for k in sorted(multi_single): add_offend("multi_not_cond_single","multi_is_cond_single",k[0],k[1],k[2],"","cond=single","not single")

    # -------------------- Step3 exact pair-set audit --------------------
    s3km = empty_key_mask(S3_s,["region","MODEL","Dataset","Cell_type"])
    check("S3_key_no_empty", int(s3km.sum())==0, f"{int(s3km.sum())} empty S3 keys")
    for _,r in S3_s[s3km].iterrows(): add_offend("S3_key_no_empty","MISSING_S3_KEY",r.region,r.Dataset,r.Cell_type,r.MODEL,"empty","non-empty")
    unexpected_regions = sorted(set(S3_s.region)-set(EXPECTED_REGIONS))
    check("S3_no_unexpected_region", len(unexpected_regions)==0, unexpected_regions)
    for rr in unexpected_regions: add_offend("S3_no_unexpected_region","UNEXPECTED_S3_REGION",rr,"","","","present","expected10")
    n_abnormal=n_self=n_within=0
    observed_pairs = {r:set() for r in EXPECTED_REGIONS}; observed_signals = {r:set() for r in EXPECTED_REGIONS}
    pair_model_map = {}
    for (r,m),g in S3_s.groupby(["region","MODEL"]):
        if r not in EXPECTED_REGIONS: continue
        if len(g)!=2:
            n_abnormal+=1; step3_conf.append(dict(relation_type="ABNORMAL_MODEL_ROWS",region=r,MODEL=m,Dataset="",Cell_type="",pair="",rows=len(g)))
            add_offend("step3_no_abnormal_model","ABNORMAL_MODEL_ROWS",r,"","",m,len(g),2); continue
        eps = list(zip(g.Dataset,g.Cell_type)); a,b = eps
        if a==b:
            n_self+=1; step3_conf.append(dict(relation_type="SELF_PAIR",region=r,MODEL=m,Dataset=a[0],Cell_type=a[1],pair="",rows=2))
            add_offend("step3_no_self_pair","SELF_PAIR",r,a[0],a[1],m,"self","cross"); continue
        if a[0]==b[0]:
            n_within+=1; step3_conf.append(dict(relation_type="WITHIN_DATASET_MODEL",region=r,MODEL=m,Dataset=a[0],Cell_type="",pair="",rows=2))
            add_offend("step3_no_within_dataset_model","WITHIN_DATASET_MODEL",r,a[0],"",m,"same_dataset","cross_dataset"); continue
        key = frozenset(eps); pair_model_map.setdefault((r,key),[]).append(m)
        observed_pairs[r].add(key); observed_signals[r].add(a); observed_signals[r].add(b)
    dup_pairs = {k:v for k,v in pair_model_map.items() if len(v)>1}
    check("step3_no_abnormal_model", n_abnormal==0, f"{n_abnormal} models !=2 rows")
    check("step3_no_self_pair", n_self==0, f"{n_self} self-pairs")
    check("step3_no_within_dataset_model", n_within==0, f"{n_within} within-dataset models")
    check("step3_no_duplicate_pair_multimodel", len(dup_pairs)==0, f"{len(dup_pairs)} dup pairs")
    for (r,key),ms in dup_pairs.items():
        step3_conf.append(dict(relation_type="DUPLICATE_PAIR_MODEL",region=r,MODEL=";".join(map(str,ms)),Dataset="",Cell_type="",
                               pair=";".join(sorted(f"{d}|{c}" for d,c in key)),rows=len(ms)))
        add_offend("step3_no_duplicate_pair_multimodel","DUPLICATE_PAIR_MODEL",r,"","",";".join(map(str,ms)),len(ms),1)

    entered_signals = {r:{(d,c) for (rr,d,c),s in Bstate.items() if rr==r and s=="entered_step3"} for r in EXPECTED_REGIONS}
    expected_pairs = {}; expected_pair_participants = {}
    for r in EXPECTED_REGIONS:
        ent = sorted(entered_signals[r]); ep=set(); part=set()
        for a_i in range(len(ent)):
            for b_i in range(a_i+1,len(ent)):
                if ent[a_i][0]!=ent[b_i][0]:
                    ep.add(frozenset([ent[a_i],ent[b_i]])); part.add(ent[a_i]); part.add(ent[b_i])
        expected_pairs[r]=ep; expected_pair_participants[r]=part
    eq_sig_part = all(entered_signals[r]==expected_pair_participants[r] for r in EXPECTED_REGIONS)
    eq_sig_obs  = all(entered_signals[r]==observed_signals[r] for r in EXPECTED_REGIONS)
    eq_pairs    = all(expected_pairs[r]==observed_pairs[r] for r in EXPECTED_REGIONS)
    check("step3_entered_eq_pair_participants", eq_sig_part, "ok" if eq_sig_part else "entered_step3 signal without valid cross-dataset pair")
    check("step3_entered_eq_observed_signals", eq_sig_obs, "ok" if eq_sig_obs else "entered_step3 != observed step3 signals")
    check("step3_expected_eq_observed_pairs", eq_pairs, "ok" if eq_pairs else "expected != observed cross-dataset pairs")
    for r in EXPECTED_REGIONS:
        for sg in entered_signals[r]-expected_pair_participants[r]:
            step3_conf.append(dict(relation_type="ENTERED_SIGNAL_WITHOUT_PAIR",region=r,MODEL="",Dataset=sg[0],Cell_type=sg[1],pair="",rows=""))
            add_offend("step3_entered_eq_pair_participants","ENTERED_SIGNAL_WITHOUT_PAIR",r,sg[0],sg[1],"","no_pair","paired")
        for sg in entered_signals[r]-observed_signals[r]:
            step3_conf.append(dict(relation_type="MISSING_EXPECTED_SIGNAL",region=r,MODEL="",Dataset=sg[0],Cell_type=sg[1],pair="",rows=""))
            add_offend("step3_entered_eq_observed_signals","MISSING_EXPECTED_SIGNAL",r,sg[0],sg[1],"","absent","present")
        for sg in observed_signals[r]-entered_signals[r]:
            step3_conf.append(dict(relation_type="UNEXPECTED_OUTPUT_SIGNAL",region=r,MODEL="",Dataset=sg[0],Cell_type=sg[1],pair="",rows=""))
            add_offend("step3_entered_eq_observed_signals","UNEXPECTED_OUTPUT_SIGNAL",r,sg[0],sg[1],"","present","absent")
        for pair in expected_pairs[r]-observed_pairs[r]:
            ps=";".join(sorted(f"{d}|{c}" for d,c in pair))
            step3_conf.append(dict(relation_type="MISSING_EXPECTED_PAIR",region=r,MODEL="",Dataset="",Cell_type="",pair=ps,rows=""))
            add_offend("step3_expected_eq_observed_pairs","MISSING_EXPECTED_PAIR",r,"","","","absent",ps)
        for pair in observed_pairs[r]-expected_pairs[r]:
            ps=";".join(sorted(f"{d}|{c}" for d,c in pair))
            step3_conf.append(dict(relation_type="UNEXPECTED_OUTPUT_PAIR",region=r,MODEL="",Dataset="",Cell_type="",pair=ps,rows=""))
            add_offend("step3_expected_eq_observed_pairs","UNEXPECTED_OUTPUT_PAIR",r,"","","",ps,"absent")
    tot_rows = len(S3_s); tot_pairs = sum(len(observed_pairs[r]) for r in EXPECTED_REGIONS)
    tot_sigs = sum(len(observed_signals[r]) for r in EXPECTED_REGIONS)
    inv.append(dict(invariant="step3_xcheck_diagnostic", result="DIAGNOSTIC",
                    detail=f"raw_rows={tot_rows} pairs={tot_pairs} signals={tot_sigs}"))

    # -------------------- source-reported NA disjointness + audit summary --------------------
    # S3 overlap uses ALL raw S3 endpoints (observed_signals skips abnormal/self/within rows, so it is too weak).
    s3_all_signal_keys = set(zip(S3_s["region"], S3_s["Dataset"], S3_s["Cell_type"]))
    ov_sig = sna_full_keys & sig_keys; ov_sum = sna_full_keys & SUM_keys
    ov_b   = sna_full_keys & B_keys;   ov_s3  = sna_full_keys & s3_all_signal_keys
    check("source_na_disjoint_sig",   len(ov_sig)==0, f"{len(ov_sig)} overlap with Step1 significant")
    check("source_na_disjoint_SUM",   len(ov_sum)==0, f"{len(ov_sum)} overlap with SUM")
    check("source_na_disjoint_B",     len(ov_b)==0,   f"{len(ov_b)} overlap with B")
    check("source_na_disjoint_S3_raw",len(ov_s3)==0,  f"{len(ov_s3)} overlap with raw Step3 endpoints")
    for k in sorted(ov_sig): add_offend("source_na_disjoint_sig","source_na_in_sig",k[0],k[1],k[2])
    for k in sorted(ov_sum): add_offend("source_na_disjoint_SUM","source_na_in_SUM",k[0],k[1],k[2])
    for k in sorted(ov_b):   add_offend("source_na_disjoint_B","source_na_in_B",k[0],k[1],k[2])
    for k in sorted(ov_s3):  add_offend("source_na_disjoint_S3_raw","source_na_in_S3",k[0],k[1],k[2])
    sna_per_region = {rg:int((sna.region==rg).sum()) for rg in EXPECTED_REGIONS}
    source_na_audit = dict(rows=int(source_reported_na.sum()), unique_dataset_celltype_keys=len(sna_dc_keys),
        expected_regions=EXPECTED_REGIONS, rows_per_region=sna_per_region,
        regions_per_key=[{"Dataset":d,"Cell_type":c,"regions":sorted(regions_per_key[(d,c)])} for d,c in sorted(sna_dc_keys)],
        partial_missing_rows=int(partial_stats_missing.sum()), malformed_non_numeric_rows=int(malformed_any.sum()),
        overlap_SUM=len(ov_sum), overlap_B=len(ov_b), overlap_S3=len(ov_s3), overlap_sig=len(ov_sig),
        exact_pattern_pass=bool(pattern_ok))

    # -------------------- write audit files (fixed schemas) --------------------
    pd.DataFrame(inv).to_csv(os.path.join(run,"cell_enrichment_input_audit.tsv"), sep="\t", index=False)
    pd.DataFrame(offend, columns=OFFENDING_COLUMNS).to_csv(os.path.join(run,"cell_enrichment_input_audit_offending_rows.tsv"), sep="\t", index=False)
    pd.DataFrame(missing_step2, columns=MISSING_STEP2_COLS).to_csv(os.path.join(run,"cell_enrichment_missing_step2_rows.tsv"), sep="\t", index=False)
    pd.DataFrame(state_conf, columns=STATE_CONFLICT_COLS).to_csv(os.path.join(run,"cell_enrichment_state_conflicts.tsv"), sep="\t", index=False)
    pd.DataFrame(step3_conf, columns=STEP3_CONFLICT_COLS).to_csv(os.path.join(run,"cell_enrichment_step3_pair_conflicts.tsv"), sep="\t", index=False)
    sna_out = A.loc[source_reported_na, ["region","Dataset","Cell_type","NGENES"]].copy()
    for raw,orig in [("BETA_raw","BETA"),("BETA_STD_raw","BETA_STD"),("SE_raw","SE"),
                     ("P_raw","P"),("P_adj_pds_raw","P.adj.pds"),("P_adj_raw","P.adj")]:
        sna_out[raw] = A_s.loc[source_reported_na, orig].values
    sna_out["step1_test_status"] = "SOURCE_REPORTED_NA"
    sna_out[SOURCE_NA_COLS].to_csv(os.path.join(run,"cell_enrichment_source_reported_na_rows.tsv"), sep="\t", index=False)
    with open(os.path.join(run,"cell_enrichment_input_audit.txt"),"w") as fh:
        fh.write("CELL-ENRICHMENT INPUT AUDIT (mapping-independent, Step A)\n"+"="*58+"\n\n")
        fh.write(f"M_BONF={M_BONF}  RAW_THRESHOLD={RAW_THRESHOLD!r}\n")
        fh.write(f"P.adj formula discrepancy (descriptive): {disc_stats}\n")
        fh.write("cond_state observed counts: "+json.dumps(cond_counts)+"\n\n")
        for a in inv: fh.write(f"[{a['result']}] {a['invariant']}: {a['detail']}\n")
        fh.write("\nSOURCE-REPORTED NA AUDIT:\n")
        fh.write(f"  rows={source_na_audit['rows']}  unique_dataset_celltype_keys={source_na_audit['unique_dataset_celltype_keys']}  "
                 f"exact_pattern_pass={source_na_audit['exact_pattern_pass']}\n")
        fh.write("  rows_per_region="+json.dumps(source_na_audit["rows_per_region"], sort_keys=True)+"\n")
        fh.write(f"  partial_missing_rows={source_na_audit['partial_missing_rows']}  "
                 f"malformed_non_numeric_rows={source_na_audit['malformed_non_numeric_rows']}\n")
        fh.write(f"  overlaps: significant={source_na_audit['overlap_sig']}  SUM={source_na_audit['overlap_SUM']}  "
                 f"B={source_na_audit['overlap_B']}  S3_raw={source_na_audit['overlap_S3']}\n")
        fh.write("  verification_note: In the consolidated Table A input, these six statistic fields are blank "
                 "after prior numeric coercion. A separate read-only comparison verified that the corresponding "
                 "original per-region Step 1 files contained literal 'NA' in all six statistic fields while NGENES "
                 "was present. This script does not re-read or re-verify those original per-region files.\n")
        fh.write("  row_file=cell_enrichment_source_reported_na_rows.tsv\n")
        fh.write(f"\nOVERALL: {'ALL PASS' if not blocking['v'] else 'BLOCKING FAILURE'}\n")
    print("== INVARIANTS =="); print(pd.DataFrame(inv).to_string(index=False))
    print(f"\nP.adj formula discrepancy (descriptive): {disc_stats}")
    if offend:
        print(f"\n(offending rows: {len(offend)} total; first 200 shown, full TSV written)")
        print(pd.DataFrame(offend, columns=OFFENDING_COLUMNS).head(200).to_string(index=False))

    if blocking["v"]:
        write_manifest(run, inputs, "AUDIT_FAILED", "one or more blocking invariants failed",
                       extra=dict(source_reported_na_audit=source_na_audit))
        print("\n!! BLOCKING AUDIT FAILURE -- only audit files + manifest written. Downstream HALTED.")
        sys.exit(1)

    # -------------------- region signal-flow --------------------
    multi_by_region = {r:len({k for k in multi_keys if k[0]==r}) for r in EXPECTED_REGIONS}
    rows = []
    for r in EXPECTED_REGIONS:
        sr = sigA[sigA.region==r]
        singl = {k for k in singleton_keys if k[0]==r}; s2ret = {k for k in explicit_step2_retained_keys if k[0]==r}
        rett = singl | s2ret
        rows.append(dict(region=r, n_tested=M_BONF, pooled_test_count=M_BONF,
            n_step1_results_available=tested_by_region[r], n_step1_source_reported_na=sna_by_region[r],
            n_step1_significant=int(sigc[r]),
            n_step1_singletons=len(singl), n_step2_entered=multi_by_region[r], n_step2_retained=len(s2ret),
            n_retained_total=len(rett), n_entered_step3=len(entered_signals[r]),
            n_step3_output_signals=len(observed_signals[r]),
            n_step3_expected_cross_dataset_pairs=len(expected_pairs[r]),
            n_step3_observed_cross_dataset_pairs=len(observed_pairs[r]),
            n_step3_raw_rows=int((S3_s.region==r).sum()),
            n_dataset_files_with_step1_signal=sr.Dataset.nunique(),
            n_dataset_files_with_retained_signal=len({k[1] for k in rett}),
            n_dataset_bases_with_retained_signal=pd.NA, n_resources_with_retained_signal=pd.NA,
            dataset_base_mapping_status="PENDING_AUTHORITATIVE_MAPPING",
            resource_mapping_status="PENDING_AUTHORITATIVE_MAPPING"))
    flow = pd.DataFrame(rows)
    assert (flow.n_retained_total==flow.n_step1_singletons+flow.n_step2_retained).all()
    flow.to_csv(os.path.join(run,"ST11_region_signal_flow_DRAFT.tsv"), sep="\t", index=False)
    print("\n== REGION SIGNAL FLOW =="); print(flow.to_string(index=False))
    cb = flow[flow.region=="cerebellum"].iloc[0]
    print("\n== CEREBELLUM (from data) ==")
    for c in ["n_step1_significant","n_step1_singletons","n_step2_entered","n_step2_retained","n_retained_total","n_entered_step3"]:
        print(f"  {c}={int(cb[c])}")
    print(f"  n_retained_total={int(cb.n_retained_total)} "
          f"= n_step1_singletons({int(cb.n_step1_singletons)}) + n_step2_retained({int(cb.n_step2_retained)})")

    # -------------------- ST10 --------------------
    st10 = A_s[["region","Dataset","Cell_type","NGENES","BETA","BETA_STD","SE","P","P.adj.pds","P.adj"]].copy()
    st10 = st10.rename(columns={"Dataset":"dataset_file","Cell_type":"cell_type_original"})
    st10["step1_test_status"] = np.where(source_reported_na.values, "SOURCE_REPORTED_NA", "TESTED")
    st10["step1_significant"] = sig_padj.astype("Int64").mask(source_reported_na, pd.NA)
    st10["pooled_test_count"] = M_BONF; st10["bonferroni_rawP_threshold"] = RAW_THRESHOLD
    st10["significance_rule"] = "pooled P.adj < 0.05"
    st10["source_file"] = "consolidated input: tableA_step1_marginal_all_regions.tsv"
    st10.to_csv(os.path.join(run,"ST10_Step1_all_DRAFT.tsv"), sep="\t", index=False, float_format="%.17g")
    st10s = [dict(region=r, pooled_test_count=M_BONF, bonferroni_rawP_threshold=RAW_THRESHOLD,
                  significance_rule="pooled P.adj < 0.05", n_step1_significant=int(sigc[r]),
                  n_step1_results_available=tested_by_region[r], n_step1_source_reported_na=sna_by_region[r],
                  min_P=float(A[A.region==r]["P_num"].min()), min_Padj=float(A[A.region==r]["P.adj_num"].min()),
                  source_file="consolidated input: tableA_step1_marginal_all_regions.tsv") for r in EXPECTED_REGIONS]
    pd.DataFrame(st10s).to_csv(os.path.join(run,"ST10_region_summary_DRAFT.tsv"), sep="\t", index=False, float_format="%.17g")

    # -------------------- legacy inventory (NOT downstream config) --------------------
    leg = read_str(LEGACY_MAP)
    legsum = (leg.groupby(["group","rule_id"]).size().reset_index(name="n_labels")
                .rename(columns={"group":"legacy_group"}).sort_values(["legacy_group","n_labels"], ascending=[True,False]))
    legsum["status"] = "LEGACY_NOT_APPROVED"; legsum["note"] = "legacy 16-group resolver; not a new 11-group config"
    legsum.to_csv(os.path.join(run,"legacy_16group_rule_summary_CANDIDATE.tsv"), sep="\t", index=False)
    leg_full = leg.copy(); leg_full["status"] = "LEGACY_NOT_APPROVED"
    leg_full.to_csv(os.path.join(run,"legacy_mapping_rule_inventory.tsv"), sep="\t", index=False)
    SPECIAL = {"AS_n":"_as_","Fibroblast":"fibroblast","VLMC":"_vlmc_","vascular_leptomeningeal":"vascular_leptomeningeal",
               "Mural":"mural","endothelial":"endothelial","pericyte":"pericyte","Vas_*":"_vas_",
               "Microglia":"microglia","Oligo":"oligodendrocyte","OPC":"_opc_"}
    src = open(LEGACY_RULE_SOURCE, encoding="utf-8").read(); tree = ast.parse(src)
    rule_rows = [dict(order=0, line=64, function="normalize",
                      condition_or_pattern=r"lower; sep->_; trailing (_\d+)+ removed",
                      match_type="NORMALIZE", returned_legacy_group="", rule_id="normalize",
                      dataset_specific=False,
                      risk_flag="AS_8->as (trailing _<digits> removed => AS_n collapses to token _as_)",
                      special_tokens="AS_n", note="normalize behaviour")]
    order = 0
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(getattr(t,"id","")=="_RAW_RULES" for t in node.targets):
            for elt in node.value.elts:
                grp = ast.literal_eval(elt.elts[0]); rid = ast.literal_eval(elt.elts[1]); pat = ast.literal_eval(elt.elts[2]); order += 1
                rule_rows.append(dict(order=order, line=elt.lineno, function="_RAW_RULES/resolve",
                    condition_or_pattern=pat, match_type=match_type(pat), returned_legacy_group=grp, rule_id=rid,
                    dataset_specific=dset_specific(pat), risk_flag=risk_flag(pat),
                    special_tokens=";".join(k for k,tok in SPECIAL.items() if tok in pat), note=""))
    rule_rows.append(dict(order=order+1, line=198, function="resolve", condition_or_pattern="norm==''",
        match_type="FALLBACK", returned_legacy_group="Other", rule_id="empty", dataset_specific=False, risk_flag="", special_tokens="fallback", note="empty label"))
    rule_rows.append(dict(order=order+2, line=203, function="resolve", condition_or_pattern="no rule matched",
        match_type="FALLBACK", returned_legacy_group="Other", rule_id="unmatched", dataset_specific=False, risk_flag="", special_tokens="fallback", note="first-match-wins fallthrough"))
    pd.DataFrame(rule_rows).to_csv(os.path.join(run,"legacy_rule_source_inventory.tsv"), sep="\t", index=False)
    pd.DataFrame([dict(canonical_group="", definition="", priority="", dataset_or_resource_scope="",
        exact_rules="", prefix_rules="", regex_rules="", manual_overrides="", status="USER_DEFINITION_REQUIRED")]
    ).to_csv(os.path.join(run,"harmonization_config_v2_TEMPLATE.tsv"), sep="\t", index=False)

    # -------------------- manifest finalization --------------------
    write_manifest(run, inputs, "MAPPING_INDEPENDENT_COMPLETE",
        "mapping-dependent outputs require user 11-group config + APPROVED map",
        extra=dict(step3_audit=dict(raw_rows=tot_rows, unique_pairs=tot_pairs, unique_signals=tot_sigs,
            entered_eq_participants=eq_sig_part, entered_eq_observed=eq_sig_obs, expected_eq_observed_pairs=eq_pairs),
            padj_formula_discrepancy=disc_stats, retained_total=len(retained_signal_set),
            entered_step3_total=sum(len(entered_signals[r]) for r in EXPECTED_REGIONS),
            source_reported_na_audit=source_na_audit))
    print("\n== OUTPUT FILES ==")
    for fn in sorted(os.listdir(run)): print(f"  {fn}  sha256={sha256(os.path.join(run,fn))[:16]}...")
    print("\nMAPPING-DEPENDENT GATE: HALTED (need user 11-group config + APPROVED map).")
    print("RUN DIR:", run)


if __name__ == "__main__":
    main()
