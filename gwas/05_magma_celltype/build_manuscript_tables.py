#!/usr/bin/env python3
"""
Manuscript integration tables (A/B) + cross-region summaries + cerebellum
asymmetry diagnostics + per-region lambda_GC, from the MAGMA cell-type output.

Input-file routing is configuration, not a literal. The manuscript-table build
read the finalized outputs for all ten BAG phenotypes from the consolidated local
STEP23 directory. This input-file routing does not change the analysis
provenance: nine phenotypes were analysed with the FUMA web implementation and
cerebellar BAG was analysed locally, as described in Supplementary Note 7.
Naming regions in MAGMA_FUMA_REGIONS reads those from a FUMA download directory
instead; it is empty by default, matching the build above.

Ground rules (per user):
- NO transformation of underlying MAGMA numbers. Schemas verified byte-identical.
- conditional p-values (step2/3) get NO extra multiple-testing correction.
- No harmonization / grouping / exclusion in the main tables.
- Cerebellum carries far more significant cell types than any other region. The inflation
  diagnostic is NOT outstanding: cerebellum lambda_GC is 1.1233, inside the range the other
  nine regions span. The caution is about interpreting a large count, not an unresolved check.

INDEPENDENT-SET DEFINITION (verified 2026-07-16):
  magma_celltype_step3.txt has NO final pass/drop flag (16 cols = pairwise stats +
  descriptive PS only; see FUMA README). The only retention flag is
  step1_2_summary.step3 = "retained FOR step3" = survival of step2 within-dataset
  forward selection. Therefore the independent-significant set == the step2-retained
  set. step3 PS/CDM.P are attached as cross-dataset annotation only (no thresholding).
  survived_at label "entered_step3" means: step2-retained AND region ran cross-dataset
  step3 (>=2 significant datasets) -- NOT a step3 pass verdict.
"""
import os, re
import numpy as np
import pandas as pd
from scipy.stats import chi2

def _env(name):
    import os as _os, sys as _sys
    v = _os.environ.get(name)
    if not v:
        _sys.exit(f"ERROR: {name} is not set (see gwas/config/paths.env.example)")
    return v


ROOT = _env("MAGMA_ROOT")
GWAS = _env("GWAS_UKB_DIR")
OUT  = os.environ.get("MAGMA_TABLES_DIR", os.path.join(ROOT, "outputs/manuscript_tables"))
os.makedirs(OUT, exist_ok=True)

# Where this build READS each region's finalized cell-type output from.
STEP23 = os.environ.get("MAGMA_STEP23_DIR", os.path.join(ROOT, "outputs/step23"))
# Regions whose files are read from a FUMA download directory instead. Empty in the
# build that produced the manuscript tables.
FUMA_DIR = os.environ.get("MAGMA_FUMA_CELLTYPE_DIR",
                          os.path.join(ROOT, "data/fuma_celltype_results"))
FUMA_REGIONS = {r for r in os.environ.get("MAGMA_FUMA_REGIONS", "").split(",") if r}
# Merged per-region summary statistics, relative to GWAS_UKB_DIR, for the lambda_GC panel.
# Default is the naming the published build used.
SUMSTATS_TMPL = os.environ.get("MAGMA_SUMSTATS_TMPL",
                               "{region}/{region}_discovery_41067_sumstats.txt")

# (M_BONF is not defined here: this script never used it. Significance comes from the
#  upstream P.adj compared against SIG_ADJ, so no Bonferroni denominator is applied.)
SIG_ADJ = 0.05
MAF_MIN, INFO_MIN = 0.01, 0.30   # lambda_GC QC filters

REGIONS = ["global","caudate","frontal_lobe","insula","occipital_lobe",
           "parietal_lobe","putamen","temporal_lobe","thalamus","cerebellum"]

# `source` records WHICH DIRECTORY THIS BUILD READ THE FILE FROM -- LOCAL for the
# consolidated STEP23 directory, FUMA for a FUMA download directory. It is not the
# analysis implementation: the FUMA web implementation analysed nine BAG phenotypes and
# cerebellar BAG was analysed locally (Supplementary Note 7), and that provenance is
# recorded there, not in this column. In the build that produced the manuscript tables
# every region was read from STEP23, so the column reads LOCAL throughout. The column is
# kept because removing it would change tableA's shape for every downstream consumer.
region_dir = lambda r: os.path.join(FUMA_DIR, r) if r in FUMA_REGIONS else os.path.join(STEP23, r)
src_of     = lambda r: "FUMA" if r in FUMA_REGIONS else "LOCAL"
readtsv    = lambda p: pd.read_csv(p, sep="\t", dtype=str)
def tonum(df, cols):
    for c in cols:
        if c in df.columns: df[c]=pd.to_numeric(df[c], errors="coerce")
    return df

# ------------------------------------------------------------------ TABLE A
print("== Table A (step1 marginal, all regions) ==")
a_frames=[]
for r in REGIONS:
    df=tonum(readtsv(os.path.join(region_dir(r),"magma_celltype_step1.txt")),
             ["NGENES","BETA","BETA_STD","SE","P","P.adj.pds","P.adj"])
    df.insert(0,"region",r); df.insert(1,"source",src_of(r))
    df["significant"]=(df["P.adj"]<SIG_ADJ).astype(int)
    a_frames.append(df)
tableA=pd.concat(a_frames, ignore_index=True)
tableA.to_csv(os.path.join(OUT,"tableA_step1_marginal_all_regions.tsv"), sep="\t", index=False)
print(f"  rows={len(tableA)}  sig={int(tableA.significant.sum())}")

# ------------------------------------------------------------------ supporting long tables
def concat_optional(fname, numcols):
    frames=[]
    for r in REGIONS:
        p=os.path.join(region_dir(r),fname)
        if os.path.exists(p):
            df=tonum(readtsv(p),numcols); df.insert(0,"region",r); df.insert(1,"source",src_of(r)); frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

step2_long=concat_optional("magma_celltype_step2.txt",["MODEL","NGENES","BETA","BETA_STD","SE","P","Marginal.P","PS"])
summary_long=concat_optional("step1_2_summary.txt",["NGENES","BETA","BETA_STD","SE","P","P.adj.pds","P.adj","step3"])
step3_long=concat_optional("magma_celltype_step3.txt",["MODEL","NGENES","BETA","BETA_STD","SE","P","CDM.BETA","CDM.BETA_STD","CDM.SE","CDM.P","Marginal.P","PS","PS.avg"])
step2_long.to_csv(os.path.join(OUT,"supp_step2_within_dataset_all.tsv"), sep="\t", index=False)
summary_long.to_csv(os.path.join(OUT,"supp_step1_2_summary_all.tsv"), sep="\t", index=False)
step3_long.to_csv(os.path.join(OUT,"supp_step3_cross_dataset_all.tsv"), sep="\t", index=False)
print(f"  supp: step2={len(step2_long)} summary={len(summary_long)} step3={len(step3_long)}")

# ------------------------------------------------------------------ TABLE B
print("== Table B (independent-significant set = step2-retained set) ==")
region_has_step3   = {r: os.path.exists(os.path.join(region_dir(r),"magma_celltype_step3.txt")) for r in REGIONS}
region_has_summary = {r: os.path.exists(os.path.join(region_dir(r),"step1_2_summary.txt")) for r in REGIONS}
s3=step3_long.copy()
def step3_annot(region,dataset,cell):
    if s3.empty: return None
    m=s3[(s3.region==region)&(s3.Dataset==dataset)&(s3.Cell_type==cell)]
    if len(m)==0: return None
    return dict(n_step3_pairs=len(m), step3_P_cond_min=np.nanmin(m["P"]), step3_P_cond_max=np.nanmax(m["P"]),
                step3_PS_min=np.nanmin(m["PS"]), step3_PS_mean=np.nanmean(m["PS"]), step3_CDM_P_min=np.nanmin(m["CDM.P"]))

SUMMARY_REQUIRED=["Dataset","Cell_type","NGENES","BETA","SE","P","P.adj","cond_state","cond_cell_type","step3"]
try:
    b_rows=[]
    for r in REGIONS:
        sig=tableA[(tableA.region==r)&(tableA.significant==1)]
        if len(sig)==0: continue
        if region_has_summary[r]:
            sm=readtsv(os.path.join(region_dir(r),"step1_2_summary.txt"))
            missing=[c for c in SUMMARY_REQUIRED if c not in sm.columns]
            if missing: raise KeyError(f"region={r} step1_2_summary missing {missing}; actual={list(sm.columns)}")
            sm=tonum(sm,["NGENES","BETA","BETA_STD","SE","P","P.adj.pds","P.adj","step3"])
            for _,row in sm.iterrows():
                retained=int(row["step3"])==1
                survived = "dropped_at_step2" if not retained else ("entered_step3" if region_has_step3[r] else "step2_retained")
                rec=dict(region=r, source=src_of(r), Dataset=row["Dataset"], Cell_type=row["Cell_type"],
                         NGENES=row["NGENES"], BETA=row["BETA"], SE=row["SE"], P=row["P"], P_adj=row["P.adj"],
                         step2_cond_state=row["cond_state"], step2_cond_celltype=row["cond_cell_type"],
                         step2_retained=int(retained), region_ran_step3=int(region_has_step3[r]), survived_at=survived)
                if survived=="entered_step3":
                    a=step3_annot(r,row["Dataset"],row["Cell_type"])
                    if a: rec.update(a)
                b_rows.append(rec)
        else:
            for _,row in sig.iterrows():
                b_rows.append(dict(region=r, source=src_of(r), Dataset=row["Dataset"], Cell_type=row["Cell_type"],
                    NGENES=row["NGENES"], BETA=row["BETA"], SE=row["SE"], P=row["P"], P_adj=row["P.adj"],
                    step2_cond_state="NA", step2_cond_celltype="NA",
                    step2_retained=1, region_ran_step3=0, survived_at="step1_single"))
    tableB=pd.DataFrame(b_rows)
    # The cerebellum exception was a "this region differs from the other nine" note, and the
    # diagnostic it pointed at (region_lambda_GC.outlier_vs_other9, cerebellum_family_decomposition)
    # is built on exactly that contrast. With all ten regions run the same way it has no referent.
    # Column kept so tableB's shape is unchanged.
    tableB["inflation_note"]=""
    col_order=["region","source","Dataset","Cell_type","NGENES","BETA","SE","P","P_adj",
               "step2_cond_state","step2_cond_celltype","step2_retained","region_ran_step3","survived_at",
               "n_step3_pairs","step3_P_cond_min","step3_P_cond_max","step3_PS_min","step3_PS_mean","step3_CDM_P_min","inflation_note"]
    for c in col_order:
        if c not in tableB.columns: tableB[c]=np.nan
    tableB=tableB[col_order]
    tableB.to_csv(os.path.join(OUT,"tableB_independent_significant_set.tsv"), sep="\t", index=False)
    print("  survived_at breakdown by region:")
    print(tableB.groupby(["region","survived_at"]).size().to_string())
except KeyError as e:
    print("  !! TABLE B KeyError ->", e); raise

# ------------------------------------------------------------------ (a) recurrence (two versions)
print("== Cross-region (a) recurrence ==")
indep=tableB[tableB.survived_at!="dropped_at_step2"]      # independent = step2-retained (+ step1_single)
def recurrence(df, tag):
    r=(df.groupby("Cell_type")
       .agg(n_regions=("region","nunique"), regions=("region",lambda x:";".join(sorted(set(x)))), n_rows=("region","size"))
       .reset_index().sort_values(["n_regions","n_rows"], ascending=False))
    r.to_csv(os.path.join(OUT,f"summaryC_a_recurrence_{tag}.tsv"), sep="\t", index=False)
    return r
rec_all=recurrence(indep, "incl_cerebellum")
rec_9  =recurrence(indep[indep.region!="cerebellum"], "excl_cerebellum")
print(f"  (a1 incl cerebellum) distinct labels={len(rec_all)}; recurring >=2 regions:")
print(rec_all[rec_all.n_regions>=2].to_string(index=False))
print(f"  (a2 EXCL cerebellum, FUMA-verified 9 region) distinct labels={len(rec_9)}; recurring >=2 regions:")
print(rec_9[rec_9.n_regions>=2].to_string(index=False))

# ------------------------------------------------------------------ (b) pooled BH-FDR
print("== Cross-region (b) pooled BH-FDR (all step1 P) ==")
from statsmodels.stats.multitest import multipletests
allp=tableA[["region","source","Dataset","Cell_type","BETA","SE","P","P.adj"]].dropna(subset=["P"]).copy()
rej,q,_,_=multipletests(allp["P"].values, alpha=0.05, method="fdr_bh")
allp["BH_q_pooled_allregions"]=q; allp["BH_sig_q05"]=rej.astype(int)
allp=allp.sort_values("P")
allp.to_csv(os.path.join(OUT,"summaryC_b_pooled_BH_FDR_all_step1.tsv"), sep="\t", index=False)
print(f"  pooled tests={len(allp)}  BH q<0.05 survivors={int(allp.BH_sig_q05.sum())}")
print(allp[allp.BH_sig_q05==1].groupby("region").size().to_string())

# ------------------------------------------------------------------ cerebellum family diagnostics
print("== Cerebellum family / dataset-dispersion diagnostics ==")
def norm(s): return re.sub(r"[._\- ]","",s.lower())
def classify(label):
    l=label.lower(); n=norm(label)
    if n.startswith("opc"):                                   return ("OPC","prefix:opc")
    if "precursor" in l and "olig" in l:                      return ("OPC","kw:oligo+precursor")
    if "olig" in l:                                           return ("oligodendrocyte","substr:olig")
    if n.startswith("nonneu"):                                return ("non_neuronal","prefix:nonneu")
    if "astro" in l:                                          return ("astrocyte","substr:astro")
    # Seeker2023 WhiteMatter level3 uses AS_<n> for astrocyte subclusters
    # (verified from that dataset's own level1='Astrocyte' / level2='astrocyte' hierarchy).
    if re.fullmatch(r"as\d+", n):                             return ("astrocyte","seeker2023:AS_n")
    if "granul" in l or n.startswith("grc"):                  return ("granule","substr:granul/grc")
    if "purkinje" in l:                                       return ("purkinje","substr:purkinje")
    if "microglia" in l:                                      return ("microglia","substr:microglia")
    if any(k in l for k in ("endothel","pericyte","fibroblast","mural","vlmc","vsmc","vascular")):
        return ("vascular","substr:vascular-token")
    return ("other","unmatched")

cer_sig=tableA[(tableA.region=="cerebellum")&(tableA.significant==1)].copy()
fm=cer_sig["Cell_type"].map(classify)
cer_sig["family"]=[f for f,_ in fm]; cer_sig["matched_rule"]=[m for _,m in fm]
dump=(cer_sig.groupby(["Cell_type","family","matched_rule"]).size().reset_index(name="n_sig")
      .sort_values(["family","n_sig"], ascending=[True,False]))
dump.to_csv(os.path.join(OUT,"cerebellum_label_family_map.tsv"), sep="\t", index=False)
print(f"  cerebellum sig cells={len(cer_sig)} across {cer_sig.Dataset.nunique()} datasets; distinct labels={cer_sig.Cell_type.nunique()}")
print("  --- full label->family dump ---")
print(dump.to_string(index=False))
print("  --- 'other' bucket (eyeball) ---")
print(dump[dump.family=="other"].to_string(index=False))
fam=(cer_sig.groupby("family").agg(n_sig_cells=("Cell_type","size"),
     n_distinct_datasets=("Dataset","nunique"), n_distinct_labels=("Cell_type","nunique")).reset_index()
     .sort_values("n_sig_cells", ascending=False))
fam["pct"]=(100*fam.n_sig_cells/len(cer_sig)).round(1)
fam.to_csv(os.path.join(OUT,"cerebellum_family_decomposition.tsv"), sep="\t", index=False)
print("  --- family decomposition ---"); print(fam.to_string(index=False))
per_ds=cer_sig.groupby("Dataset").size()
print(f"  per-dataset sig cells: min={per_ds.min()} median={int(per_ds.median())} max={per_ds.max()} | 1sig={int((per_ds==1).sum())}ds >1sig={int((per_ds>1).sum())}ds")
ol=cer_sig[cer_sig.family.isin(["oligodendrocyte","OPC"])]
print(f"  oligo lineage: {len(ol)}/{len(cer_sig)} cells in {ol.Dataset.nunique()} distinct datasets (same lineage replicated, not {len(cer_sig)} independent signals)")

# ------------------------------------------------------------------ lambda_GC per region
print("== lambda_GC per region (MAF>=%.2f, INFO>=%.2f) ==" % (MAF_MIN, INFO_MIN))
DENOM=chi2.ppf(0.5,1)   # 0.4549364
lam_rows=[]
for r in REGIONS:
    f=os.path.join(GWAS, SUMSTATS_TMPL.format(region=r))
    d=pd.read_csv(f, sep="\t", usecols=["AF_Allele2","imputationInfo","p.value","N"],
                  dtype={"AF_Allele2":"float64","imputationInfo":"float64","p.value":"float64","N":"float64"})
    n_raw=len(d)
    maf=np.minimum(d["AF_Allele2"], 1-d["AF_Allele2"])
    keep=(maf>=MAF_MIN)&(d["imputationInfo"]>=INFO_MIN)&d["p.value"].notna()&(d["p.value"]>0)
    dk=d[keep]
    med_p=float(np.median(dk["p.value"].values))
    lam=float(chi2.isf(med_p,1)/DENOM)
    lam_rows.append(dict(region=r, source=src_of(r), N_SNP_raw=n_raw, N_SNP_QC=int(keep.sum()),
                         N_samples=int(np.nanmax(d["N"].values)), median_p=round(med_p,5), lambda_GC=round(lam,4)))
lam=pd.DataFrame(lam_rows)
lam_ex_cer=lam[lam.region!="cerebellum"]["lambda_GC"]
lo,hi=float(lam_ex_cer.min()), float(lam_ex_cer.max())
cer_lam=float(lam[lam.region=="cerebellum"]["lambda_GC"])
lam["outlier_vs_other9"]=np.where(lam.region=="cerebellum",
    f"cerebellum={cer_lam} vs other9 range [{lo:.3f},{hi:.3f}]", "")
lam.to_csv(os.path.join(OUT,"region_lambda_GC.tsv"), sep="\t", index=False)
print(lam.to_string(index=False))
print(f"  cerebellum lambda_GC={cer_lam} ; other 9 range [{lo:.3f},{hi:.3f}] ; "
      f"cerebellum {'OUTSIDE(higher)' if cer_lam>hi else ('OUTSIDE(lower)' if cer_lam<lo else 'WITHIN')} that range")

# ------------------------------------------------------------------ Methods note
with open(os.path.join(OUT,"METHODS_NOTES.md"),"w") as fh:
    fh.write(
"# Manuscript-table method notes\n\n"
"## Independent-significant set definition\n"
"`magma_celltype_step3.txt` contains NO final pass/drop flag (16 columns = pairwise\n"
"conditional statistics + descriptive PS only; FUMA README confirmed). The only\n"
"retention flag is `step1_2_summary.step3` = 'retained FOR step3' = survival of the\n"
"step2 within-dataset forward selection. **Independent-significant set == step2-retained set.**\n"
"step3 columns (CD-conditional P, PS, PS.avg) are attached as cross-dataset annotation\n"
"only; no threshold is applied to them (step2/3 are themselves the selection).\n\n"
"`survived_at` values: `dropped_at_step2`, `entered_step3` (step2-retained AND region ran\n"
"cross-dataset step3, i.e. >=2 significant datasets -- NOT a step3 pass verdict),\n"
"`step2_retained` (retained but region did not run step3), `step1_single`.\n\n"
"## Cerebellum\n"
"Cerebellum step1/2/3 are from a LOCAL MAGMA run (byte-identical schema to FUMA; same\n"
"12,531 pooled-Bonferroni denominator). Significant count NOT confirmed as a headline\n"
"number pending inflation diagnostics (region_lambda_GC.tsv, cerebellum_family_decomposition.tsv).\n")

print("\nDONE ->", OUT)
