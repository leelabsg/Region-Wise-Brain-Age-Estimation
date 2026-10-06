#!/usr/bin/env python3
"""Build the genetic/phenotypic correlation table from LDSC --rg logs.

- rg comes from the 'Genetic Correlation:' line and p from the standalone 'P:'
  line. The value in parentheses on the rg line is its SE, not a p-value.
- Benjamini-Hochberg within three independent families as first run:
    (1) inter-regional 45 pairs, (2) Telomere x 10, (3) Protein x 10.
  The published table's external-aging family (21 tests) and the ProtAge-204
  rows are produced by 07_protage/11_ldsc_downstream.py.
- Phenotypic columns are Pearson correlations on individual-level
  inverse-normal-transformed brain-age gaps; they require access to the cohort
  phenotype tables and are NA when those are absent.

    set -a; . gwas/config/paths.env; set +a
    python gwas/04_ldsc/build_rg_table.py
"""
import os, re, sys
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests


def env(name, default=None):
    v = os.environ.get(name, default)
    if not v:
        sys.exit(f"ERROR: {name} is not set (see gwas/config/paths.env.example)")
    return v


BASE = env("LDSC_WORK_DIR")
BRAINS_DIR = os.environ.get("RG_BRAINS_DIR", os.path.join(BASE, "results", "gene_corr_brains"))
RES_DIR = os.environ.get("RG_RESULTS_DIR", os.path.join(BASE, "results"))

# display name -> region code
NAME2CODE = {
    "Whole Brain": "global", "Caudate": "caudate", "Cerebellum": "cerebellum",
    "Frontal Lobe": "frontal_lobe", "Insula": "insula", "Occipital Lobe": "occipital_lobe",
    "Parietal Lobe": "parietal_lobe", "Putamen": "putamen", "Temporal Lobe": "temporal_lobe",
    "Thalamus": "thalamus",
}
REGION_NAMES = list(NAME2CODE.keys())

# float incl scientific notation
FLT = r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?"
RE_RG = re.compile(r"^Genetic Correlation:\s*(" + FLT + r")\s*\(", re.M)
# P line: broad token so LDSC underflow prints like 'P: 0.' still parse (float('0.')==0.0)
RE_P  = re.compile(r"^P:\s*([0-9.eE+-]+)\s*$", re.M)
RE_Z  = re.compile(r"^Z-score:\s*(" + FLT + r")\s*$", re.M)

def parse_log(path):
    """Return (rg, p, z) parsing rg line and standalone P: / Z-score: lines separately."""
    with open(path) as f:
        txt = f.read()
    mrg = RE_RG.search(txt)
    mp  = RE_P.search(txt)
    mz  = RE_Z.search(txt)
    if not mrg:
        raise ValueError(f"no 'Genetic Correlation:' in {path}")
    if not mp:
        raise ValueError(f"no 'P:' line in {path}")
    z = float(mz.group(1)) if mz else np.nan
    return float(mrg.group(1)), float(mp.group(1)), z

def find_pair_log(c1, c2):
    """Find inter-regional log for code pair in either ordering."""
    for a, b in ((c1, c2), (c2, c1)):
        cand = os.path.join(BRAINS_DIR, f"gene_corr_{a}_{b}.log")
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(f"no inter-regional log for {c1} x {c2}")

# ---- Build the exact 65-row order (user-specified) ----
rows = []

# family 1: inter-regional 45 pairs (upper triangle in REGION_NAMES order)
for i in range(len(REGION_NAMES)):
    for j in range(i + 1, len(REGION_NAMES)):
        t1, t2 = REGION_NAMES[i], REGION_NAMES[j]
        log = find_pair_log(NAME2CODE[t1], NAME2CODE[t2])
        rg, p, z = parse_log(log)
        rows.append({"Trait1": t1, "Trait2": t2, "family": "inter", "Genetic_r": rg,
                     "Genetic_Z": z, "Genetic_p": p, "log": os.path.basename(log)})

# family 2: Telomere Length PC1 x 10
for t2 in REGION_NAMES:
    log = os.path.join(RES_DIR, f"gene_corr_telomere_{NAME2CODE[t2]}.log")
    rg, p, z = parse_log(log)
    rows.append({"Trait1": "Telomere Length PC1", "Trait2": t2, "family": "telomere", "Genetic_r": rg,
                 "Genetic_Z": z, "Genetic_p": p, "log": os.path.basename(log)})

# family 3: Protein Age x 10
for t2 in REGION_NAMES:
    log = os.path.join(RES_DIR, f"gene_corr_protein_{NAME2CODE[t2]}.log")
    rg, p, z = parse_log(log)
    rows.append({"Trait1": "Protein Age", "Trait2": t2, "family": "protein", "Genetic_r": rg,
                 "Genetic_Z": z, "Genetic_p": p, "log": os.path.basename(log)})

df = pd.DataFrame(rows)

# ---- Underflowed p-values: LDSC printed 'P: 0.' when p < smallest double (~5e-324) ----
DENORM_MIN = 5e-324          # np.nextafter(0, 1); smallest positive double
underflow = df["Genetic_p"] == 0.0
# p used for BH: replace exact 0 with the smallest double so ranking is unchanged but FDR != 0
df["p_for_bh"] = df["Genetic_p"].where(~underflow, DENORM_MIN)
# display string for the table
df["Genetic_p_str"] = df["Genetic_p"].map(lambda v: "%.6g" % v)
df.loc[underflow, "Genetic_p_str"] = "<2.2e-308"

# ---- BH-FDR within each family independently (using p_for_bh) ----
df["Genetic_fdr"] = np.nan
for fam in ["inter", "telomere", "protein"]:
    m = df["family"] == fam
    q = multipletests(df.loc[m, "p_for_bh"].values, alpha=0.05, method="fdr_bh")[1]
    df.loc[m, "Genetic_fdr"] = q

# ---- Phenotypic correlation: Pearson on individual-level regional BAG ----
# Data: user-supplied GWAS input (British baseline).
# Column: {region}_corrected_delta_age_int (inverse-normal-transformed BAG; matches GWAS phenoCol).
from scipy.stats import pearsonr
PHENO_CSV = env("UKB_PHENO_FILE")
# GWAS (SAIGE) ran on the inverse-normal-transformed phenotype (phenoCol = <region>,
# confirmed = {code}_corrected_delta_age_int). Use _int so rp matches rg's scale.
COL_INT = "{code}_corrected_delta_age_int"
COL_RAW = "{code}_corrected_delta_age"
pheno = pd.read_csv(PHENO_CSV)

# Protein Age individual-level data (proteomic-age gap), merged to BAG on IID.
# Same phenotype construction as regions: corrected_delta_age(_int). rp on the overlap sample.
PROT_CSV = env("PROTEIN_AGE_PHENO_FILE")
prot = pd.read_csv(PROT_CSV, usecols=["IID", "corrected_delta_age", "corrected_delta_age_int"]) \
         .rename(columns={"corrected_delta_age": "prot_raw", "corrected_delta_age_int": "prot_int"})
bagprot = pheno.merge(prot, on="IID", how="inner")   # overlap with both measures

df["Phenotypic_r"] = np.nan     # <- primary (int); goes to the table
df["Phenotypic_p"] = np.nan
df["rp_raw"] = np.nan
df["Phenotypic_N"] = np.nan
pheno_N_inter = set()
pheno_N_prot = set()
for idx, r in df.iterrows():
    fam = r["family"]
    if fam == "inter":
        a, b = NAME2CODE[r["Trait1"]], NAME2CODE[r["Trait2"]]
        i1, i2 = pheno[COL_INT.format(code=a)], pheno[COL_INT.format(code=b)]
        r1, r2 = pheno[COL_RAW.format(code=a)], pheno[COL_RAW.format(code=b)]
        mask = i1.notna() & i2.notna() & r1.notna() & r2.notna()
        n = int(mask.sum())
        rr, pp = pearsonr(i1[mask], i2[mask])
        df.at[idx, "Phenotypic_r"] = rr
        df.at[idx, "Phenotypic_p"] = pp
        df.at[idx, "rp_raw"] = pearsonr(r1[mask], r2[mask])[0]
        df.at[idx, "Phenotypic_N"] = n
        pheno_N_inter.add(n)
    elif fam == "protein":
        # Protein Age (Trait1) x region (Trait2) over the overlap sample
        code = NAME2CODE[r["Trait2"]]
        reg_i, reg_r = bagprot[COL_INT.format(code=code)], bagprot[COL_RAW.format(code=code)]
        mask = reg_i.notna() & reg_r.notna() & bagprot["prot_int"].notna() & bagprot["prot_raw"].notna()
        n = int(mask.sum())
        rr, pp = pearsonr(bagprot["prot_int"][mask], reg_i[mask])
        df.at[idx, "Phenotypic_r"] = rr
        df.at[idx, "Phenotypic_p"] = pp
        df.at[idx, "rp_raw"] = pearsonr(bagprot["prot_raw"][mask], reg_r[mask])[0]
        df.at[idx, "Phenotypic_N"] = n
        pheno_N_prot.add(n)
    # telomere: no individual-level telomere data supplied -> NA

# ---- Phenotypic FDR: Benjamini-Hochberg within each family (same structure as Genetic_fdr) ----
# inter-regional (45) and Protein Age (10) separately; Telomere rows stay NA.
df["Phenotypic_fdr"] = np.nan
for fam in ["inter", "protein"]:
    m = (df["family"] == fam) & df["Phenotypic_p"].notna()
    p = df.loc[m, "Phenotypic_p"].to_numpy()
    p_bh = np.where(p == 0.0, DENORM_MIN, p)   # substitute underflow 0 for ranking (inter)
    df.loc[m, "Phenotypic_fdr"] = multipletests(p_bh, alpha=0.05, method="fdr_bh")[1]

# ---- Write TSV: ...|Genetic_fdr|Phenotypic_r|Phenotypic_p|Phenotypic_fdr ----
df["Genetic_r"] = df["Genetic_r"].map(lambda v: "%.4g" % v)
df["Genetic_Z"] = df["Genetic_Z"].map(lambda v: "%.6g" % v)
df["Genetic_fdr"] = df["Genetic_fdr"].map(lambda v: "%.4g" % v)
df["Phenotypic_r_str"] = df["Phenotypic_r"].map(lambda v: "NA" if pd.isna(v) else "%.4g" % v)
# Phenotypic_p / Phenotypic_fdr: NA for telomere; '<2.2e-308' where value underflowed below double precision
def _pstr(v):
    if pd.isna(v): return "NA"
    return "<2.2e-308" if v < 2.2e-308 else "%.4g" % v
df["Phenotypic_p_str"]   = df["Phenotypic_p"].map(_pstr)
df["Phenotypic_fdr_str"] = df["Phenotypic_fdr"].map(_pstr)

# select the string versions explicitly (avoid rename collision with numeric Genetic_p)
out = df[["Trait1", "Trait2", "Genetic_r", "Genetic_Z",
          "Genetic_p_str", "Genetic_fdr", "Phenotypic_r_str", "Phenotypic_p_str", "Phenotypic_fdr_str"]].copy()
out.columns = ["Trait1", "Trait2", "Genetic_r", "Genetic_Z",
               "Genetic_p", "Genetic_fdr", "Phenotypic_r", "Phenotypic_p", "Phenotypic_fdr"]
OUT_PATH = env("RG_TABLE_OUT", os.path.join(BASE, "outputs", "tables", "rg_table.tsv"))
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
# dynamic N strings from complete-case counts (no hardcoding)
def _nstr(s):
    return str(next(iter(s))) if len(s) == 1 else f"{min(s)}-{max(s)}"
n_inter, n_prot = _nstr(pheno_N_inter), _nstr(pheno_N_prot)
with open(OUT_PATH, "w") as fh:
    # footnote FIRST (comment lines above the header, matching ST5-8 format)
    fh.write("# ST13. Genetic and phenotypic correlations among regional brain-age gaps and with biomarkers.\n")
    fh.write("# Genetic_r/Genetic_Z/Genetic_p: LDSC cross-trait rg ('Genetic Correlation:' line; P from the separate 'P:' line, not SE).\n")
    fh.write("# Genetic_p '<2.2e-308' = LDSC printed 'P: 0.' (underflow); FDR used p=5e-324 for that row (rank unchanged).\n")
    fh.write("# Genetic_fdr: Benjamini-Hochberg within each family separately (inter-regional 45; Telomere x10; Protein x10).\n")
    fh.write("# Phenotypic_r/Phenotypic_p: Pearson r and two-sided p on inverse-normal-transformed phenotypes ({region}_corrected_delta_age_int).\n")
    fh.write("# Phenotypic_fdr: Benjamini-Hochberg within each family separately (inter-regional 45; Protein x10), same structure as Genetic_fdr. '<2.2e-308' = below double precision.\n")
    fh.write(f"#   inter-regional: both regional BAGs, British baseline N={n_inter}; Phenotypic_p/fdr '<2.2e-308' = below double precision at this N (all 45 pairs).\n")
    fh.write(f"#   Protein Age x region: proteomic-age gap x regional BAG on the overlap sample N={n_prot}. Telomere Length PC1 rows = NA (no individual-level telomere data supplied).\n")
    fh.write(f"# Phenotypic source (BAG): {PHENO_CSV}\n")
    fh.write(f"# Phenotypic source (Protein Age): {PROT_CSV}\n")
    out.to_csv(fh, sep="\t", index=False)
print(f"WROTE {OUT_PATH}  ({len(out)} rows)")

# ---- Summary ----
inter = df[df.family == "inter"]
tel = df[df.family == "telomere"]
pro = df[df.family == "protein"]
inter_r = inter["Genetic_r"].astype(float)
inter_fdr = inter["Genetic_fdr"].astype(float)
print("\n================ SUMMARY ================")
print(f"Total rows: {len(df)}  (inter={len(inter)}, telomere={len(tel)}, protein={len(pro)})")

# 1) p==0 underflow report
uf = df[underflow]
print(f"\n[p underflow] {len(uf)} row(s) had LDSC 'P: 0.' (true p < ~5e-324):")
for _, r in uf.iterrows():
    print(f"    {r.Trait1} x {r.Trait2:16s} Z={r.Genetic_Z}  -> Genetic_p reported as '<2.2e-308', BH used {DENORM_MIN:g}")

print(f"\n[inter-regional] rg range: {inter_r.min():.4f} .. {inter_r.max():.4f}")
print(f"[inter-regional] FDR<0.05: {(inter_fdr<0.05).sum()} / 45")
tel_sig = tel[tel['Genetic_fdr'].astype(float) < 0.05]["Trait2"].tolist()
pro_sig = pro[pro['Genetic_fdr'].astype(float) < 0.05]["Trait2"].tolist()
print(f"\n[Telomere] FDR<0.05: {len(tel_sig)} / 10 -> {tel_sig}")
print(f"[Protein]  FDR<0.05: {len(pro_sig)} / 10 -> {pro_sig}")

ph_int = inter["Phenotypic_r"].astype(float)   # int (table value)
ph_raw = inter["rp_raw"].astype(float)
print(f"\n[Phenotypic] method: Pearson on INT ('{{region}}_corrected_delta_age_int')")
print(f"[Phenotypic][inter] N (all 45): {sorted(pheno_N_inter)}")
print(f"[Phenotypic][inter] rp_int range: {ph_int.min():.4f} .. {ph_int.max():.4f}   rp_raw range: {ph_raw.min():.4f} .. {ph_raw.max():.4f}")

# Protein Age x region phenotypic vs genetic (overlap sample)
pro_ph = pro["Phenotypic_r"].astype(float)
print(f"\n[Phenotypic][Protein Age x region] N (overlap BAG n Protein): {sorted(pheno_N_prot)}")
print(f"[Phenotypic][Protein Age x region] rp range: {pro_ph.min():.4f} .. {pro_ph.max():.4f}")
prot_tab = pro[["Trait2","Genetic_r","Phenotypic_r","Phenotypic_p","Genetic_fdr"]].copy()
prot_tab["rp_int"] = pro_ph.round(4).values
prot_tab["rp_p"]   = pro["Phenotypic_p"].astype(float).map(lambda v: f"{v:.2e}").values
print(pro[["Trait2","Genetic_r","Genetic_fdr"]].assign(
        rp_int=pro_ph.round(4).values,
        rp_p=pro["Phenotypic_p"].astype(float).map(lambda v: f"{v:.2e}").values,
        rp_fdr=pro["Phenotypic_fdr"].astype(float).map(lambda v: f"{v:.2e}").values,
        rp_fdr_sig=(pro["Phenotypic_fdr"].astype(float) < 0.05).values).to_string(index=False))
n_pro_ph_sig = int((pro["Phenotypic_fdr"].astype(float) < 0.05).sum())
print(f"[Phenotypic][Protein Age x region] FDR<0.05: {n_pro_ph_sig} / 10")

# raw vs int side-by-side (numbers only; no interpretation)
cmp = inter[["Trait1", "Trait2"]].copy()
cmp["rp_raw"] = ph_raw.values
cmp["rp_int"] = ph_int.values
cmp["diff"] = ph_int.values - ph_raw.values
maxabs = cmp["diff"].abs().max()
rm = cmp.loc[cmp["diff"].abs().idxmax()]
n2 = int((cmp["rp_raw"].round(2) == cmp["rp_int"].round(2)).sum())
print(f"\n[raw vs int] max |diff| = {maxabs:.5f}  at  {rm.Trait1} x {rm.Trait2}")
print(f"[raw vs int] pairs agreeing to 2 decimals: {n2} / 45")
print(cmp.to_string(index=False, formatters={
    "rp_raw": lambda v: f"{v:.4f}", "rp_int": lambda v: f"{v:.4f}", "diff": lambda v: f"{v:+.5f}"}))

pd.set_option("display.width", 240); pd.set_option("display.max_rows", 200)
disp = inter[["Trait1", "Trait2", "Genetic_r"]].copy()
disp["Phenotypic_r"] = ph_int.round(4).values
print("\n---- inter-regional (all 45): Trait1 | Trait2 | Genetic_r | Phenotypic_r ----")
print(disp.to_string(index=False))
print("\n---- inter-regional full (rg | Z | p | fdr | rp) ----")
print(inter[["Trait1","Trait2","Genetic_r","Genetic_Z","Genetic_p_str","Genetic_fdr","Phenotypic_r"]].to_string(index=False))
print("\n---- Telomere (10) ----")
print(tel[["Trait2","Genetic_r","Genetic_Z","Genetic_p_str","Genetic_fdr"]].to_string(index=False))
print("\n---- Protein (10) ----")
print(pro[["Trait2","Genetic_r","Genetic_Z","Genetic_p_str","Genetic_fdr"]].to_string(index=False))
