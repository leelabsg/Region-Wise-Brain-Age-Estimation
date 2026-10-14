"""ST21 / ST22 / ST23: the regional BAG polygenic scores in ADNI.

    ST21  Pearson r between PRS_region and the measured ADNI BAG, baseline
          visits only. The BAG is the stored CN-based bias-corrected value
          ({region}_corrected_delta_age of PRS_ADNI_BAG_WIDE); it is read, not
          refitted. Also reports the orientation check below.
    ST22  Logit(label ~ PRS + AGE + SEX + PC1..PC10 + APOE4)     primary
    ST23  Logit(label ~ PRS + AGE + SEX + PC1..PC10)             sensitivity

    label: AD vs CN and MCI vs CN. Benjamini-Hochberg FDR across all tests
    within each table. PRS is raw SCORE1_SUM -- no residualisation, no
    z-scaling -- matching the UK Biobank PheWAS.

Orientation check: a correctly oriented weight set must correlate POSITIVELY
with the measured brain-age gap in an independent cohort. It is reported, not
enforced, because it is a property of the weights rather than of this script;
a negative result means step 01's allele mapping should be re-examined before
anything here is interpreted.

    set -a; . gwas/config/paths.env; set +a
    python gwas/06_prs/06_analyze_adni.py
"""
import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

from _regions import LABEL, REGIONS, require

warnings.filterwarnings("ignore")
PRS_WORK_DIR, ADNI_COVAR, ADNI_APOE, ADNI_PRED_TMPL, ADNI_BAG_WIDE = require(
    "PRS_WORK_DIR", "PRS_ADNI_COVAR", "PRS_ADNI_APOE_RAW", "PRS_ADNI_PRED_TMPL",
    "PRS_ADNI_BAG_WIDE")
# The brain-age model writes the whole-brain predictions under the token "imgs",
# not "global"; every other region uses its own code. This is a filename quirk of
# the model outputs, not a different phenotype.
PRED_TOKEN = dict(zip(REGIONS, REGIONS))
PRED_TOKEN["global"] = os.environ.get("PRS_ADNI_PRED_GLOBAL_TOKEN", "imgs")
OUT_DIR = os.environ.get("ST21_23_OUT_DIR", os.path.join(PRS_WORK_DIR, "tables"))
APOE_VARIANT = os.environ.get("PRS_ADNI_APOE_VARIANT", "19:45411941")
PCS = [f"PC{i}" for i in range(1, 11)]
os.makedirs(OUT_DIR, exist_ok=True)

prs = None
for r in REGIONS:
    s = pd.read_csv(f"{PRS_WORK_DIR}/scores_adni/{r}.sscore", sep="\t").rename(
        columns={"#IID": "IID", "SCORE1_SUM": f"PRS_{r}"})[["IID", f"PRS_{r}"]]
    prs = s if prs is None else prs.merge(s, on="IID", how="outer")

# ---------------- ST21 and the orientation check --------------------------
# Baseline rows of the prediction files link Subject to the BAG table, whose key is
# eid = "sub-" + Subject without "_" + "X" + baseline acquisition date (YYYYMMDD).
# Any unlinked baseline row, a prediction that disagrees with the BAG table, or a
# participant set that differs between BAGs stops the run.
print("=== ST21: corr(PRS, measured ADNI BAG, baseline) ===", flush=True)
wide = pd.read_csv(ADNI_BAG_WIDE)
if wide["eid"].duplicated().any():
    raise SystemExit("ERROR: duplicate eid in PRS_ADNI_BAG_WIDE")
rows, npos, global_r, samples = [], 0, None, {}
for r in REGIONS:
    col = f"{r}_corrected_delta_age"
    if col not in wide.columns or f"{r}_predicted_age" not in wide.columns:
        raise SystemExit(f"ERROR: PRS_ADNI_BAG_WIDE lacks {col} / {r}_predicted_age")
    b = pd.read_csv(ADNI_PRED_TMPL.format(region=PRED_TOKEN[r]))
    b = b[b["is_baseline"] == True].drop_duplicates("Subject").copy()
    n_baseline = len(b)
    b["eid"] = ("sub-" + b["Subject"].str.replace("_", "") + "X"
                + b["Acq Date"].str.replace("-", ""))
    b = b.merge(wide[["eid", col, f"{r}_predicted_age"]], on="eid", how="inner", validate="1:1")
    if len(b) != n_baseline:
        raise SystemExit(f"ERROR: {r}: baseline rows not linked to PRS_ADNI_BAG_WIDE")
    if (b["predicted_brain_age"] - b[f"{r}_predicted_age"]).abs().max() > 1e-4:
        raise SystemExit(f"ERROR: {r}: predictions disagree with PRS_ADNI_BAG_WIDE")
    b = b[["Subject", col]].rename(columns={"Subject": "IID", col: "BAG"})
    m = prs[["IID", f"PRS_{r}"]].merge(b, on="IID", how="inner").dropna()
    samples[r] = frozenset(m["IID"])
    rr, pp = pearsonr(m[f"PRS_{r}"], m["BAG"])
    npos += rr > 0
    if r == "global":
        global_r = rr
    rows.append({"BAG phenotype": LABEL[r], "Pearson r": round(rr, 4), "p value": pp,
                 "N": len(m), "direction": "positive" if rr > 0 else "negative"})
    print(f"  {LABEL[r]:16s} r={rr:+.4f} p={pp:.3e} n={len(m)}  "
          f"{'POS' if rr > 0 else 'NEG'}", flush=True)
if len(set(samples.values())) != 1:
    raise SystemExit("ERROR: ST21 participant sets differ between BAGs")
pd.DataFrame(rows).to_csv(f"{OUT_DIR}/ST21.tsv", sep="\t", index=False)
oriented = npos >= 9 and global_r is not None and global_r >= 0.03
print(f"\nORIENTATION CHECK: {npos}/{len(REGIONS)} positive, global r={global_r:+.4f} "
      f"-> {'consistent' if oriented else 'INCONSISTENT, re-check step 01'}", flush=True)

# ---------------- ST22 / ST23 --------------------------------------------
cov = pd.read_csv(ADNI_COVAR, sep="\t")
for c in ["AGE", "SEX"] + PCS:
    cov[c] = pd.to_numeric(cov[c], errors="coerce")
cov = cov[cov["DX"].isin(["CN", "MCI", "AD"])].copy()
ap = pd.read_csv(ADNI_APOE, sep=r"\s+")
e4 = [c for c in ap.columns if c.startswith(APOE_VARIANT)][0]
ap["APOE4"] = ap[e4].round().clip(0, 2)
df = cov.merge(ap[["IID", "APOE4"]], on="IID", how="left").merge(prs, on="IID", how="left")


def run(cols, tag):
    res = []
    for lab in ["AD", "MCI"]:
        sub = df[df["DX"].isin(["CN", lab])].copy()
        sub["y"] = (sub["DX"] == lab).astype(int)
        for r in REGIONS:
            d = sub[["y", f"PRS_{r}"] + cols].dropna()
            if d["y"].sum() < 30:
                continue
            m = sm.Logit(d["y"], sm.add_constant(d[[f"PRS_{r}"] + cols])).fit(disp=0, maxiter=100)
            if not m.mle_retvals["converged"]:
                continue
            b = m.params[f"PRS_{r}"]
            res.append({"Region": LABEL[r], "label": lab, "n_case": int(d["y"].sum()),
                        "n_control": int((d["y"] == 0).sum()), "coef": b,
                        "odds_ratio": np.exp(b), "p_value": m.pvalues[f"PRS_{r}"]})
    o = pd.DataFrame(res)
    o["fdr"] = multipletests(o["p_value"], method="fdr_bh")[1]
    o.to_csv(f"{OUT_DIR}/{tag}.tsv", sep="\t", index=False)
    print(f"\n=== {tag} ({len(o)} tests): FDR<0.05 -> {(o.fdr < 0.05).sum()} ===", flush=True)
    print(o.sort_values("p_value").to_string(index=False), flush=True)


run(["AGE", "SEX"] + PCS + ["APOE4"], "ST22")
run(["AGE", "SEX"] + PCS, "ST23")
print(f"\nwrote {OUT_DIR}/ST21.tsv, ST22.tsv, ST23.tsv")
