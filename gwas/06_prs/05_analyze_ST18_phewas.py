"""ST18: PheWAS of the regional BAG polygenic scores over UK Biobank phenotypes.

Three analysis samples are reported. Each is a separate run of this script with
its own ST18_SAMPLE and its own output directory; they are never written to the
same place and are never pooled.

  ST18_SAMPLE                  analysis sample                              n
  ukb_nonimaging_wb            UK Biobank non-imaging, White British    389,045   primary
  ukb_nonimaging               UK Biobank non-imaging, unrestricted     442,295   sensitivity
  ukb_imaging_heldout_10267    UK Biobank imaging held-out               10,267

Models
  A  non-sex-specific phecode
         logit(case) ~ PRS_region + Age + Sex + PC1..PC10
  B  sex-specific phecode, restricted to the relevant sex, cases and controls
     alike; Sex is dropped from the design because it is constant after the
     restriction
         logit(case) ~ PRS_region + Age + PC1..PC10
  C  phecode 751.1, fitted unrestricted with the model-A design: its constituent
     subcodes carry conflicting sex designations, so no single restriction
     applies

Test units. A phecode with two decimals is folded into its one-decimal unit
(XX174.11 -> XX174.1); a phecode with no decimal is its own unit (XX185). One-
decimal and integer codes are therefore SEPARATE units -- 174.1 is not merged
into 174. Within a unit the child columns are combined with a skipna maximum, so
a participant who is a case at any child is a case for the unit, [0, NaN] is a
control and only [NaN, ...NaN] stays NaN.

A NaN in the phecode table is the phecode exclusion range for a related disease
-- neither case nor control -- and is dropped from that unit's test after the
aggregation above, rather than being recoded as a control. After the sex
restriction and complete-case selection on the covariates and the PRS, a test is
fitted when it has at least 30 cases and both outcome classes. Fits that raise,
that do not converge at maxiter=50, or that return a non-finite coefficient, SE
or p-value are excluded and recorded in ST18_fit_status.tsv with their reason.

Benjamini-Hochberg FDR is applied within one BAG PRS within one analysis sample,
over every valid test in that family. Families are never pooled across PRS or
across samples, and the correction is computed on the full set of valid tests,
never on a pre-selected subset.

The PRS is raw SCORE1_SUM summed over chromosomes 1-22, with no standardisation
and no sign flip.

Usage -- one invocation per sample:

    set -a; . gwas/config/paths.env; set +a

    # primary: White British non-imaging
    ST18_SAMPLE=ukb_nonimaging_wb \
    PRS_PHEWAS_TABLE=$WORK/phewas/ukb_nonimaging.csv \
    PRS_PHEWAS_KEEP=$WORK/phewas/white_british.ids \
    ST18_EXPECTED_N=389045 \
    ST18_OUT_DIR=$WORK/tables/ukb_nonimaging_wb \
      python gwas/06_prs/05_analyze_ST18_phewas.py

    # sensitivity: the same table with no ancestry restriction
    ST18_SAMPLE=ukb_nonimaging \
    PRS_PHEWAS_TABLE=$WORK/phewas/ukb_nonimaging.csv \
    ST18_ALLOW_UNRESTRICTED=1 \
    ST18_EXPECTED_N=442295 \
    ST18_OUT_DIR=$WORK/tables/ukb_nonimaging \
      python gwas/06_prs/05_analyze_ST18_phewas.py

    # imaging held-out
    ST18_SAMPLE=ukb_imaging_heldout_10267 \
    PRS_PHEWAS_TABLE=$WORK/phewas/ukb_imaging_heldout_10267.csv \
    ST18_ALLOW_UNRESTRICTED=1 \
    ST18_EXPECTED_N=10267 \
    ST18_OUT_DIR=$WORK/tables/ukb_imaging_heldout_10267 \
      python gwas/06_prs/05_analyze_ST18_phewas.py
"""
import os
import sys
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from _regions import LABEL, REGIONS, require

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))


def die(msg):
    raise SystemExit(f"ERROR: {msg}")


PRS_WORK_DIR, PHEWAS_TABLE, PHECODE_CATALOG = require(
    "PRS_WORK_DIR", "PRS_PHEWAS_TABLE", "PRS_PHECODE_CATALOG")

# --- which sample this run is, and where it writes -------------------------
# ST18_SAMPLE is mandatory: it labels every output row and, unless ST18_OUT_DIR
# is given explicitly, gives the run its own directory. Two samples can never
# land on the same tables/ST18_full.tsv.
SAMPLE = os.environ.get("ST18_SAMPLE", "")
if not SAMPLE:
    die("ST18_SAMPLE is not set. Name the analysis sample, e.g. "
        "ukb_nonimaging_wb / ukb_nonimaging / ukb_imaging_heldout_10267; it "
        "selects the output directory and labels every result row.")
OUT_DIR = os.environ.get("ST18_OUT_DIR", os.path.join(PRS_WORK_DIR, "tables", SAMPLE))
NPROC = int(os.environ.get("PRS_NPROC", "24"))

# --- analysis-sample restriction -------------------------------------------
# PRS_PHEWAS_KEEP holds one participant id per line and selects the ANALYSIS
# sample. It is unrelated to which participants were scored (03_score.sh) and to
# the derivation / held-out split of the weights.
#
# Omitting it is a decision, not a default: a primary analysis that loses its
# keep file would otherwise become the broader sensitivity analysis in silence.
# Either pass PRS_PHEWAS_KEEP, or -- when the table has already been restricted
# upstream -- state that with ST18_ALLOW_UNRESTRICTED=1. A keep file that is
# missing, empty or matches nothing is an input error and never falls back to
# the whole table.
KEEP_FILE = os.environ.get("PRS_PHEWAS_KEEP", "")
ALLOW_UNRESTRICTED = os.environ.get("ST18_ALLOW_UNRESTRICTED", "") == "1"
# Optional but recommended: the analysis sample size this run must produce.
EXPECTED_N = int(os.environ.get("ST18_EXPECTED_N", "0"))

# --- sex designations for the aggregated test units ------------------------
# PRS_PHECODE_SEX_MAP is a two-column TSV (phecode, sex) shipped with the
# release; it is the mapping the published run applied. PRS_PHECODE_SEX_DEFS,
# when set, is the phecode definitions table the map was derived from -- a table
# with `phecode` and `sex` columns, where a unit takes its sex from an exact key
# match if one exists and otherwise from unanimous agreement among its child
# codes. Setting it re-derives the mapping and STOPS on any disagreement with
# the shipped file; it never silently replaces one with the other.
SEX_MAP = os.environ.get("PRS_PHECODE_SEX_MAP",
                         os.path.join(HERE, "sex_specific_phecodes.tsv"))
SEX_DEFS = os.environ.get("PRS_PHECODE_SEX_DEFS", "")
# Units whose child codes disagree on sex: fitted unrestricted (model C).
AMBIGUOUS = [s.strip() for s in os.environ.get("PRS_PHECODE_SEX_AMBIGUOUS", "751.1").split(",")
             if s.strip()]

COVS_A = ["Age", "Sex"] + [f"PC{i}" for i in range(1, 11)]
COVS_B = ["Age"] + [f"PC{i}" for i in range(1, 11)]
SEXCODE = {"Male": 1, "Female": 2}


def norm_phecode(x):
    """Phecode as text: 401.0 -> '401', 250.20 -> '250.2'. Never a float."""
    s = str(x).strip()
    try:
        return f"{float(s):.2f}".rstrip("0").rstrip(".")
    except ValueError:
        return ""


def agg_key(col):
    """Test unit of a raw phecode column: XX174.11 -> XX174.1 ; XX185 -> XX185.

    Integer and one-decimal codes stay separate units; only the second decimal
    place is folded away.
    """
    s = col.replace("XX", "")
    return f"XX{s.split('.')[0]}.{s.split('.')[1][:1]}" if "." in s else col


def prs(region):
    total = None
    for i in range(1, 23):
        s = pd.read_csv(f"{PRS_WORK_DIR}/scores_ukb/{region}/chr{i}.sscore",
                        sep="\t")[["IID", "SCORE1_SUM"]].rename(columns={"SCORE1_SUM": f"c{i}"})
        total = s if total is None else total.merge(s, on="IID", how="inner")
    # str on both sides: the phecode table's id column and plink2's IID are read
    # with different dtypes, and an int64/object merge raises.
    return pd.DataFrame({"eid": total["IID"].astype(str),
                         f"PRS_{region}": total[[f"c{i}" for i in range(1, 23)]].sum(axis=1)})


def read_sex_map(path):
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")
    for col in ("phecode", "sex"):
        if col not in t.columns:
            die(f"{path} has no `{col}` column")
    out = {}
    for r in t.itertuples():
        key = "XX" + norm_phecode(r.phecode)
        if key == "XX":
            die(f"{path}: unparseable phecode {r.phecode!r}")
        if r.sex not in SEXCODE:
            die(f"{path}: {r.phecode} has sex {r.sex!r}; expected one of {sorted(SEXCODE)}")
        if key in out and out[key] != r.sex:
            die(f"{path}: {r.phecode} appears twice with conflicting sex "
                f"({out[key]} vs {r.sex})")
        out[key] = r.sex
    return out


def derive_sex_map(path, units):
    defs = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    for col in ("phecode", "sex"):
        if col not in defs.columns:
            die(f"{path} has no `{col}` column")
    defs = defs[defs["sex"].isin(SEXCODE)].copy()
    defs["key"] = defs["phecode"].map(norm_phecode)
    conflict = defs.groupby("key")["sex"].nunique()
    conflict = sorted(conflict[conflict > 1].index)
    if conflict:
        die(f"{path}: these phecodes carry conflicting sex labels: {', '.join(conflict)}")
    sexdef = defs.drop_duplicates("key").set_index("key")["sex"]
    out = {}
    for unit, children in units.items():
        k = norm_phecode(unit.replace("XX", ""))
        if k in sexdef.index:
            out[unit] = sexdef[k]
            continue
        cs = {sexdef[norm_phecode(c.replace("XX", ""))] for c in children
              if norm_phecode(c.replace("XX", "")) in sexdef.index}
        if len(cs) == 1:
            out[unit] = cs.pop()
    return out


def restrict(mapping, units):
    """Both sides of a comparison must cover the same units and drop the same
    ambiguous ones before they can be compared at all."""
    return {u: s for u, s in mapping.items()
            if u in units and norm_phecode(u.replace("XX", "")) not in AMBIGUOUS}


def resolve_sex_map(units):
    shipped = restrict(read_sex_map(SEX_MAP), units)
    if not SEX_DEFS:
        print(f"sex map: {SEX_MAP} ({len(shipped)} units applicable to this table)", flush=True)
        return shipped
    derived = restrict(derive_sex_map(SEX_DEFS, units), units)
    disagree = sorted(u for u in set(derived) | set(shipped) if derived.get(u) != shipped.get(u))
    print(f"sex map: shipped {SEX_MAP} ({len(shipped)}) vs derived from {SEX_DEFS} "
          f"({len(derived)}); disagreements {len(disagree)}", flush=True)
    if disagree:
        for u in disagree:
            print(f"  {u}: shipped={shipped.get(u)} derived={derived.get(u)}", flush=True)
        die("the shipped sex map and the definitions table disagree. Resolve the "
            "difference before fitting; this run does not choose between them.")
    return shipped


# --- load ------------------------------------------------------------------
print(f"[{SAMPLE}] loading the phenotype table ...", flush=True)
df = pd.read_csv(PHEWAS_TABLE, low_memory=False)
print(f"  shape {df.shape}", flush=True)
df["eid"] = df["eid"].astype(str)
if not df["eid"].is_unique:
    die("the phecode table has duplicate participant ids")

if KEEP_FILE:
    if not os.path.exists(KEEP_FILE):
        die(f"PRS_PHEWAS_KEEP is set but does not exist: {KEEP_FILE}")
    keep = {ln.strip() for ln in open(KEEP_FILE) if ln.strip()}
    if not keep:
        die(f"PRS_PHEWAS_KEEP is empty: {KEEP_FILE}")
    matched = df["eid"].isin(keep)
    if not matched.any():
        die(f"PRS_PHEWAS_KEEP matched no participant in the phecode table: {KEEP_FILE}")
    df = df[matched].reset_index(drop=True)
    restriction = f"PRS_PHEWAS_KEEP={KEEP_FILE}"
    print(f"  analysis sample restricted: {len(df):,} rows, "
          f"{len(keep):,} ids listed", flush=True)
elif ALLOW_UNRESTRICTED:
    restriction = "none (ST18_ALLOW_UNRESTRICTED=1)"
    print(f"  no keep file; the whole table is the analysis sample: {len(df):,} rows", flush=True)
else:
    die("no PRS_PHEWAS_KEEP and no ST18_ALLOW_UNRESTRICTED=1. Pass the keep file "
        "for a restricted analysis sample, or set ST18_ALLOW_UNRESTRICTED=1 to "
        "state that PRS_PHEWAS_TABLE is already the intended sample. A missing "
        "keep file must not turn a restricted analysis into an unrestricted one.")

if EXPECTED_N and len(df) != EXPECTED_N:
    die(f"analysis sample is {len(df):,} rows, ST18_EXPECTED_N says {EXPECTED_N:,}. "
        "Check PRS_PHEWAS_TABLE and PRS_PHEWAS_KEEP; participants are never added "
        "or dropped to reach a target count.")

# Refuse to overwrite a previous run's results.
MAIN_OUT = os.path.join(OUT_DIR, "ST18_full.tsv")
if os.path.exists(MAIN_OUT) and os.environ.get("ST18_OVERWRITE", "") != "1":
    die(f"{MAIN_OUT} already exists. Point ST18_OUT_DIR at a new directory, or set "
        "ST18_OVERWRITE=1 to replace it deliberately.")

raw = [c for c in df.columns if c.startswith("XX")]
for c in raw:
    # NaN = phecode exclusion range; kept as NaN and dropped per test below
    df[c] = df[c].astype(np.float32)

for region in REGIONS:
    p = prs(region)
    if f"PRS_{region}" in df.columns:
        df = df.drop(columns=[f"PRS_{region}"])
    df = df.merge(p, on="eid", how="left")
    print(f"  PRS_{region}: matched {df[f'PRS_{region}'].notna().sum():,}/{len(df):,}", flush=True)

units = {}
for col in raw:
    units.setdefault(agg_key(col), []).append(col)
sub = df[raw]
# skipna maximum: a case at any child makes the unit a case; [0, NaN] -> 0;
# only an all-NaN row stays NaN and is dropped from the test below.
agg = pd.DataFrame({k: (sub[v[0]] if len(v) == 1 else sub[v].max(axis=1))
                    for k, v in units.items()}, index=df.index)

UNIT_SEX = resolve_sex_map(units)
print(f"phecodes: {len(raw)} raw -> {len(units)} aggregated units; sex-specific "
      f"{len(UNIT_SEX)} (Female {sum(v == 'Female' for v in UNIT_SEX.values())}, "
      f"Male {sum(v == 'Male' for v in UNIT_SEX.values())}); "
      f"fitted unrestricted despite conflicting child designations: "
      f"{','.join(AMBIGUOUS) or 'none'}", flush=True)

SEX = df["Sex"].to_numpy()
COV = {"A": df[COVS_A].to_numpy(np.float64), "B": df[COVS_B].to_numpy(np.float64)}
PRS = {r: df[f"PRS_{r}"].to_numpy(np.float64) for r in REGIONS}
Y = {u: agg[u].to_numpy(np.float32) for u in units}
N_ANALYSIS = len(df)
del df, sub, agg


def one(task):
    region, unit = task
    pc = unit.replace("XX", "")
    sex_rule = UNIT_SEX.get(unit, "")
    model = "B" if sex_rule else ("C" if pc in AMBIGUOUS else "A")
    cov = COV["B"] if model == "B" else COV["A"]
    base = (SEX == SEXCODE[sex_rule]) if model == "B" else np.ones(len(SEX), bool)

    y, x = Y[unit], PRS[region]
    ok = base & np.isfinite(x) & np.isfinite(cov).all(axis=1) & (~np.isnan(y))
    y2 = y[ok].astype(np.int8)
    n, n_case = int(ok.sum()), int(y2.sum())
    rec = {"sample": SAMPLE, "Region": LABEL[region], "region": region, "phecode": pc,
           "model": model, "sex_rule": sex_rule, "sex_restricted": model == "B",
           "n_total": n, "n_case": n_case, "n_control": n - n_case,
           # The two executed outputs spelled the control count differently
           # (N_control / n_control); both are written so either consumer works.
           "N_control": n - n_case}

    if n_case < 30:
        return dict(rec, status="excluded", reason="fewer than 30 cases")
    if len(np.unique(y2)) < 2:
        return dict(rec, status="excluded", reason="single outcome class")
    X = np.column_stack([np.ones(n), x[ok], cov[ok]])
    try:
        m = sm.Logit(y2, X).fit(disp=0, maxiter=50)
    except Exception as e:
        return dict(rec, status="excluded", reason=f"fit raised {type(e).__name__}")
    if not m.mle_retvals["converged"]:
        return dict(rec, status="excluded", reason="did not converge at maxiter=50")
    b, se, p = m.params[1], m.bse[1], m.pvalues[1]
    # A rank-deficient fit can still report converged=True; its bse and p_value
    # come back NaN and, left in, would blank out the whole BH family.
    if not (np.isfinite(b) and np.isfinite(se) and np.isfinite(p)):
        return dict(rec, status="excluded", reason="non-finite coefficient, SE or P")
    lo, hi = m.conf_int()[1]
    return dict(rec, status="fitted", reason="",
                coef=float(b), se=float(se), odds_ratio=float(np.exp(b)),
                ci95_lower=float(lo), ci95_upper=float(hi),
                or_ci95_lower=float(np.exp(lo)), or_ci95_upper=float(np.exp(hi)),
                p_value=float(p))


if __name__ == "__main__":
    tasks = [(r, u) for r in REGIONS for u in sorted(units)]
    print(f"attempting {len(tasks):,} logistic regressions "
          f"({len(units):,} units x {len(REGIONS)} BAG PRS) ...", flush=True)
    with Pool(NPROC) as pool:
        status = list(pool.imap_unordered(one, tasks, chunksize=8))

    ST = pd.DataFrame(status).sort_values(["Region", "phecode"]).reset_index(drop=True)
    out = ST[ST.status == "fitted"].drop(columns=["status", "reason"]).reset_index(drop=True)
    print(f"  attempted {len(ST):,} -> fitted {len(out):,}, "
          f"excluded {len(ST) - len(out):,}", flush=True)
    print("  exclusion reasons: "
          + str(ST[ST.status == "excluded"].reason.value_counts().to_dict()), flush=True)
    print("  by model: "
          + str(ST.groupby(["model", "status"]).size().unstack(fill_value=0).to_dict()), flush=True)

    cat = pd.read_csv(PHECODE_CATALOG)
    cat["phecode"] = cat["Phecode"].map(norm_phecode)
    out = out.merge(cat[["phecode", "Phenotype", "Category"]].drop_duplicates("phecode"),
                    on="phecode", how="left")
    if out.duplicated(subset=["Region", "phecode"]).any():
        die("duplicate (Region, phecode) among the fitted tests")

    # BH within one BAG PRS, over every valid test in this sample.
    out["fdr"] = np.nan
    for region, g in out.groupby("Region"):
        out.loc[g.index, "fdr"] = multipletests(g["p_value"].to_numpy(float),
                                                method="fdr_bh")[1]
    out = out.sort_values(["Region", "p_value"], kind="mergesort").reset_index(drop=True)

    sig = out[out.fdr < 0.05]
    summary = pd.DataFrame([
        {"sample": SAMPLE, "Region": region, "n_valid_tests": len(g),
         "n_sig": int((g.fdr < 0.05).sum()),
         "n_sig_positive": int(((g.fdr < 0.05) & (g.coef > 0)).sum()),
         "n_sig_negative": int(((g.fdr < 0.05) & (g.coef < 0)).sum()),
         "n_distinct_sig_phecodes": int(g.loc[g.fdr < 0.05, "phecode"].nunique())}
        for region, g in out.groupby("Region")])

    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(MAIN_OUT, sep="\t", index=False)
    sig.to_csv(f"{OUT_DIR}/ST18.tsv", sep="\t", index=False)
    ST.to_csv(f"{OUT_DIR}/ST18_fit_status.tsv", sep="\t", index=False)
    summary.to_csv(f"{OUT_DIR}/ST18_summary_by_prs.tsv", sep="\t", index=False)

    print(f"\n=== ST18 [{SAMPLE}] ===")
    print(f"  analysis sample  : {N_ANALYSIS:,}  (restriction: {restriction})")
    print(f"  valid tests      : {len(out):,}")
    print(f"  FDR<0.05 rows    : {len(sig)}")
    print(f"  distinct phecodes: {sig.phecode.nunique()}")
    print(f"  positive coef    : {(sig.coef > 0).sum()}/{len(sig)}")
    print(f"wrote {OUT_DIR}/ST18.tsv, ST18_full.tsv, ST18_fit_status.tsv, "
          f"ST18_summary_by_prs.tsv")
