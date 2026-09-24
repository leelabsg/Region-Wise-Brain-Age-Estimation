#!/usr/bin/env python3
"""Update the existing ST13 table IN PLACE with the ProtAge-204 correlations and a revised
external-aging genetic-correlation FDR family.

Genetic FDR families:
  - inter-regional (45 tests): UNCHANGED (already BH over 45 in the current table).
  - external-aging (21 tests, BH over all 21 together):
        Telomere Length PC1 x 10 BAG  +  Protein Age x 10 BAG  +  Protein Age x Telomere Length PC1.
    Genetic_fdr is RECALCULATED and replaced for all 21 rows (Genetic_p stays = log P).
Row changes:
  - 10 'Protein Age' x BAG rows: Genetic_r/Z/p are KEPT as build_rg_table.py estimated them from this
    run's own LDSC --rg logs. Setting PROTAGE_SUBSTITUTE_GENETIC=1 instead replaces all ten from the
    pre-computed PROTEIN_AGE_RG_TSV; that is opt-in because it discards the run's own estimates, and
    the published table was built with it OFF. Phenotypic_r/p/q are taken from
    PROTEIN_AGE_PHENO_CORR_TSV, the ProtAge-204 phenotypic correlation output (BH over 10, unchanged
    family); build_rg_table.py's inline phenotypic columns come from PROTEIN_AGE_PHENO_FILE, which is
    a different input and is not the 204-panel result.
  - 1 'Protein Age' x 'Telomere Length PC1' row appended (Genetic from ldsc_verified; Phenotypic NA);
  - 10 'Telomere Length PC1' x BAG rows: Genetic_r/Z/p and Phenotypic (NA) kept; only Genetic_fdr
    recomputed within the 21-family.
Phenotypic FDR families unchanged: inter-regional x45; Protein Age x BAG x10; Protein x Telomere NA.
Atomic in-place write (temp + os.replace). NOTE: the target ST13 file is NOT git-tracked, so the
launcher makes a one-time timestamped backup before this script overwrites it.
"""
import os, csv
import numpy as np



def env(name):
    v = os.environ.get(name)
    if not v:
        raise SystemExit(f"ERROR: {name} is not set (see gwas/config/paths.env.example)")
    return v


ST13 = env("RG_TABLE_OUT")                    # written by build_rg_table.py; updated in place
# Supplies the Protein Age x Telomere Length PC1 row, which this run does not estimate and the table
# has no other source for. It ALSO carries ten Protein Age x BAG rows; those are used only when
# PROTAGE_SUBSTITUTE_GENETIC=1.
LDSC_VERIFIED = env("PROTEIN_AGE_RG_TSV")
PHENO = env("PROTEIN_AGE_PHENO_CORR_TSV")     # phenotypic r/p/q for Protein Age x BAG
# Opt-in. OFF (the default, and how the published table was built) keeps the Protein Age x BAG
# genetic estimates that build_rg_table.py read from this run's own LDSC logs.
SUBSTITUTE_GENETIC = os.environ.get("PROTAGE_SUBSTITUTE_GENETIC", "") == "1"

REGIONS = ["global","caudate","cerebellum","frontal_lobe","insula",
           "occipital_lobe","parietal_lobe","putamen","temporal_lobe","thalamus"]
DISPLAY = {"global":"Whole Brain","caudate":"Caudate","cerebellum":"Cerebellum",
           "frontal_lobe":"Frontal Lobe","insula":"Insula","occipital_lobe":"Occipital Lobe",
           "parietal_lobe":"Parietal Lobe","putamen":"Putamen","temporal_lobe":"Temporal Lobe",
           "thalamus":"Thalamus"}
TELDISP = os.environ.get("TELOMERE_DISPLAY", "Telomere Length PC1")
# Display name for the proteomic-age trait, e.g. the specific model used.
PROTDISP = os.environ.get("PROTEIN_AGE_DISPLAY", "Protein Age")
UNDERFLOW = 5e-324   # table convention for '<2.2e-308'

def pfloat(s):
    s = s.strip()
    if s.startswith("<") or s in ("0", "0.0"):
        return UNDERFLOW
    return float(s)

def bh_fdr(pvals):
    p = np.asarray(pvals, float); m = len(p)
    order = np.argsort(p); q = np.empty(m); prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank-1]; prev = min(prev, p[i]*m/rank); q[i] = min(prev, 1.0)
    return q

def fmtq(q):
    return "%.4g" % q

# --- Genetic (Protein rows): rg, Z, P_log from ldsc_verified ---
gen = {}
with open(LDSC_VERIFIED) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        gen[row["Trait2"]] = dict(rg=row["rg"], Z=row["Z"], p=row["P_log"])
assert "telomere" in gen, "missing genetic row: telomere (Protein Age x Telomere Length PC1)"
if SUBSTITUTE_GENETIC:
    for k in REGIONS:
        assert k in gen, f"missing genetic row: {k}"
print("Protein Age x BAG genetic rows: "
      + (f"REPLACED from {LDSC_VERIFIED} (PROTAGE_SUBSTITUTE_GENETIC=1)" if SUBSTITUTE_GENETIC
         else "kept from this run's LDSC logs as build_rg_table.py wrote them"))
print(f"Protein Age x BAG phenotypic rows: from {PHENO}")

# --- Phenotypic (Protein x BAG): unadjusted INT-INT Pearson (r, p, q over 10) ---
phe = {}
with open(PHENO) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["family"] == "primary_unadjusted_int_ST13compat":
            phe[row["trait"]] = dict(r=row["r"], p=row["p_value"], q=row["q_bh10"])
for r in REGIONS:
    assert r in phe, f"missing phenotypic row: {r}"

# Phenotypic values for the Protein x BAG rows, keyed by display name. The genetic
# columns are filled per row below: kept from the table, or replaced when opted in.
new_phe = {DISPLAY[r]: phe[r] for r in REGIONS}
new_gen = {DISPLAY[r]: gen[r] for r in REGIONS} if SUBSTITUTE_GENETIC else {}
tel = gen["telomere"]
tel_row = [PROTDISP, TELDISP, tel["rg"], tel["Z"], tel["p"], None, "NA", "NA", "NA"]

# --- read + transform ST13 ---
with open(ST13) as f:
    lines = f.read().splitlines()

out, header_seen, replaced = [], False, set()
tel_bag_rows = []   # references into `out` for Telomere x BAG rows (to set fdr)
prot_rows = []      # references for Protein x BAG rows
genfdr_c, src_c = False, False
for ln in lines:
    if ln.startswith("#"):
        if ln.startswith("# Genetic_fdr:"):
            out.append("# Genetic_fdr: Benjamini-Hochberg within each family separately "
                       "(inter-regional 45; external-aging 21 = Telomere Length PC1 x BAG 10 "
                       "+ Protein Age x BAG 10 + Protein Age x Telomere Length PC1 1).")
            genfdr_c = True
        elif "Phenotypic source (Protein Age)" in ln:
            out.append("# Phenotypic source (Protein Age): ProtAge-204 age-bias-corrected proteomic-age gap "
                       "(protage204_gap_int) x {region}_corrected_delta_age_int, overlap N=5350; "
                       "unadjusted Pearson on inverse-normal-transformed phenotypes; BH-FDR over 10 BAG.")
            out.append("# Protein Age model = ProtAge-204 (Argentieri 2024 replication; repo commit eda0b4e). "
                       "Genetic_p = authoritative log P. 'Protein Age x Telomere Length PC1' Phenotypic = NA.")
            src_c = True
        else:
            out.append(ln)
        continue
    if ln.strip() == "":
        out.append(ln)
        continue
    cols = ln.split("\t")
    if not header_seen:
        assert cols[:9] == ["Trait1","Trait2","Genetic_r","Genetic_Z","Genetic_p","Genetic_fdr",
                            "Phenotypic_r","Phenotypic_p","Phenotypic_fdr"], f"unexpected header: {cols}"
        out.append(cols); header_seen = True; continue
    if cols[0] in ("Protein Age", PROTDISP):
        disp = cols[1]
        assert disp in new_phe, f"no phenotypic values for Protein Age '{disp}'"
        p = new_phe[disp]
        if SUBSTITUTE_GENETIC:
            g = new_gen[disp]
            row = [PROTDISP, disp, g["rg"], g["Z"], g["p"], None, p["r"], p["p"], p["q"]]
        else:
            # keep Genetic_r / Genetic_Z / Genetic_p from the table; Genetic_fdr is
            # recomputed below within the 21-test external-aging family.
            row = [PROTDISP, disp, cols[2], cols[3], cols[4], None, p["r"], p["p"], p["q"]]
        out.append(row); prot_rows.append(row); replaced.add(disp)
    elif cols[0] == TELDISP:
        out.append(cols); tel_bag_rows.append(cols)
    else:
        out.append(cols)   # inter-regional, unchanged

assert genfdr_c, "Genetic_fdr comment not updated"
assert src_c, "Protein Age source comment not updated"
assert replaced == set(new_phe.keys()), f"updated {sorted(replaced)} != {sorted(new_phe.keys())}"
assert len(replaced) == 10, f"updated {len(replaced)} Protein x BAG rows, expected 10"
assert len(tel_bag_rows) == 10, f"found {len(tel_bag_rows)} Telomere x BAG rows, expected 10"
assert out[-1] is prot_rows[-1] if prot_rows else False, "Protein Age block not last; refusing to append"
out.append(tel_row); prot_rows.append(tel_row)

# --- external-aging 21-family BH-FDR over log P (Telomere10 + Protein10 + Protein-telomere1) ---
ext_rows = tel_bag_rows + prot_rows          # 10 + 11 = 21
assert len(ext_rows) == 21, f"external-aging family = {len(ext_rows)}, expected 21"
ext_p = [pfloat(r[4]) for r in ext_rows]     # col 4 = Genetic_p
ext_q = bh_fdr(ext_p)
for r, q in zip(ext_rows, ext_q):
    r[5] = fmtq(q)                           # col 5 = Genetic_fdr

# --- validate row counts: 45 inter + 10 tel-BAG + 10 prot-BAG + 1 prot-tel = 66 ---
data = [r for r in out if isinstance(r, list) and r[0] != "Trait1"]
n_tel = sum(1 for r in data if r[0] == TELDISP)
n_prot = sum(1 for r in data if r[0] == PROTDISP)
n_inter = len(data) - n_tel - n_prot
assert len(data) == 66, f"expected 66 data rows, got {len(data)}"
assert n_inter == 45 and n_tel == 10 and n_prot == 11, f"blocks: inter={n_inter} tel={n_tel} prot={n_prot}"

# --- atomic in-place write ---
tmp = ST13 + ".tmp.%d" % os.getpid()
with open(tmp, "w") as f:
    for item in out:
        f.write((item if isinstance(item, str) else "\t".join(item)) + "\n")
os.replace(tmp, ST13)

print(f"updated in place: {ST13}")
print(f"data rows = {len(data)} (inter-regional {n_inter}, Telomere x BAG {n_tel}, Protein Age {n_prot})")
print(f"external-aging 21-family BH applied to {len(ext_rows)} rows")
print("\n=== external-aging rows (Trait1  Trait2  Genetic_r  Genetic_p  Genetic_fdr) ===")
for r in ext_rows:
    print(f"{r[0]:<20} {r[1]:<18} rg={r[2]:<9} p={r[4]:<12} fdr={r[5]}")
print("DONE.")
