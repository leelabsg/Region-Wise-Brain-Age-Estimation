# Genetic analysis of region-wise brain-age gaps

Everything downstream of the brain-age model: the GWAS of each regional
brain-age gap (BAG) in UK Biobank, its replication in ADNI, LD score regression
(SNP heritability and genetic correlation), and MAGMA cell-type enrichment.

<img width="4320" height="3060" alt="fig1_bcd_overview_preview" src="https://github.com/user-attachments/assets/39782365-1ab4-403e-a6dd-a337461ef952" />


The model that produces the phenotype lives in the repository root; this
directory starts from the per-region BAG table it writes.

Paths, host names and container references were lifted out into
`config/paths.env`. Analysis settings were not touched: every SAIGE, PLINK,
LDSC and MAGMA flag, the bias-correction formula, the point at which the
inverse-normal transform is applied, the covariate lists, and the cell-type
retention rule are byte-for-byte what produced the published numbers.

## Setup

```bash
cp gwas/config/paths.env.example gwas/config/paths.env
$EDITOR gwas/config/paths.env          # every site-specific path lives here
```

Shell scripts read it through `config/common.sh`. Python scripts read the same
variables from the environment:

```bash
set -a; . gwas/config/paths.env; set +a
```

`config/regions.tsv` maps region codes to display names and records which cohort
each region was run in.

## Order of execution

```
  RegionBAE model  ──►  per-region brain-age predictions
          │
          ▼
  01_ukb_gwas          SAIGE step 1 → step 2 → merge → FUMA input
          │
          ├──────────────► FUMA SNP2GENE (web)  ──►  magma.genes.raw per region
          │                                              │
          ▼                                              ▼
  02_adni_replication                          05_magma_celltype
  build_pheno → GRM → BGEN → step 1 → step 2     cell-type enrichment
          │
          ▼
  03_postgwas          λGC / QQ, Manhattan, UKB→ADNI concordance

  04_ldsc              munge → merge → h² and rg      (runs off 01's sumstats)
                       external traits: 00_preprocess_telomere.py (telomere) or
                       07_protage/11_ldsc_downstream.py (ProtAge-204) → 00_munge_external.sh
                       (telomere) → 03_merge_munged.py, then the same h² and rg steps

  06_prs               SBayesRC weights → scoring → ST17 / ST18 / ST21-23
                       (runs off 01's TRAINING-split sumstats)

  07_protage           ProtAge-204 proteomic-age gap: protein QC → KNN imputation →
                       LightGBM clock (out-of-fold) → gap → SAIGE GWAS → LDSC h²,
                       the 21-test external-ageing rg family and the phenotypic
                       correlations with the BAGs
```

`05_magma_celltype` cannot start until the FUMA SNP2GENE job has finished, and
`03_postgwas`'s concordance step needs both cohorts' summary statistics.

`06_prs` is the one branch that must not see the full sample: its weights are
fitted on the training split only, so the held-out validation in ST17 stays held
out. SBayesRC is an MCMC sampler, and multi-threaded runs from identical inputs
and seed did not always return identical weights in the whole-brain
reproducibility test; the weights are deposited so that the reported scores can
be reproduced from them. `06_prs/README.md` gives the scope of that test.

## What each stage needs

| stage | needs |
|---|---|
| `01_ukb_gwas` | UK Biobank imputed BGEN + `.bgi` + sample file, an LD-pruned PLINK1 set for the GRM, and the per-region BAG phenotype table. Approved UK Biobank application. |
| `02_adni_replication` | ADNI imputed dosage pgen, the QC sample list, in-sample PCA eigenvectors, genotyping-batch indicators, and the model's ADNI predictions. Approved ADNI data access. |
| `03_postgwas` | merged summary statistics from both cohorts. |
| `04_ldsc` | a clone of [bulik/ldsc](https://github.com/bulik/ldsc) and its EUR LD-score panels (`eur_ref_ld_chr`, `eur_w_ld_chr`). The two external ageing traits are built here as well, from `LDSC_TELOMERE_RAW_SUMSTATS` (GWAS Catalog accession GCST90435144) and `LDSC_PROTEIN_AGE_RAW_SUMSTATS` (the merged SAIGE output for the proteomic-age gap); neither file is redistributed here. The phenotypic-correlation columns additionally need individual-level BAG tables. |
| `07_protage` | UK Biobank Olink proteomics (instance 0) with plate/batch fields, birth year/month and assessment date, the published ProtAge-204 model repository (protein panel and Olink names), the same genotype files as `01_ukb_gwas`, and the BAG tables for the phenotypic correlations. Approved UK Biobank application. |
| `05_magma_celltype` | MAGMA v1.10, the 1000 Genomes phase 3 EUR reference (GRCh37) with an ENSG gene annotation, FUMA's pre-processed single-cell expression matrices, FUMA's `magma_celltype.R`, and one `magma.genes.raw` per region from a FUMA SNP2GENE job. |

### Phenotype table

One row per participant, `IID` plus the covariates, plus one column per region
named `{region}_corrected_delta_age_int` — the bias-corrected brain-age gap
after a rank-based inverse-normal transform. UK Biobank covariates are
`Sex, Age, PC1..PC10`; ADNI adds genotyping-batch dummies, because ADNI was
genotyped in several waves.

`02_adni_replication/build_pheno.py` builds the ADNI table from the model's
predictions. The order is fixed and matters: fit the age-bias regression on
cognitively normal baseline scans, apply it to all baseline scans, restrict to
the final genotyped sample, and only then inverse-normal transform, per region.
Fitting the correction on the analysis sample, or transforming before
restricting, changes the phenotype.

### FUMA dependency

Two distinct things come from FUMA and neither is redistributed here.

1. **Gene analysis.** `magma.genes.raw` per region comes from a FUMA SNP2GENE
   job (FUMA v1.8.2; the job's own MAGMA is v1.10). Upload the files written by
   `01_ukb_gwas/make_fuma_input.sh`, then download the MAGMA results and point
   `FUMA_GENES_RAW_TMPL` at them. Running MAGMA's gene analysis locally is *not*
   equivalent — see [Local vs FUMA gene analysis](#local-vs-fuma-gene-analysis).
2. **Cell-type expression matrices and the cell-type script.** The
   pre-processed single-cell resource (`celltype/*.txt`, downloaded 2026-07-12;
   876 files, of which 679 are the human-brain datasets FUMA itself uses) and
   `storage/scripts/magma_celltype.R` from the
   [FUMA-webapp](https://github.com/Kyoko-wtnb/FUMA-webapp) repository. Steps 1,
   2 and 3 of the cell-type analysis are that script, run unmodified; the only
   part re-implemented here is step 3's pairwise loop, described under
   [The step-3 parallel driver](#the-step-3-parallel-driver).
   `05_magma_celltype/checks/` contains the comparisons that established our runs
   reproduce FUMA's own output exactly.

Record the FUMA version and access date for anything you re-run; FUMA's
pre-processed resource is updated over time and the dataset list is not stable
across versions.

## One release throughout

Every branch of this pipeline — LD score regression, the MAGMA cell-type
analysis and the ADNI replication — runs on the same UK Biobank release,
N = 41,067. `04_ldsc/02_munge.sh` still takes N as an explicit argument; pass
the N of whatever release you feed it.

## Local vs FUMA gene analysis

An early pass ran MAGMA's gene analysis locally against a full-Ensembl
annotation (49,329 genes) and the results were discarded: FUMA's annotation is
protein-coding only (19,011 genes), which changes `NGENES` (22,933 vs 16,777)
and flips significance for 21 cell types. The published
analysis uses FUMA's `magma.genes.raw`. `05_magma_celltype/01_gene_analysis_local.sh`
and `deprecated_pass1_gene_property.sh` are kept for the record and are not part
of the published pipeline; their headers say so.

## Effect-allele orientation in the LDSC branch

Worth understanding before adding a trait, because it is invisible in the output
and only shows up as a sign.

LDSC reads `A1` as the effect allele. SAIGE reports `BETA` and `AF_Allele2` for
`Allele2`, so `01_preprocess.py` writes `A1 = Allele2` into the merge-alleles
list and `02_munge.sh` passes `--a1 Allele2 --a2 Allele1`. The munged `Z`
therefore refers to the allele named in `A1`, which is what a downstream reader
assumes. `munge_sumstats.py` never re-signs `Z` — its `allele_merge` only drops
SNPs that fail to match the merge-alleles list — so this setting controls the
`A1` label and nothing else.

Earlier runs passed `--a1 Allele1`, which put the *other* allele in `A1` and
inverted every sign relative to the label. That did not corrupt any published
number. `validation/allele_convention_check.txt` is the measurement behind that
sentence: re-munging all 220 region-chromosome files of the N=41,067 input set
under the corrected convention leaves the SNP set identical in 220 of 220 files
and the `Z` column bit-identical (max |dZ| = 0 over 72,963,690 SNPs), while
moving `A1` onto the effect allele for 100% of rows. On chromosome 22 of the
whole-brain GWAS, where `--h2` was run both ways, heritability is 0.3566 (0.126)
either way. Because LDSC reads only `SNP`, `Z` and `N` from a munged file and
all three are unchanged, no h2 or rg estimate can move. The same inversion was
carried by every file this pipeline correlates — `pad_<region>`, the
proteomic-age munge, and the telomere munge (its `A1` came from
`other_allele`). The summary statistics released with the paper were produced
under the old convention and are left as they are; the change above applies to
any future run.

A separate matter, which does bite, is **allele order**. Whether a GWAS reports a
variant as `Allele1=C, Allele2=T` or the reverse is a property of the genotype
source, not of the phenotype, and it differs between runs:

| pair | allele order |
|---|---|
| `pad_<region>` (N=41,067) vs proteomic age (N=38,806) | same for 7,269,915 of 7,269,941 shared SNPs |
| `pad_<region>` (N=41,067) vs telomere | **opposite** for essentially every shared SNP |
| `pad_<region>` (N=45,076) vs telomere | same — this is why the original strict align worked |

The first aligner written for this step was a strict-match one: it kept only SNPs
whose `A1`/`A2` matched the reference exactly. On the released N=41,067 sumstats
that returns an empty telomere file rather than an error, because the allele order
is opposite for essentially every shared SNP. `align_flip.py` applies the same
restriction but flips `Z` on the order-swapped SNPs, which is the correct
alignment in both cases: run against the N=45,076 inputs it reproduces the
published `telomere_<region>_aligned.sumstats.gz` row for row with `max|dZ| = 0`,
and against the released inputs it produces a file with the same meaning. It is
the only aligner shipped here, and every external correlation in the correlation
table uses it.

The ProtAge-204 x telomere row of that table uses the same two steps as the ten
ProtAge-204 x BAG rows, with telomere length in place of a BAG phenotype
(`07_protage/11_ldsc_downstream.py`): `align_flip.py --source_file <telomere munged>
--reference_file <ProtAge-204 munged>`, then `ldsc.py --rg <ProtAge-204 munged>,<aligned telomere>`.
The published row (rg = -0.1886, Z = -4.0158, P = 5.92e-05) came from an alignment
reporting `shared=6815695 exact=0 swapped=6815692 other_dropped=3`.

Add a trait and the thing to check is its allele order against `pad_<region>`;
`align_flip.py` prints `exact=` and `swapped=` counts so the answer is in the log.

## Directory map

```
config/
  paths.env.example    every site-specific path, host and container reference
  common.sh            config loading, region list, SAIGE docker/conda dispatch
  regions.tsv          region code, display name, cohort

01_ukb_gwas/
  saige_step1.sh          null model per region
  saige_step2.sh          association per (region, chromosome), resumable
  check_step2_status.sh   per-chromosome completion report
  merge_chromosomes.sh    chr1-22 → one sumstats file, with completeness checks
  make_fuma_input.sh      column subset + gzip for FUMA upload
  variant_qc_filter.py    mfi INFO >= 0.8 and analysis-sample MAF >= 0.01; the
                          filter applied to the deposited summary statistics

02_adni_replication/
  build_pheno.py       CN bias correction → INT → covariate merge
  00_prep_grm.sh       LD-pruned hardcalls → PLINK1 GRM set
  01_convert_bgen.sh   dosage pgen → unphased BGEN v1.2, ref-first
  02_saige_step1.sh    null model per region
  03_saige_step2.sh    association per (region, chromosome), resumable

03_postgwas/
  lambda_qq.py              λGC and QQ panel
  manhattan_per_region.py   one Manhattan row per region
  manhattan_multi_region.py all regions overlaid on one panel
  concordance_audit.py      UKB→ADNI sign concordance on a common-SNP-pruned
                            lead set, with a lead-recovery diagnostic

04_ldsc/
  00_preprocess_telomere.py    external lane: GWAS Catalog telomere columns →
                               per-chromosome split
  00_munge_external.sh         external lane: per-trait --N, then 02_munge.sh
  01_preprocess.py             common-variant filter, per-chromosome split
  02_munge.sh                  munge_sumstats.py, N passed explicitly
  03_merge_munged.py           chr1-22 → one .sumstats.gz
  align_flip.py                restrict to a reference file's SNPs and alleles,
                               keeping A1/A2-swapped SNPs by flipping Z
  05_heritability.sh           ldsc.py --h2
  06_genetic_correlation.sh    ldsc.py --rg for one pair
  build_h2_table.py            h² table from the standalone --h2 logs
  h2_table_traits.tsv          display name → log path, in table order
  build_rg_table.py            correlation table, BAG-BAG family BH-FDR
                               (the external-ageing rows: 07_protage/11_ldsc_downstream.py)

06_prs/
  01_preprocess_gwas.py       SAIGE sumstats → GCTB .ma, with --check-orientation
  02_gctb_commands.sh         --impute-summary per LD block, merge, --sbayes RC
  03_score.sh                 plink2 --score, UK Biobank or ADNI
  04_analyze_ST17.py          held-out UKB validation, OLS BAG ~ PRS
  04a_build_adni_weights.py   rewrite the weights against ADNI variant IDs
  05_analyze_ST18_phewas.py   UKB PheWAS, one run per analysis sample
  sex_specific_phecodes.tsv   sex designations applied to sex-specific phecodes
  06_analyze_adni.py          ADNI correlation and case-control models
  07_export_weights.py        deposit package: weights + sha256 manifest
  08_phewas_additional_corrections.py  sample-wide and phecode-level BH (ST18 sensitivity)
  09_relatedness_prune.py     KING > 0.0884 pruning for the relatedness sensitivity
  _regions.py                 region codes and display names

05_magma_celltype/
  config.py                     paths, region list, model options
  inventory.py                  classify the expression matrices; human-brain only
  datasets_fuma_679.txt         the dataset list FUMA itself used
  02_celltype_step1.sh          per-dataset gene-property against FUMA's genes.raw
  03_celltype_step23.sh         runs FUMA's magma_celltype.R unmodified
  04_step3_parallel.py          parallel port of the step-3 pairwise loop
  cellgroup_map.py              cell-label harmonization engine
  build_mapping_table.py        harmonization table + rule audit
  build_manuscript_tables.py    step1/2/3 outputs → analysis tables
  mns_step_A_mapping_independent.py   significant-set table
  mns_step_B_retained_signals.py      conditional-retention table
  mns_step_C_harmonize.py             harmonized cell-group table
  cell_class_map.tsv            curated broad class per cell-type label
  05_build_supplementary_table.py  retained signals → the published table
  checks/                       comparisons against FUMA's own output
  reference/                    app.config template

07_protage/
  01_inputs.py            NPX pool, data-use exclusion, panel mapping, batch, granular age
  02_qc.py                assay / participant missingness QC (50% / 54% / 50%)
  03_knn_impute.R         impute::impute.knn, k = 10, rng.seed 362436069
  04_scale_cohort.py      MinMax + median centring, batch 0/7/NA removed, splits
  05_tune.py              Optuna TPE, 200 trials, LightGBM DART, 5-fold CV
  06_holdout_oof.py       30% hold-out and 5-fold out-of-fold predictions
  07_gap_phenotype.py     raw gap, age-residualised gap, rank-INT phenotype
  08_build_gwas_pheno.py  GWAS sample and phenotype/covariate file
  09_saige_step1.sh       SAIGE 1.1.9 null model
  10_saige_step2.sh       SAIGE 1.1.9 association, per chromosome
  11_ldsc_downstream.py   munge, h², rg with BAGs and telomere, 21-test BH,
                          phenotypic partial correlations

validation/
  allele_convention_check.txt   chr22 measurement behind the allele-labelling
                                note, with the commands to re-run it
```

### The step-3 parallel driver

For most regions FUMA's R script runs step 3 serially without trouble. Cerebellum
had 158 retained cell types spread over 158 datasets, i.e. C(158,2) = 12,403
cross-dataset pairs and roughly 83 hours of serial MAGMA. `04_step3_parallel.py`
ports that pairwise loop only — the forward selection is untouched — and was
validated to zero difference against the R output on five regions.

One detail there is load-bearing: R's `write.table` re-formats the covariate
matrix at 15 significant digits, and MAGMA's p-value for a collinear-sensitive
pair depends on that. The driver writes `%.15g` to match, which makes its
intermediate files byte-identical to R's.

## Notes on the published tables

`validation/allele_convention_check.txt` records the chromosome-22 comparison
behind the allele-labelling note above, and the commands to re-run it.

Three details about how the tables were produced:

- The h² printed inside an `--rg` log is a *different* estimate — LDSC re-fits
  it on the two-trait SNP intersection — and is not what the heritability table
  reports. `build_h2_table.py` refuses a log that is not a standalone `--h2`
  run.
- The ten brain rows of the heritability table went through a builder that
  fixed the number of decimals; the two external-trait rows were transcribed
  from their logs and kept the log's own precision. `h2_table_traits.tsv` marks
  those rows `verbatim`.
- The correlation table's FDR families are the inter-regional 45 and the
  external-aging 21 (`07_protage/11_ldsc_downstream.py`), computed separately. Where LDSC printed `P: 0.` the table
  shows `<2.2e-308` and the FDR used 5e-324, which preserves the ranking.

## What cannot be reproduced from this repository

No genotypes, imaging, phenotypes, participant identifiers, summary statistics or
model weights are included. UK Biobank and ADNI data are available to approved
researchers through their own access processes, and the analyses above cannot be
run without them.

Two further limits are worth stating plainly. The phenotypic-correlation columns
of the correlation table need individual-level brain-age gaps, so they cannot be
recomputed from summary statistics alone. And the broad cell class attached to
each significant cell type was assigned by hand rather than derived by rule,
which is why it ships as a lookup table (`cell_class_map.tsv`) instead of as
code.
