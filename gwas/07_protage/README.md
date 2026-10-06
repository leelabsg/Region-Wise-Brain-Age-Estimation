# ProtAge-204 proteomic-age gap

The external proteomic ageing measure correlated with the regional BAGs
(Supplementary Tables 14 and 15). The 204-protein panel and the model
specification follow the published ProtAge model (Argentieri et al., 2024;
its repository supplies the panel, the Olink names and the LightGBM helper
settings referred to in the code). Proteomic age is re-estimated here in UK
Biobank with out-of-fold predictions, and its age-bias-corrected gap is the
GWAS phenotype.

## Order

```
01_inputs.py          Olink NPX (instance 0), data-use exclusion, panel mapping, batch,
                      granular age = (assessment date - 1st of birth month) / 365.25
02_qc.py              assays >= 50% missing dropped; participants >= 54% missing dropped;
                      assays re-checked at 50%
03_knn_impute.R       impute::impute.knn(t(X), k = 10, rowmax = 0.5, colmax = 0.95,
                      maxp = ncol(X), rng.seed = 362436069); protein columns only
04_scale_cohort.py    per-protein MinMax + median centring on the imputed pool;
                      Olink batch 0, 7 and missing removed (clock cohort);
                      70/30 split and 5-fold folds, seed 3456
05_tune.py            Optuna TPE (seed 3456), 200 trials, LightGBM DART, 5,000 trees,
                      mean 5-fold R2 on the 70% development set
06_holdout_oof.py     one fit on 70% -> 30% hold-out; 5-fold out-of-fold prediction
                      for every clock-cohort participant
07_gap_phenotype.py   gap = OOF prediction - age; residualised on age; rank-INT
08_build_gwas_pheno.py  self-reported British, in the step-1 PLINK set and the BGEN sample
09_saige_step1.sh     SAIGE 1.1.9, Age, Sex, PC1-PC10, LOCO, --invNormalize=TRUE
10_saige_step2.sh     SAIGE 1.1.9, per chromosome, minMAC 10, alt-first (Allele2 = effect)
11_ldsc_downstream.py LDSC h2; rg with the ten BAGs and telomere length PC1; BH over the
                      21 external-ageing tests; phenotypic partial correlations
```

```bash
set -a; . gwas/config/paths.env; set +a
cd gwas/07_protage
python 01_inputs.py && python 02_qc.py && OMP_NUM_THREADS=1 Rscript 03_knn_impute.R
python 04_scale_cohort.py && python 05_tune.py && python 06_holdout_oof.py && python 07_gap_phenotype.py
python 08_build_gwas_pheno.py
bash 09_saige_step1.sh
for c in $(seq 1 22); do bash 10_saige_step2.sh $c; done     # then concatenate chr1-22 into PROTAGE_SUMSTATS
python 11_ldsc_downstream.py
```

The `--chrom` value in `10_saige_step2.sh` is written zero-padded, matching
the chromosome labels of the UK Biobank BGEN index files used here; adjust it
if your `.bgi` labels chromosomes differently.

## Published run

| step | N |
|---|---|
| clock cohort (out-of-fold) | 44,927 |
| self-reported British | 39,114 |
| in the step-1 PLINK set | 38,857 |
| in the imputed BGEN sample = GWAS sample (steps 1 and 2) | 38,806 |

Out-of-fold accuracy over the clock cohort: r = 0.936, R2 = 0.877,
MAE = 2.26 years. Quality control, imputation and scaling preceded data
partitioning, and tuning was not nested within the final out-of-fold
partitions. These performance estimates are therefore descriptive rather than
fully nested validation results. The GWAS tested 35,199,738 variants (MAC >= 10, GRCh37).
LDSC: h2 = 0.0985 (s.e. 0.0151) over 1,176,243 regression SNPs.

## Not included

The covariate source for Sex and the genetic principal components
(`PROTAGE_COVARIATES`) and all individual-level inputs require an approved
UK Biobank application. No participant-level output is shipped.
