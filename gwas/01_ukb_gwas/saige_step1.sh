#!/bin/bash
# SAIGE Step 1 (null model) -- UK Biobank region-wise brain-age-gap GWAS.
#
# One null model per region. The phenotype is the inverse-normal-transformed,
# bias-corrected brain-age gap produced by the RegionBAE model
# (column {region}_corrected_delta_age_int in $UKB_PHENO_FILE).
#
# Inverse normalization is retained to match the executed SAIGE command.
#
# Visit was included in the original covariate argument list but was constant in
# the final discovery and PRS-derivation samples and was excluded from the fitted
# models. It is retained here to match the original command. The fitted
# covariates were sex, age and the first ten genetic principal components.
#
# SAIGE version: step 1 of the published run used SAIGE 1.1.8 and step 2 used
# SAIGE 1.1.9. SAIGE_RUNNER_STEP1 / SAIGE_IMAGE_STEP1 / SAIGE_ENV_STEP1 select the
# step-1 build; each falls back to the corresponding SAIGE_* value when unset.
#
# Usage: bash saige_step1.sh [region ...]      (default: all UK Biobank regions)
set -uo pipefail
. "$(dirname "$0")/../config/common.sh"
gwas_require GWAS_UKB_DIR UKB_GRM_PLINK UKB_PHENO_FILE

COVARS="${UKB_COVARS:-Sex,Age,Visit,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10}"
SAIGE_RUNNER="${SAIGE_RUNNER_STEP1:-${SAIGE_RUNNER:-docker}}"
SAIGE_IMAGE="${SAIGE_IMAGE_STEP1:-${SAIGE_IMAGE:-}}"
SAIGE_ENV="${SAIGE_ENV_STEP1:-${SAIGE_ENV:-}}"

if [ "$#" -ge 1 ]; then regions=("$@"); else mapfile -t regions < <(gwas_regions ukb); fi

i=0
for region in "${regions[@]}"; do
    out_dir="$GWAS_UKB_DIR/$region"
    log_dir="$out_dir/logs"
    mkdir -p "$out_dir" "$log_dir"
    node="$(gwas_node $i)"

    echo "===== Step 1: $region ${node:+on $node} ====="
    cmd="$(gwas_saige_cmd step1_fitNULLGLMM.R \
        --plinkFile="$UKB_GRM_PLINK" \
        --phenoFile="$UKB_PHENO_FILE" \
        --phenoCol="${region}_corrected_delta_age_int" \
        --covarColList="$COVARS" \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --useSparseGRMtoFitNULL=FALSE \
        --outputPrefix="${out_dir}/step1" \
        --nThreads="${SAIGE_NTHREADS:-4}" \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE)"
    gwas_run "$node" "$cmd" > "$log_dir/step1_${region}.log" 2>&1 &
    i=$((i + 1))
done
wait
echo "All Step 1 jobs finished. Check $GWAS_UKB_DIR/*/logs/step1_*.log and */step1.rda"
