#!/bin/bash
# SAIGE 1.1.9 step 1 for the ProtAge-204 gap. Quantitative trait, precomputed rank-INT phenotype with
# --invNormalize=TRUE, covariates Age, Sex, PC1-PC10, LOCO, on the same LD-pruned PLINK set as the BAG GWAS.
#   set -a; . gwas/config/paths.env; set +a; bash gwas/07_protage/09_saige_step1.sh
set -euo pipefail
: "${SAIGE_ENV:?}" "${UKB_GRM_PLINK:?}" "${PROTAGE_GWAS_PHENO:?}" "${PROTAGE_GWAS_DIR:?}"
OUTP=$PROTAGE_GWAS_DIR/step1/protage204
mkdir -p "$(dirname "$OUTP")"
[ -e "$OUTP.rda" ] && { echo "refusing to overwrite $OUTP.rda" >&2; exit 1; }
ver=$("$SAIGE_ENV/bin/Rscript" -e 'cat(as.character(packageVersion("SAIGE")))')
[ "$ver" = "1.1.9" ] || { echo "SAIGE version $ver != 1.1.9" >&2; exit 2; }
"$SAIGE_ENV/bin/Rscript" "$SAIGE_ENV/bin/step1_fitNULLGLMM.R" \
  --plinkFile="$UKB_GRM_PLINK" \
  --phenoFile="$PROTAGE_GWAS_PHENO" \
  --phenoCol=protage204_gap_int \
  --covarColList=Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --sampleIDColinphenoFile=IID --traitType=quantitative --invNormalize=TRUE \
  --outputPrefix="$OUTP" --nThreads=4 --LOCO=TRUE --IsOverwriteVarianceRatioFile=TRUE
