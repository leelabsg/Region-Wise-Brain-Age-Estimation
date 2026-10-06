#!/bin/bash
# SAIGE 1.1.9 step 2 for one chromosome:  bash gwas/07_protage/10_saige_step2.sh <chr 1-22> [range-file]
# Imputed BGEN (alt-first allele order, so Allele2 is the effect allele), minMAF 0, minMAC 10, LOCO, Firth
# for P < 0.05. An optional --rangestoIncludeFile (tab-separated "chrom start end", chromosome written as in
# the .bgi) splits one chromosome into contiguous blocks; per-variant results are identical to an unsplit run
# and the blocks are concatenated in position order.
set -euo pipefail
: "${SAIGE_ENV:?}" "${UKB_BGEN_TMPL:?}" "${UKB_BGEN_INDEX_TMPL:?}" "${UKB_SAMPLE_FILE:?}" "${PROTAGE_GWAS_DIR:?}"
CHR=$1; RANGE=${2:-}
BGEN=${UKB_BGEN_TMPL//\{chr\}/$CHR}; BGI=${UKB_BGEN_INDEX_TMPL//\{chr\}/$CHR}
S1=$PROTAGE_GWAS_DIR/step1/protage204
OUT=$PROTAGE_GWAS_DIR/step2/chr$CHR${RANGE:+.$(basename "$RANGE" .range)}.txt
mkdir -p "$(dirname "$OUT")"
[ -e "$OUT" ] && { echo "refusing to overwrite $OUT" >&2; exit 1; }
"$SAIGE_ENV/bin/Rscript" "$SAIGE_ENV/bin/step2_SPAtests.R" \
  --bgenFile="$BGEN" --bgenFileIndex="$BGI" --sampleFile="$UKB_SAMPLE_FILE" --AlleleOrder=alt-first \
  --SAIGEOutputFile="$OUT" --chrom="$(printf '%02d' "$CHR")" --minMAF=0 --minMAC=10 --LOCO=TRUE \
  --GMMATmodelFile="$S1.rda" --varianceRatioFile="$S1.varianceRatio.txt" \
  --is_Firth_beta=TRUE --is_imputed_data=TRUE --pCutoffforFirth=0.05 --is_output_moreDetails=TRUE \
  ${RANGE:+--rangestoIncludeFile="$RANGE"}
