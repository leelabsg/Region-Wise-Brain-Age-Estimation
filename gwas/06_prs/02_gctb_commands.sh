#!/bin/bash
# SBayesRC weights for one region: .ma per chromosome -> LD-block imputation ->
# merge -> posterior sampling.
#
# This is the method, not our job scheduler. Every gctb flag below is what
# produced the published weights; the loops are plain and serial so that the
# parameters are readable. On our cluster the LD-block step was fanned out over
# many cores and the ten regions were run concurrently, which changes nothing
# about the result -- see README.md, "Reproducibility of the weights", before
# assuming a re-run reproduces our numbers bit for bit.
#
# Usage: bash 02_gctb_commands.sh <region>
set -uo pipefail
. "$(dirname "$0")/../config/common.sh"
gwas_require PRS_WORK_DIR PRS_SUMSTATS_DIR GCTB PLINK2 GCTB_LDM_EIGEN GCTB_ANNOT SBAYESRC_SEED

REGION="${1:?usage: bash 02_gctb_commands.sh <region>}"
PPR="$PRS_WORK_DIR/ma/$REGION"
IMP="$PRS_WORK_DIR/imputed/$REGION"
RUN="$PRS_WORK_DIR/weights/$REGION"
mkdir -p "$PPR" "$IMP" "$RUN"

# --- 1. SAIGE -> .ma, per chromosome --------------------------------------
# --check-orientation asserts A1 is the effect allele. Keep it on: an
# unnoticed flip here inverts every downstream sign and gctb will not complain.
for chr in $(seq 1 22); do
    python "$(dirname "$0")/01_preprocess_gwas.py" \
        --file_path "$PRS_SUMSTATS_DIR/$REGION/chr${chr}.txt" \
        --save_path "$PPR/chr${chr}.ma" \
        --check-orientation || exit 1
done

# --- 2. impute the summary statistics to the LD reference, block by block --
# One call per (chromosome, LD block); ldm.info lists the blocks and their
# chromosomes. The eigen-decomposed ukbEUR_HM3 reference has 591 blocks.
tail -n +2 "$GCTB_LDM_EIGEN/ldm.info" | awk -F'\t' '{print $1"\t"$2}' |
while IFS=$'\t' read -r block chr; do
    "$GCTB" --ldm-eigen "$GCTB_LDM_EIGEN" \
            --gwas-summary "$PPR/chr${chr}.ma" \
            --impute-summary --block "$block" --seed "$SBAYESRC_SEED" \
            --out "$IMP/chr${chr}_block${block}" \
            > "$IMP/chr${chr}_block${block}.log" 2>&1 || echo "IMPUTE FAIL $chr $block"
done

# --- 3. concatenate the imputed blocks into one .ma ------------------------
files=("$IMP"/*.imputed.ma)
[ "${#files[@]}" -eq 591 ] || { echo "ABORT: ${#files[@]}/591 imputed blocks"; exit 1; }
head -n 1 "${files[0]}" > "$IMP/merged.ma"
for f in "${files[@]}"; do tail -n +2 "$f" >> "$IMP/merged.ma"; done
echo "merged.ma rows=$(wc -l < "$IMP/merged.ma") md5=$(md5sum "$IMP/merged.ma" | cut -d' ' -f1)"

# --- 4. SBayesRC ----------------------------------------------------------
# Defaults not shown on the command line and not changed: chain length 3,000,
# burn-in 1,000, the baseline2.2 annotation set. In the whole-brain
# reproducibility test, two single-threaded runs with identical inputs and seed
# gave byte-identical weights, while four-thread runs did not (README.md). The
# deposited weights were produced at four threads; use them to reproduce the
# reported scores.
"$GCTB" --ldm-eigen "$GCTB_LDM_EIGEN" \
        --gwas-summary "$IMP/merged.ma" \
        --sbayes RC --annot "$GCTB_ANNOT" \
        --seed "$SBAYESRC_SEED" --thread "${SBAYESRC_THREADS:-4}" \
        --out "$RUN/SBayesR" \
        > "$RUN/gctb.log" 2>&1
echo "exit=$?  snpRes=$RUN/SBayesR.snpRes"
[ -s "$RUN/SBayesR.snpRes" ] || { echo "ABORT: no snpRes for $REGION"; exit 1; }
echo "snpRes md5=$(md5sum "$RUN/SBayesR.snpRes" | cut -d' ' -f1)"
