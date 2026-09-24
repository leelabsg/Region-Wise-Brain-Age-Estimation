#!/bin/bash
# SAIGE Step 2 (single-variant association) -- UK Biobank, imputed BGEN.
#
# One job per (region, chromosome). Resumable: a chromosome whose result exists
# and whose log reports "Analysis done" is skipped.
#
# UK Biobank BGEN files label autosomes 01..09 / 10..22; --chrom is zero-padded
# to match, while output files keep the unpadded number.
#
# The option list is the one the published run passed, written out in full rather
# than left to SAIGE's defaults, so a future SAIGE release changing a default
# cannot silently change the analysis.
#
# Usage: bash saige_step2.sh [chr ...]                  (default 1-22)
#        REGIONS="caudate thalamus" bash saige_step2.sh 22
set -uo pipefail
. "$(dirname "$0")/../config/common.sh"
gwas_require GWAS_UKB_DIR UKB_BGEN_TMPL UKB_BGEN_INDEX_TMPL UKB_SAMPLE_FILE

if [ -n "${REGIONS:-}" ]; then regions=($REGIONS); else mapfile -t regions < <(gwas_regions ukb); fi
if [ "$#" -ge 1 ]; then chrs=("$@"); else chrs=($(seq 1 22)); fi

launch() {
    local region=$1 chr=$2 chrom=$3 node=$4
    local res="$GWAS_UKB_DIR/$region/results" log="$GWAS_UKB_DIR/$region/logs"
    mkdir -p "$res" "$log"
    local cmd
    cmd="$(gwas_saige_cmd step2_SPAtests.R \
        --bgenFile="${UKB_BGEN_TMPL//\{chr\}/$chr}" \
        --bgenFileIndex="${UKB_BGEN_INDEX_TMPL//\{chr\}/$chr}" \
        --sampleFile="$UKB_SAMPLE_FILE" \
        --AlleleOrder=alt-first \
        --SAIGEOutputFile="$res/chr${chr}.txt" \
        --chrom="${chrom}" \
        --minMAF=0 \
        --minMAC=10 \
        --minInfo=0 \
        --maxMissing=0.15 \
        --impute_method=best_guess \
        --LOCO=TRUE \
        --GMMATmodelFile="$GWAS_UKB_DIR/$region/step1.rda" \
        --varianceRatioFile="$GWAS_UKB_DIR/$region/step1.varianceRatio.txt" \
        --SPAcutoff=2 \
        --dosage_zerod_cutoff=0.2 \
        --dosage_zerod_MAC_cutoff=10 \
        --markers_per_chunk=10000 \
        --max_MAC_for_ER=4 \
        --is_fastTest=FALSE \
        --is_Firth_beta=TRUE \
        --is_imputed_data=TRUE \
        --pCutoffforFirth=0.05 \
        --is_output_moreDetails=TRUE \
        --is_overwrite_output=TRUE)"
    gwas_run "$node" "$cmd" > "$log/step2_chr${chr}.log" 2>&1
}

i=0; launched=0; skipped=0
for region in "${regions[@]}"; do
    for chr in "${chrs[@]}"; do
        out="$GWAS_UKB_DIR/$region/results/chr${chr}.txt"
        if [ -s "$out" ] && grep -q "Analysis done" "$GWAS_UKB_DIR/$region/logs/step2_chr${chr}.log" 2>/dev/null; then
            skipped=$((skipped + 1)); continue
        fi
        if [ "$chr" -lt 10 ]; then chrom="0$chr"; else chrom="$chr"; fi
        gwas_throttle
        node="$(gwas_node $i)"
        echo "[$(date +%H:%M:%S)] launch $region chr$chr ${node:+-> $node}"
        launch "$region" "$chr" "$chrom" "$node" &
        i=$((i + 1)); launched=$((launched + 1))
    done
done
wait
echo "=== step2 done: launched=$launched skipped=$skipped ==="
