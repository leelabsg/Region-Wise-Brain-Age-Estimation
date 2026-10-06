#!/bin/bash
# Munge the two external ageing traits that ST15 correlates against the regional
# BAG GWAS: telomere length PC1 and the ProtAge-204 proteomic-age gap.
#
# The munge command is NOT repeated here. This script only resolves the
# trait-specific sample size and hands off to 02_munge.sh, which is the same
# invocation used for pad_<region>, so there is exactly one place in the release
# where munge_sumstats.py is called and one place where --a1/--a2 are set.
#
# What differs per trait is the input, the preprocess step that produced the
# per-chromosome files under $LDSC_WORK_DIR/files/<trait>, and --N:
#
#   telomere_pc1   N = 438351   00_preprocess_telomere.py
#                  Input: LDSC_TELOMERE_RAW_SUMSTATS (GWAS Catalog GCST90435144).
#                  N is the file's own `n` column. It is constant at 438351 over
#                  all 9,589,468 rows of the 22 per-chromosome preprocess outputs
#                  that were munged. The raw Catalog download carries a second
#                  value for a minority of variants (438351 for 14,587,804 rows;
#                  438170 for 434,898), all of which fall outside the munged set.
#                  The released sumstats were munged with --N 438351.0.
#
#   protage204     N = 38806    01_preprocess.py
#                  Input: LDSC_PROTEIN_AGE_RAW_SUMSTATS, the merged SAIGE step-2
#                  output for the ProtAge-204 proteomic-age gap. N is the SAIGE
#                  analysed sample size and equals the N column (field 14) of
#                  that file, which is checked to be constant before munging.
#                  The published run (07_protage/11_ldsc_downstream.py, which
#                  makes the same munge call) used --N 38806.
#
# ProtAge-204 needs no external-specific preprocess script. Its input is ordinary
# SAIGE output (CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 ... BETA
# SE ... p.value N), so 01_preprocess.py applies to it unchanged.
#
# Usage: bash 00_munge_external.sh telomere_pc1
#        bash 00_munge_external.sh protage204
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
. "$HERE/../config/common.sh"
gwas_require LDSC_WORK_DIR

TRAIT="${1:?Usage: $0 <telomere_pc1|protage204>}"
case "$TRAIT" in
    telomere_pc1)
        gwas_require LDSC_TELOMERE_N
        N="$LDSC_TELOMERE_N"
        ;;
    protage204)
        gwas_require LDSC_PROTEIN_AGE_N
        N="$LDSC_PROTEIN_AGE_N"
        ;;
    *)
        echo "ERROR: unknown external trait '$TRAIT' (expected telomere_pc1 or protage204)" >&2
        exit 1
        ;;
esac

# The per-chromosome inputs must already exist; run the matching preprocess step
# first (00_preprocess_telomere.py for telomere_pc1, 01_preprocess.py otherwise).
IN_DIR="$LDSC_WORK_DIR/files/$TRAIT"
for chrom in $(seq 1 22); do
    for path in \
        "${IN_DIR}/chr${chrom}_${TRAIT}_GWAS_SummaryStatistics.txt" \
        "${IN_DIR}/chr${chrom}_${TRAIT}_GWAS_SummaryStatistics_snplist.txt"
    do
        if [ ! -s "$path" ]; then
            echo "ERROR: missing preprocess output: $path" >&2
            exit 1
        fi
    done
done

echo "[00_munge_external] $TRAIT --N $N -> $LDSC_WORK_DIR/munge_output/$TRAIT"
bash "$HERE/02_munge.sh" "$TRAIT" "$N"
