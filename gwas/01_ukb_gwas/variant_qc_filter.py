#!/usr/bin/env python
"""Post-GWAS variant filter for the regional BAG summary statistics (the deposited files).

A ROW FILTER on the merged SAIGE output; BETA, SE and P are copied unchanged. A variant is retained when
  (a) it matches exactly one UK Biobank imputation mfi v3 record on GRCh37 chromosome, position and the
      unordered allele pair, and that record has INFO >= 0.8 (INFO 'NA' fails), and
  (b) the analysis-sample MAF = min(AF_Allele2, 1 - AF_Allele2) from the SAIGE output is >= 0.01.
mfi MAF is NOT a criterion; it is appended as the reference column mfi_MAF together with mfi_INFO. Unmatched and
ambiguous (several mfi records) variants are excluded. In the published run every region retained 9,458,475 of
35,784,891 tested variants. CHR is written without zero padding.

    python gwas/01_ukb_gwas/variant_qc_filter.py <merged_sumstats.txt> <mfi_dir> <out.tsv.gz>

<mfi_dir> holds ukb_mfi_chr{1..22}_v3.txt (tab-separated: id, rsid, position, allele1, allele2, MAF, minor allele, INFO).
"""
import gzip, json, sys

HEADER = ['CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2', 'AC_Allele2', 'AF_Allele2', 'imputationInfo', 'BETA', 'SE',
          'Tstat', 'var', 'p.value', 'N']
INFO_MIN, MAF_MIN = 0.8, 0.01


def nan_or_float(v):
    return float('nan') if v == 'NA' else float(v)


def load_mfi(path):
    d = {}
    with open(path) as f:
        for raw in f:
            c = raw.rstrip('\n').split('\t')
            if len(c) < 8 or not c[2].isdigit():
                continue
            d.setdefault((int(c[2]), tuple(sorted((c[3], c[4])))), []).append((nan_or_float(c[5]), nan_or_float(c[7])))
    return d


def main(src, mfi_dir, out):
    counts = dict(input_rows=0, unmatched=0, ambiguous=0, fail_info_only=0, fail_maf_only=0, fail_both=0, retained=0)
    with open(src) as f, gzip.open(out, 'xt', compresslevel=6) as o:
        if f.readline().rstrip('\n').split('\t') != HEADER:
            sys.exit('FAIL: unexpected header')
        o.write('\t'.join(HEADER + ['mfi_MAF', 'mfi_INFO']) + '\n')
        ch, d = None, None
        for line in f:
            c = line.rstrip('\n').split('\t')
            cc = int(c[0])
            if cc != ch:
                if ch is not None and cc < ch:
                    sys.exit('FAIL: input not in chromosome order')
                ch, d = cc, load_mfi(f'{mfi_dir}/ukb_mfi_chr{cc}_v3.txt')
            counts['input_rows'] += 1
            af = float(c[6]); maf = min(af, 1 - af)
            hits = d.get((int(c[1]), tuple(sorted((c[3], c[4])))), [])
            if len(hits) != 1:
                counts['unmatched' if not hits else 'ambiguous'] += 1
                continue
            mfi_maf, info = hits[0]
            ok_i = info == info and info >= INFO_MIN
            ok_m = maf >= MAF_MIN
            if ok_i and ok_m:
                counts['retained'] += 1
                o.write('\t'.join([str(cc)] + c[1:]) + f'\t{mfi_maf!r}\t{info!r}\n')
            elif ok_m:
                counts['fail_info_only'] += 1
            elif ok_i:
                counts['fail_maf_only'] += 1
            else:
                counts['fail_both'] += 1
    print(json.dumps(counts, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:4])
