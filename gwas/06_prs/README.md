# Polygenic scores for the regional brain-age gaps

SBayesRC weights from the UK Biobank training split, scored in the held-out UK
Biobank sample and in ADNI. Produces Supplementary Tables 17, 18 and 21-23, and
the deposited weight files.

GCTB and PLINK are third-party tools and are not included; the scripts here are
the preprocessing, the exact tool invocations, and the analyses.

## Reproducibility of the weights

**`--seed` is not enough. `--seed` with `--thread 1` is.**

SBayesRC is an MCMC sampler, and `--thread` is the variable that matters. Run
the whole-brain region from a byte-identical `.ma` input with the same seed:

| threads | runs | distinct outputs |
|---|---|---|
| `--thread 4` | 5 (across three nodes) | **3** |
| `--thread 1` | 2 | **1** — byte-identical |

Every deposited weight set was itself produced by a `--thread 4` run at
`--seed 20260821`, one run per region. For the whole-brain region both
single-threaded runs reproduce that deposited file byte for byte, so there the
deposit is the deterministic answer rather than one parallel run's draw.
Single-threaded cost is about twice: 1 h 51 m against 56 m for four threads, per
region.

Two things the parallel variation is **not**, because the obvious explanations
are wrong:

* not a hardware or library effect. Every node involved has the same CPU model,
  no AVX-512, no threading environment variables set, and runs the same
  statically linked GCTB binary;
* not across-node versus within-node. One weight set appeared on three different
  nodes, and a single node produced two different weight sets. The variation is
  between runs, and it disappears at one thread.

So: **run `--thread 1` if you need to reproduce the weights**, and use more
threads only when the downstream aggregate results are what you care about. How
far the parallel variation propagates when it does occur, measured on the
whole-brain region:

| level | quantity | value |
|---|---|---|
| weights | corr(BETA) between distinct chains | 0.639 – 0.739 |
| weights | Jaccard of the top 1,000 by \|BETA\| | 0.20 – 0.23 |
| hyperparameters | SNP h² across chains | 0.173 – 0.187 (within 0.4 posterior SD) |
| PRS, held-out UKB (n = 10,267) | corr(PRS) between chains | **0.950 – 0.959** |
| PRS | max abs. difference after standardising | 1.15 – 1.37 SD |
| PRS | individuals keeping their decile | 45 – 48 % |
| ST17 | β | +0.886 to +0.901 |
| ST17 | P | 3.8e-60 to 4.9e-62 |
| ST17 | adjusted R² | 0.0256 – 0.0265 |
| ST17 | associations changing direction or significance | **0 of 3 chains** |

So the reported associations are insensitive to which chain is used, while an
individual's score is not: about half of the sample moves decile between chains,
and roughly a third of the top 10 % is replaced. Per-participant scores should
be taken from the deposited weights rather than regenerated.

This was measured for the whole-brain region. The remaining regional weights are
the archived canonical outputs of one `--thread 4` run each; their single-thread
bitwise reproducibility was not independently tested, and neither was the spread
above.

`07_export_weights.py` remains the practical answer: deposit the weights and
nothing downstream depends on re-running the sampler at all. Scoring, ST17,
ST18 and ST21-23 are deterministic given the weights. It is create-only: it
refuses to overwrite an existing weight file or an existing `MANIFEST.txt`, so a
second run against a directory that already holds an export stops with exit 1
rather than replacing a package that may already have been deposited.

## Order

```
  regional BAG GWAS (01_ukb_gwas, training split)
          │
          ▼
  01_preprocess_gwas.py      SAIGE sumstats → GCTB .ma, with the allele check
          │
          ▼
  02_gctb_commands.sh        impute to the LD reference → SBayesRC   [per region]
          │
          ├─────────────► 07_export_weights.py   deposit package + sha256 manifest
          │
          ▼
  03_score.sh ukb            plink2 --score, chr1-22
          │
          ├──► 04_analyze_ST17.py          held-out validation   (ST17)
          └──► 05_analyze_ST18_phewas.py   PheWAS, once per sample (ST18)

  04a_build_adni_weights.py  → 03_score.sh adni → 06_analyze_adni.py  (ST21-23)
```

```bash
set -a; . gwas/config/paths.env; set +a
for r in $(awk -F'\t' '!/^#/ && $3 ~ /ukb/ {print $1}' gwas/config/regions.tsv); do
    bash gwas/06_prs/02_gctb_commands.sh "$r"
    bash gwas/06_prs/03_score.sh ukb "$r"      # PRS_UKB_KEEP unset scores everyone
done
python gwas/06_prs/07_export_weights.py
python gwas/06_prs/04_analyze_ST17.py
```

## The allele mapping

GCTB's `.ma` wants `A1` to be the effect allele and `freq` to be the frequency of
`A1`. SAIGE reports `BETA` and `AF_Allele2` against `Allele2`, so `Allele2` is
the effect allele and the mapping is `A1 <- Allele2`, `A2 <- Allele1`.

An earlier version of this pipeline set `A1 <- Allele1` while leaving `b` and
`freq` on `Allele2`. GCTB accepts that without complaint, the run finishes, and
every downstream sign is inverted. `01_preprocess_gwas.py --check-orientation`
re-reads the SAIGE file and asserts the mapping held; `02_gctb_commands.sh`
passes it on every chromosome. `06_analyze_adni.py` reports the independent
check on the other side: correctly oriented weights correlate positively with
the measured brain-age gap in ADNI.

## Not included

The cluster orchestration — node assignment, ssh fan-out, job counting, progress
watchers — is not here. It is a property of one cluster, not of the method. The
loops in `02_gctb_commands.sh` and `03_score.sh` are serial and state the
parameters plainly; parallelise them however suits the machine, and the results
do not change.

## Files

| file | what it does |
|---|---|
| `01_preprocess_gwas.py` | SAIGE summary statistics → GCTB `.ma`, with `--check-orientation` |
| `02_gctb_commands.sh` | `--impute-summary` per LD block, merge, `--sbayes RC` |
| `03_score.sh` | `plink2 --score` for UK Biobank (per chromosome) or ADNI |
| `04_analyze_ST17.py` | held-out UKB validation, OLS BAG ~ PRS (ST17) |
| `04a_build_adni_weights.py` | rewrite the weights against ADNI variant IDs |
| `05_analyze_ST18_phewas.py` | UK Biobank PheWAS (ST18), one run per analysis sample; sex-specific phecodes are restricted to that sex in cases and controls alike, a phecode's exclusion range is dropped from that phecode's test rather than counted as a control, and BH runs within one BAG PRS within one sample |
| `sex_specific_phecodes.tsv` | the sex designations ST18 applied: 134 aggregated test units, 109 Female and 25 Male |
| `06_analyze_adni.py` | ADNI correlation and case-control models (ST21-23) |
| `07_export_weights.py` | deposit package: per-region weights + sha256 manifest |
| `_regions.py` | region codes and display names, from `config/regions.tsv` |

## ST18: three analysis samples

`05_analyze_ST18_phewas.py` is run once per sample, each with its own
`ST18_SAMPLE` and its own output directory; the samples are never pooled and BH
is computed within one BAG PRS within one sample.

| `ST18_SAMPLE` | analysis sample | n | role |
|---|---|---|---|
| `ukb_nonimaging_wb` | UK Biobank non-imaging, White British | 389,045 | primary |
| `ukb_nonimaging` | UK Biobank non-imaging, no ancestry restriction | 442,295 | sensitivity |
| `ukb_imaging_heldout_10267` | UK Biobank imaging held-out | 10,267 | |

The ancestry restriction is supplied as `PRS_PHEWAS_KEEP`, a list of participant
ids. Omitting it is a decision, not a default: a run with no keep file must set
`ST18_ALLOW_UNRESTRICTED=1` to state that `PRS_PHEWAS_TABLE` is already the
intended sample, and a keep file that is missing, empty or matches nothing is an
error rather than a fallback to the whole table. `ST18_EXPECTED_N` asserts the
resulting sample size.

Three model forms are fitted:

| model | phecodes | design |
|---|---|---|
| A | not sex-specific | `logit(case) ~ PRS_region + Age + Sex + PC1..PC10` |
| B | sex-specific | `logit(case) ~ PRS_region + Age + PC1..PC10`, fitted in that sex only, cases and controls alike; Sex is dropped because it is constant after the restriction |
| C | 751.1 | model A, unrestricted: its constituent subcodes carry conflicting sex designations |

Test units are the integer part of a phecode plus at most one decimal digit, so
174.11 folds into 174.1 while 174.1 and 174 stay separate units; within a unit the
child columns are combined with a skipna maximum. A test is fitted when it has at
least 30 cases and both outcome classes after the restriction and complete-case
selection; every attempt, fitted or excluded, is written to `ST18_fit_status.tsv`
with its reason.

`PRS_UKB_KEEP` in `03_score.sh` is a different thing: it selects who gets
*scored*, not who is analysed. The published run left it unset and scored every
participant in the pgen, because ST18 needs scores outside the imaging cohort.
