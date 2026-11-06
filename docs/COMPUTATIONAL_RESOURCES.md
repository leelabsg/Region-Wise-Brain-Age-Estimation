# Computational resources

This page gives the hardware used for the analyses in the paper and the compute
each stage required. Times are wall-clock times and include data loading.
GPU-hours and CPU core-hours are summed over all runs of a stage.

## Hardware

| Component | Specification |
|---|---|
| GPU | NVIDIA RTX A5000, 24 GB; one GPU per training or inference run |
| CPU | 2 × AMD EPYC 7542 (32 cores each; 64 cores / 128 threads per node) |
| Memory | 256 GB per node |

The genetic analyses ran on CPU nodes with the same specification and did not
use a GPU.

## Software

| Stage | Software |
|---|---|
| GWAS, null model (step 1) | SAIGE 1.1.8 |
| GWAS, association testing (step 2) | SAIGE 1.1.9 |
| Polygenic scores | SBayesRC (GCTB 2.5.2) |

## Brain-age models

Each region has its own 3D CNN with 2,984,801 trainable parameters, trained with 4-fold
cross-validation. That makes ten regions × four folds = 40 models per cohort.
Inputs are 128 × 128 × 128 volumes.

| Stage | Data | Time per fold | Time per region | Total |
|---|---|---|---|---|
| UK Biobank training (random initialization) | 25,642 participants; ≈19,200 training / ≈6,400 validation per fold | 7.4–10.6 h (median 7.7 h) | 29.7–33.1 h | **311 GPU-hours** |
| ADNI fine-tuning (from the UK Biobank model of the same region and fold) | 1,179 cognitively normal participants; ≈884 training / ≈295 validation per fold | 44.8–56.6 min (median 48.8 min) | 3.1–3.4 h | **32.7 GPU-hours** |
| ADNI inference (four fold models per scan) | 10,823 scans | — | 13–21 min | **2.7 GPU-hours** |

Regions are independent and can run in parallel on separate GPUs. With one GPU
per region and the four folds run in sequence, all 40 UK Biobank models take
about 1.4 days, and fine-tuning all ten regions takes about 3.5 hours. Training,
fine-tuning and inference each fit on one 24 GB GPU.

## Genome-wide association studies

There are ten regional brain-age-gap phenotypes, each analysed in two UK
Biobank samples: discovery (N = 41,067) and polygenic-score derivation
(N = 30,800). Association testing runs as one job per phenotype per autosome,
using one CPU thread per job.

| Stage | Jobs | Threads per job | Time per job | CPU core-hours |
|---|---|---|---|---|
| SAIGE step 1, null model (215,514 LD-pruned markers) | 10 per sample | 8 | 10–50 min | 56 (discovery), 34 (derivation) |
| SAIGE step 2, discovery N = 41,067 | 220 (10 phenotypes × 22 chromosomes) | 1 | 1.5–10.8 h (median 5.2 h) | **1,147** |
| SAIGE step 2, derivation N = 30,800 | 220 | 1 | 1.4–10.4 h (median 5.0 h) | **1,092** |

Step 2 dominates the cost. Each job's time scales with the number of imputed
variants on its chromosome. The jobs are independent of one another, so the
wall-clock time is set by how many cores are available.

## Polygenic-score derivation

SBayesRC was run with the UK Biobank European HapMap3 LD reference and
4 threads per region. Each region took 59–68 min, for a total of 41.5 CPU
core-hours across the ten regions.

## Other analyses

LD score regression, MAGMA, PheWAS and Mendelian randomization each run on a
single CPU node and need no GPU.

## Summary

| Resource | Total |
|---|---|
| GPU, UK Biobank training | 311 GPU-hours |
| GPU, ADNI fine-tuning and inference | 35 GPU-hours |
| CPU, GWAS (two samples × ten phenotypes) | ≈2,330 core-hours |
| CPU, SBayesRC | 42 core-hours |

To fine-tune the released UK Biobank models on a new cohort of about 1,200
participants, plan on about 3.3 GPU-hours per region. Inference with all four
fold models takes 13–21 minutes per region for about 10,800 scans on one GPU.
