# Region-Wise Brain Age Estimation

Region-wise brain age estimation from 3D T1-weighted MRI. Regions are defined by
an atlas label map and masked **on the fly** at training and inference time, so a
single set of whole-brain volumes serves every region. One independently
parameterized 3D CNN is trained per region with k-fold cross-validation, and
brain-age gaps (BAG) are reported raw, bias-corrected, and optionally
inverse-normal transformed.

This repository is the cleaned-up release of the research code. Model
architecture, masking, optimization, cross-validation and the BAG formulas are
unchanged; only hard-coded paths and the fixed region table were generalized.

## Pipeline overview

This repository implements the core model-development workflow shown below: atlas-based anatomical masking of three-dimensional T1-weighted MRI, separate training of independently parameterized whole-brain and regional 3D-CNN models, cross-validation, brain-age prediction, and BAG post-processing.

<p align="center">
  <img width="4320" height="1440" alt="fig1_a_bag_derivation" src="https://github.com/user-attachments/assets/48cd559e-36e3-44cb-93ad-d5145fbec393" />
</p>

*Overview of the model-development workflow implemented in this repository. Each whole-brain or region-masked 3D MRI input is used to train a separate 3D-CNN brain-age model.*

## What the pipeline does

1. **Region definition** — a JSON config maps region names to atlas label values
   (`configs/regions_mni_structural.json`). Any atlas and any label combination
   can be used; a region may take a single label or the union of several.
2. **On-the-fly masking** — the label mask is binary-dilated on the atlas grid,
   nearest-neighbour resampled onto the image grid, and applied as
   `image[~mask] = 0` inside the data loader. No region images are written to
   disk. `"global": null` means whole brain, no masking.
3. **Per-region training** — one run trains one region: a fresh CNN, its own
   checkpoints, its own results directory. No weights are shared between regions.
4. **Cross-validation** — k-fold (default 4) with per-fold checkpoints; the epoch
   with the lowest validation MAE becomes `cv-{fold}-best.pth.tar`.
5. **Prediction** — every fold model scores every row; the final value is either
   the subject's own out-of-fold model or the fold ensemble.
6. **Brain-age gaps** — raw, bias-corrected (optional), INT (optional), with the
   calibration coefficients written to a JSON sidecar.

## Installation

```bash
pip install -r requirements.txt
```

ANTsPy (`antspyx`) is required for atlas mask construction; PyTorch with CUDA is
strongly recommended for training.

## Input data expectations

The data loader deliberately performs **no** preprocessing: it reads a NIfTI
volume, applies the region mask, and feeds raw float voxels to the network.
Cohort images must therefore already be preprocessed and mutually aligned before
training. In the paper, this meant:

- skull-stripped T1w, N4 bias-field corrected,
- non-linearly registered to an MNI152 **brain** template (brain-to-brain, so the
  moving and fixed images have the same field of view),
- resampled to a common isotropic grid (128×128×128),
- ADNI: intensity-normalized per volume after the 128×128×128 resampling and
  before region masking (intensities scaled so that the 99th percentile of
  positive voxels equals 105; non-positive voxels set to 0; no upper clipping).

The atlas does not need to be on the same grid as the images — it is resampled to
the image grid with nearest-neighbour interpolation — but it **must be in the
same space** (same template, same orientation). Mismatched shapes raise an error;
a mismatched *space* would silently mask the wrong tissue, so verify the overlay
once before a full run.

Cohort tables are CSV files with at least:

| column | meaning |
|---|---|
| `imgs` | absolute path to the preprocessed whole-brain volume |
| `age` | chronological age at scan (years) |
| `subjectID` | optional; carried into the prediction outputs |

See `examples/example_cohort.csv`. Optional columns used by `predict.py`:
`prediction_mode` (`oof_fold_K` or `ensemble`).

## Usage

### Train one region (4-fold CV)

```bash
python scripts/train_cv.py \
    --cohort_csv  /path/to/cohort_train.csv \
    --region      temporal_lobe \
    --region_config configs/regions_mni_structural.json \
    --atlas       /path/to/MNI-maxprob-thr0-1mm.nii.gz \
    --model_dir   /path/to/output/models/temporal_lobe \
    --results_dir /path/to/output/results/temporal_lobe \
    --n_splits 4 --seed 7 --epochs 40 --batch_size 8 --lr 1e-4
```

`examples/run_training.sh` loops this over every region in the config. Regions are
independent, so they can also be run in parallel on separate GPUs.

Outputs per fold `F`:

```
{model_dir}/cv-{F}-{epoch}.pth.tar       every epoch
{model_dir}/cv-{F}-best.pth.tar          lowest validation MAE
{results_dir}/{F}/best/{region}.pkl      predictions of the best epoch
{results_dir}/{region}_oof_predictions.csv   concatenated out-of-fold predictions
```

### Fine-tune pretrained regional models on another cohort

```bash
python scripts/finetune_fold.py \
    --train_csv /path/to/folds/cv_fold_0_train.csv \
    --val_csv   /path/to/folds/cv_fold_0_val.csv \
    --region    temporal_lobe \
    --region_config configs/regions_mni_structural.json \
    --fold 0 \
    --model_load_folder /path/to/pretrained/temporal_lobe \
    --load_epoch 40 \
    --model_dir   /path/to/output/models_finetuned/temporal_lobe \
    --results_dir /path/to/output/results_finetuned/temporal_lobe/fold_0
```

Fold `F` of the target cohort starts from fold `F` of the pretrained model of the
**same region**. Every layer is fine-tuned — nothing is frozen — and the optimizer
is rebuilt after loading so no pretraining optimizer state carries over.
`--from_scratch` runs the same recipe from random initialisation.

`train_cv.py --mode finetune` does the same thing when the folds come from the
internal KFold instead of external fold CSVs.

### Predict and compute brain-age gaps

```bash
python scripts/predict.py \
    --cohort_csv /path/to/cohort_all.csv \
    --region temporal_lobe \
    --region_config configs/regions_mni_structural.json \
    --model_dir /path/to/output/models_finetuned/temporal_lobe \
    --n_folds 4 \
    --output_csv /path/to/output/predictions/temporal_lobe.csv

python scripts/compute_bag.py \
    --predictions_csv /path/to/output/predictions/temporal_lobe.csv \
    --calibration_filter "prediction_mode=oof*" \
    --int \
    --output_csv /path/to/output/bag/temporal_lobe_bag.csv
```

## Brain-age gap definitions

Let `y` be chronological age and `ŷ` the predicted brain age.

**Raw**

```
raw_bag = ŷ − y
```

**Bias-corrected** (de Lange & Cole style, division form). Fit
`ŷ = a·y + b` by ordinary least squares on a reference sample — in the paper, the
**out-of-fold predictions of cognitively normal controls** — then

```
bias_corrected_bag = (ŷ − b) / a − y
```

Estimating the coefficients and applying them are separate steps. Within each
cohort-specific calibration scheme, `(a, b)` are estimated once on the designated
reference sample and then applied unchanged to the corresponding target predictions.
For UK Biobank, the reference comprises pooled out-of-fold predictions from the
25,642 training participants, and the resulting coefficients are applied unchanged
to all UK Biobank prediction sets. ADNI uses a separate reference comprising
out-of-fold predictions from cognitively normal participants; UK Biobank and ADNI
calibration coefficients are not mixed. Refitting the correction within a target
sample or patient subgroup would therefore differ from the procedure used in the
paper.

`compute_bag.py` requires the reference to be specified explicitly using
`--calibration_csv`, `--calibration_filter COLUMN=VALUE`, or
`--calibration_scope self` when the input itself is the reference table. If none
is supplied, the script exits rather than fitting the correction on the target
input by default. The fitted `(a, b)`, reference-sample size and reference source
are written to a `.calibration.json` sidecar.

MAE and the Pearson correlation between predicted and chronological age are
reported on the **uncorrected** predicted age.

**INT** (optional), applied to the bias-corrected BAG:

```
int_bias_corrected_bag = Φ⁻¹((rank − 0.5) / n),  rank = rankdata(x, method='average')
```

Ranks are computed over the finite values of the file being processed; non-finite
entries stay `NaN`. Order is always raw → bias correction → INT.

## Repository layout

```
regionbae/
  regions.py        atlas + region-label config (the only generalized component)
  dataset.py        on-the-fly region masking data loader
  CNN.py            3D CNN (unchanged)
  CNN_Trainer.py    training / validation / checkpointing (unchanged; W&B optional)
  lr_scheduler.py   cosine annealing with warm-up restarts (unchanged)
  utils.py          seeding, scheduler factory, logger shim
  postprocess.py    raw / bias-corrected / INT brain-age gaps
scripts/
  train_cv.py       k-fold training of one region
  finetune_fold.py  fine-tune one pretrained region model on one fold
  predict.py        fold-wise / ensemble inference + raw BAG
  compute_bag.py    bias correction + INT with calibration provenance
configs/            example region configuration
examples/           example cohort table and end-to-end shell scripts
gwas/               genetic analysis of the brain-age gaps (see gwas/README.md)
```

## Genetic analysis

Everything downstream of the model — the GWAS of each regional brain-age gap in
UK Biobank, its replication in ADNI, LD score regression, MAGMA cell-type
enrichment, the BAG polygenic scores and their PheWAS, and the ProtAge-204
proteomic-age gap used as an external ageing measure — is in [`gwas/`](gwas/), with its own
[README](gwas/README.md) covering the order of execution, what data access each
stage needs, the FUMA dependency, and how to verify that the released scripts
reproduce the published tables.

## Data availability

No imaging data, participant identifiers, phenotypes, or trained model weights are
included in this repository.

GWAS summary statistics (ten BAGs, N = 41,067; ProtAge-204, N = 38,806) and the
BAG polygenic-score weights (derivation n = 30,800) are available from
https://storage.googleapis.com/leelabsg-regional-brain-age/README.md, which
lists the files, their columns, filters and checksums. UK Biobank and ADNI data are available to approved
researchers through the respective data-access processes. The atlas file is not
redistributed here — supply your own (e.g. the FSL MNI structural atlas) and note
its own licence terms.

## Related prior work

This repository builds on methodological work developed in the following prior studies:

- Kim J, Lee J, Lee S. Investigation of Genetic Variants and Causal Biomarkers Associated with Brain Aging. *medRxiv*. 2022.03.04.22271813. https://doi.org/10.1101/2022.03.04.22271813

- *A Multitask Deep Learning Model for Voxel-level Brain Age Estimation.* MLMI Workshop at MICCAI 2023. https://arxiv.org/abs/2310.11385

## License

MIT — see [LICENSE](LICENSE).
