#!/usr/bin/env python
"""compute_bag.py - raw, bias-corrected and INT brain-age gaps from predictions.

Input: one or more prediction CSVs written by ``predict.py`` (or the
``{region}_oof_predictions.csv`` written by ``train_cv.py``). Every file needs
``chronological_age`` and ``predicted_age``; ``region`` is used to label the
output when present.

Bias correction (the paper's definition, optional via ``--no_bias_correction``):
fit ``predicted_age = a * age + b`` on a reference sample — the out-of-fold
predictions of healthy controls — then

    bias_corrected_bag = (predicted_age - b) / a - chronological_age

Estimating the coefficients and applying them are separate steps. The reference
sample the coefficients are ESTIMATED on must be named explicitly, one of:
  1. ``--calibration_csv`` (e.g. the CN out-of-fold prediction table), or
  2. the rows of the input selected by ``--calibration_filter COLUMN=VALUE``
     (a trailing ``*`` makes it a prefix match, e.g. ``prediction_mode=oof*``), or
  3. ``--calibration_scope self``, which fits on all rows of the input and is
     correct only when the input already *is* the reference table.
Once estimated, the coefficients are APPLIED unchanged to every row of every
predictions file given. Passing none of the three is an error: a missing
reference must never fall through to fitting the correction on the target sample
itself, which would shrink its own brain-age gap towards zero.

The UK Biobank and ADNI analyses use different references and must not be mixed:
the UK Biobank coefficients were estimated on the pooled out-of-fold predictions
of the 25,642 training participants and applied unchanged to training and to both
test sets, while ADNI has its own cognitively-normal out-of-fold reference.
(This guard is a release addition; the published runs named their reference
explicitly, so it changes no published number.)

INT (optional via ``--int``) is applied to the bias-corrected BAG within the
output file: ``norm.ppf((rankdata(x, 'average') - 0.5) / n)``.

Outputs the per-row table plus a JSON sidecar recording (a, b), the reference
sample size and its source, so every number is traceable.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from regionbae.postprocess import (  # noqa: E402
    apply_bias_correction, fit_bias_correction, inverse_normal_transformation, raw_bag,
)

AGE_COL = 'chronological_age'
PRED_COL = 'predicted_age'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--predictions_csv', nargs='+', required=True,
                   help='prediction CSV(s); processed independently, then concatenated')
    p.add_argument('--output_csv', required=True)
    p.add_argument('--calibration_csv', default='',
                   help='reference sample for the bias-correction fit '
                        '(control out-of-fold predictions)')
    p.add_argument('--calibration_filter', default='',
                   help='COLUMN=VALUE selecting the reference rows inside the input')
    p.add_argument('--calibration_scope', default='', choices=['', 'self'],
                   help="'self' fits the correction on all rows of the input; only correct "
                        'when the input already is the reference table')
    p.add_argument('--no_bias_correction', action='store_true',
                   help='report raw BAG only')
    p.add_argument('--int', dest='apply_int', action='store_true',
                   help='also report the inverse-normal transformed bias-corrected BAG')
    return p.parse_args()


def reference_rows(df, calibration_csv, calibration_filter, calibration_scope):
    """Return (reference_dataframe, provenance_string).

    Exactly one reference must be named. There is deliberately no default: falling
    back to the target sample would silently replace the paper's bias correction
    with a per-sample refit.
    """
    named = [bool(calibration_csv), bool(calibration_filter), calibration_scope == 'self']
    if sum(named) == 0:
        raise SystemExit(
            'ERROR: no bias-correction reference given. Pass --calibration_csv, or '
            '--calibration_filter COLUMN=VALUE, or --calibration_scope self if this '
            'predictions file already is the reference sample; or --no_bias_correction '
            'to report raw BAG only. This never defaults to fitting on the input.')
    if sum(named) > 1:
        raise SystemExit('ERROR: --calibration_csv, --calibration_filter and '
                         '--calibration_scope self are mutually exclusive')
    if calibration_csv:
        ref = pd.read_csv(calibration_csv)
        for col in (AGE_COL, PRED_COL):
            if col not in ref.columns:
                raise KeyError(f'{calibration_csv} is missing required column "{col}"')
        return ref, f'calibration_csv={calibration_csv}'
    if calibration_filter:
        if '=' not in calibration_filter:
            raise ValueError('--calibration_filter must look like COLUMN=VALUE')
        col, value = calibration_filter.split('=', 1)
        if col not in df.columns:
            raise KeyError(f'--calibration_filter column "{col}" not in the predictions table')
        if value.endswith('*'):  # prefix match, e.g. prediction_mode=oof*
            ref = df[df[col].astype(str).str.startswith(value[:-1])]
        else:
            ref = df[df[col].astype(str) == value]
        if ref.empty:
            raise ValueError(f'--calibration_filter {calibration_filter} selected no rows')
        return ref, f'filter={calibration_filter}'
    return df, 'all rows of the predictions table (--calibration_scope self)'


def main():
    args = parse_args()
    frames, provenance = [], []

    for path in args.predictions_csv:
        df = pd.read_csv(path)
        for col in (AGE_COL, PRED_COL):
            if col not in df.columns:
                raise KeyError(f'{path} is missing required column "{col}"')

        df['raw_bag'] = raw_bag(df[AGE_COL], df[PRED_COL])
        region = str(df['region'].iloc[0]) if 'region' in df.columns and len(df) else Path(path).stem
        record = {'predictions_csv': path, 'region': region, 'n_rows': int(len(df))}

        if args.no_bias_correction:
            print(f'[{region}] raw BAG only (n={len(df)}, MAE={df["raw_bag"].abs().mean():.4f})')
        else:
            ref, source = reference_rows(df, args.calibration_csv, args.calibration_filter,
                                         args.calibration_scope)
            if 'region' in ref.columns and ref['region'].nunique() == 1:
                ref_region = str(ref['region'].iloc[0])
                if ref_region != region:
                    print(f'[WARN] calibration region {ref_region!r} != prediction region {region!r}; '
                          'bias-correction coefficients are region-specific')
            a, b = fit_bias_correction(ref[AGE_COL], ref[PRED_COL])
            df['bias_corrected_bag'] = apply_bias_correction(df[AGE_COL], df[PRED_COL], a, b)
            record.update({'bias_correction': {'a': a, 'b': b, 'n_reference': int(len(ref)),
                                               'reference_source': source}})
            print(f'[{region}] bias correction a={a:.6f} b={b:.6f} '
                  f'(n_ref={len(ref)}, {source}); '
                  f'corrected BAG mean={df["bias_corrected_bag"].mean():.4f}')

            if args.apply_int:
                df['int_bias_corrected_bag'] = inverse_normal_transformation(df['bias_corrected_bag'])
                record['int'] = {'scope': 'rows of this predictions file', 'n': int(len(df))}

        frames.append(df)
        provenance.append(record)

    out = pd.concat(frames, ignore_index=True)
    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    sidecar = out_path.with_suffix('.calibration.json')
    sidecar.write_text(json.dumps(provenance, indent=2))

    print(f'\nWritten: {out_path} ({len(out)} rows)')
    print(f'Calibration provenance: {sidecar}')


if __name__ == '__main__':
    main()
