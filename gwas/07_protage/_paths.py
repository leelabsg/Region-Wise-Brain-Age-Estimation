"""Paths for the ProtAge-204 pipeline, read from the environment (see gwas/config/paths.env.example)."""
import os, sys


def env(name, default=None):
    v = os.environ.get(name, default)
    if not v:
        sys.exit(f'ERROR: {name} is not set (see gwas/config/paths.env.example)')
    return v


WORK = env('PROTAGE_WORK_DIR')
PRE, QC, SPL, LOGS = (os.path.join(WORK, d) for d in ('preprocessing', 'qc', 'splits', 'logs'))
for _d in (PRE, QC, SPL, LOGS, os.path.join(WORK, 'optuna'), os.path.join(WORK, 'models'), os.path.join(WORK, 'results')):
    os.makedirs(_d, exist_ok=True)
