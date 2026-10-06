#!/usr/bin/env python
"""Stage 05 - Optuna tuning on the 70% development set (lgbm_functions.optuna_lgbm_cv, dart=True, no categorical).
Search space copied from the helper; boosting_type fixed to 'dart' (the helper lets Optuna choose gbdt/dart).
Objective: mean R2 over KFold(5, shuffle, seed 3456). TPESampler(seed 3456), maximize, 200 trials, sequential.
The 5 folds of one trial are fitted in parallel processes (independent fits; compute setting only).
No eval_set / early stopping: DART fits use all n_estimators, as in the helper (callbacks=None for dart).
Persistent SQLite study; completed trials are never re-run. Stale RUNNING trials from a crash are marked FAIL.
"""
import os, sys, json, time, datetime as dt
import numpy as np, pandas as pd, optuna, lightgbm as lgb
from sklearn.metrics import r2_score
from concurrent.futures import ProcessPoolExecutor

from _paths import WORK as RUN
SEED, NTRIALS = 3456, 200
FOLD_WORKERS = int(os.environ.get('FOLD_WORKERS', '5'))
THREADS_PER_FIT = int(os.environ.get('THREADS_PER_FIT', '16'))
BASE = dict(verbosity=-1, objective='regression', metric='None', n_estimators=5000, random_state=SEED)
STORAGE = f'sqlite:///{RUN}/optuna/protage204_dart.db'
STUDY = 'protage204_dart_tuning'

D = pd.read_parquet(f'{RUN}/preprocessing/clock_cohort.parquet').set_index('eid')
FEAT = json.load(open(f'{RUN}/preprocessing/feature_order.json'))
tf = pd.read_csv(f'{RUN}/splits/tuning_folds.tsv', sep='\t')
X = D.loc[tf.eid, FEAT]; y = D.loc[tf.eid, 'age_target']; fold = tf.tuning_fold.to_numpy()

def fit_fold(args):
    k, params = args
    tr, va = fold != k, fold == k
    m = lgb.LGBMRegressor(**params).fit(X.iloc[tr], y.iloc[tr])
    return k, float(r2_score(y.iloc[va], m.predict(X.iloc[va])))

def objective(trial):
    p = {
        'boosting_type': 'dart',
        'subsample': trial.suggest_float('subsample', 0.1, 1.0),
        'num_leaves': trial.suggest_int('num_leaves', 2, 256),
        'max_depth': trial.suggest_int('max_depth', 1, 4),
        'min_child_samples': trial.suggest_int('min_child_samples', 5, 300),
        'learning_rate': trial.suggest_float('learning_rate', 1e-5, 1, log=True),
        'min_child_weight': trial.suggest_float('min_child_weight', 1e-5, 100, log=True),
        'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1, log=True),
        'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1, log=True),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.1, 1.0),
    }
    p.update(BASE)
    p.update(drop_rate=trial.suggest_float('drop_rate', 0.1, 1.0), skip_drop=trial.suggest_float('skip_drop', 0.1, 1.0),
             max_drop=trial.suggest_int('max_drop', 1, 100), drop_seed=SEED, n_jobs=THREADS_PER_FIT)
    t0 = time.time()
    with ProcessPoolExecutor(FOLD_WORKERS) as ex:
        res = dict(ex.map(fit_fold, [(k, p) for k in range(5)]))
    scores = [res[k] for k in range(5)]
    trial.set_user_attr('fold_r2', scores); trial.set_user_attr('seconds', round(time.time() - t0, 1))
    print(f'[{dt.datetime.utcnow():%FT%TZ}] trial {trial.number}: mean R2 {np.mean(scores):.5f} ({time.time() - t0:.0f}s)', flush=True)
    return float(np.mean(scores))

if __name__ == '__main__':
    study = optuna.create_study(direction='maximize', storage=STORAGE, study_name=STUDY, load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(seed=SEED))
    for t in study.get_trials(deepcopy=False, states=(optuna.trial.TrialState.RUNNING,)):
        study._storage.set_trial_state_values(t._trial_id, state=optuna.trial.TrialState.FAIL)
        print(f'marked stale RUNNING trial {t.number} as FAIL', flush=True)
    done = len(study.get_trials(deepcopy=False, states=(optuna.trial.TrialState.COMPLETE,)))
    remaining = NTRIALS - done
    print(f'completed trials {done}; running {remaining} more (fold workers {FOLD_WORKERS} x threads {THREADS_PER_FIT})', flush=True)
    if remaining > 0:
        study.optimize(objective, n_trials=remaining, gc_after_trial=True, catch=(Exception,))
    st = study.get_trials(deepcopy=False)
    n_ok = sum(t.state == optuna.trial.TrialState.COMPLETE for t in st)
    n_fail = sum(t.state == optuna.trial.TrialState.FAIL for t in st)
    if n_ok < NTRIALS:
        print(f'incomplete: {n_ok} complete, {n_fail} failed; rerun to continue', flush=True); sys.exit(1)
    study.trials_dataframe().to_csv(f'{RUN}/optuna/trial_history.tsv', sep='\t', index=False)
    best = dict(study.best_params); best.update(boosting_type='dart', drop_seed=SEED)
    json.dump(dict(best_params=best, best_value_mean_cv_r2=study.best_value, best_trial=study.best_trial.number,
                   n_complete=n_ok, n_failed=n_fail, sampler='TPESampler(seed=3456)', base_params=BASE,
                   note='boosting_type fixed to dart; n_jobs is a compute setting'),
              open(f'{RUN}/optuna/selected_hyperparameters.json', 'w'), indent=1)
    open(f'{RUN}/logs/s04.done', 'w').write(dt.datetime.utcnow().isoformat() + 'Z\n')
    print('stage 04 done', flush=True)
