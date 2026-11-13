"""Frozen measured-fluorescence transfer screen; no latent-rate truth is assumed."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

SOURCE_SHA = '5f9e49088e2218d6f35bca61f3039da72864906de711740437949e11cb47ec0b'
MODELS = ('single', 'plateau', 'mixture', 'raw_contrast_mixture', 'empirical_mean')


def survival(t, parameters, model):
    if model == 'single':
        return np.exp(-parameters[0]*t)
    if model == 'plateau':
        k, plateau = parameters
        return plateau+(1-plateau)*np.exp(-k*t)
    slow, gap, weight = parameters
    return weight*np.exp(-slow*t)+(1-weight)*np.exp(-(slow+gap)*t)


def load_movies(root):
    path = root/'source/combined_dataframe.csv'
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_SHA:
        raise ValueError('Source checksum mismatch')
    source = pd.read_csv(path)
    metadata = pd.read_csv(root/'extrusion/movies.csv', dtype={'date_raw': str})
    curves = {}
    for row in metadata.itertuples(index=False):
        group = source[source.filename == row.filename]
        means = group.groupby('frame')[['mean_list_interp_unbleached', 'mean_list_interp_bleached']].mean()
        difference = (means.iloc[:, 0]-means.iloc[:, 1]).to_numpy()
        control = group.drop_duplicates('frame').sort_values('frame')
        reference = (control.unfrap_cell-control.background_int).to_numpy()
        if difference[5] == 0 or np.any(reference <= 0):
            raise ValueError('Invalid normalization; retain and resolve before evaluation')
        x, z = difference[5:]/difference[5], reference[5:]/reference[5]
        curves[row.movie_id] = dict(x=x, z=z, y=x/z)
    return metadata, curves


def fit(t, x, z, model):
    y = x/z
    if model == 'empirical_mean':
        return dict(prediction=y.mean(axis=0).tolist(), parameters=None, starts=[], success=True)
    if model == 'single':
        starts = [[k] for k in [.001, .01, .1]]
        bounds = ([1e-6], [5.])
    elif model == 'plateau':
        starts = [[k, .25] for k in [.001, .01, .1]]
        bounds = ([1e-6, 0.], [5., 1.])
    else:
        starts = [[k, gap, w] for k in [.001, .01, .1] for gap in [.01, .1] for w in [.25, .75]]
        bounds = ([1e-6, 1e-6, 0.], [5., 5., 1.])
    def residual(parameters):
        predicted = survival(t, parameters, model)
        return ((x-z*predicted) if model == 'raw_contrast_mixture' else (y-predicted)).ravel()
    records = []
    for start in starts:
        result = least_squares(residual, start, bounds=bounds, max_nfev=1500,
                               ftol=1e-10, xtol=1e-10, gtol=1e-10)
        records.append(dict(parameters=result.x.tolist(), objective=float(np.mean(result.fun**2)),
                            success=bool(result.success), status=int(result.status), evaluations=int(result.nfev),
                            active_bounds=result.active_mask.tolist()))
    valid = [r for r in records if r['success']]
    if not valid:
        return dict(parameters=None, prediction=None, starts=records, success=False)
    best = min(valid, key=lambda r:r['objective'])
    return dict(parameters=best['parameters'], prediction=survival(t, best['parameters'], model).tolist(),
                starts=records, success=True, objective=best['objective'])


def run(root, output):
    metadata, curves = load_movies(root)
    t = np.arange(245)*.5
    fits, sensitivity, rows, predictions = {}, [], [], []
    conditions = sorted(metadata.condition.unique())
    maximum_identity_error = 0.
    for condition in conditions:
        development = metadata[(metadata.condition == condition)&(metadata.split == 'development')]
        ids = development.movie_id.tolist()
        x = np.stack([curves[m]['x'] for m in ids]); z = np.stack([curves[m]['z'] for m in ids])
        fits[condition] = {}
        for model in MODELS:
            fitted = fit(t, x, z, model)
            fits[condition][model] = fitted
            if not fitted['success']:
                raise RuntimeError(f'All starts failed: {condition}, {model}; no test evaluation permitted')
            s = np.array(fitted['prediction'])
            maximum_identity_error = max(maximum_identity_error, float(np.max(np.abs((x-z*s)**2-z*z*(x/z-s)**2))))
        for drift in [-.05, .05]:
            altered = z*(1+drift*t/t[-1])
            for model in ['mixture', 'raw_contrast_mixture']:
                alternative = fit(t, x, altered, model)
                baseline = np.array(fits[condition][model]['prediction'])
                sensitivity.append(dict(condition=condition, model=model, imposed_end_drift=drift,
                                        **alternative,
                                        max_curve_change=None if not alternative['success'] else float(np.max(np.abs(np.array(alternative['prediction'])-baseline)))))
    # All development choices above are fixed before these records are evaluated.
    for movie in metadata[metadata.split == 'reserved_test'].itertuples(index=False):
        target = curves[movie.movie_id]['y']
        for model in MODELS:
            predicted = np.array(fits[movie.condition][model]['prediction'])
            residual = target-predicted
            rows.append(dict(movie_id=movie.movie_id, condition=movie.condition, date_raw=movie.date_raw,
                             model=model, mse=float(np.mean(residual**2)), mae=float(np.mean(abs(residual))),
                             max_error=float(np.max(abs(residual)))))
            for frame, (observed, prediction) in enumerate(zip(target, predicted), start=5):
                predictions.append(dict(movie_id=movie.movie_id, condition=movie.condition, model=model,
                                        frame=frame, observed=float(observed), prediction=float(prediction)))
    table = pd.DataFrame(rows)
    balanced = table.groupby(['condition', 'model']).mse.mean().unstack('model')
    comparisons = {}
    for baseline in ['mixture', 'empirical_mean']:
        improvement = 1-balanced.raw_contrast_mixture.mean()/balanced[baseline].mean()
        wins = int((balanced.raw_contrast_mixture < balanced[baseline]).sum())
        comparisons[baseline] = dict(relative_mse_reduction=float(improvement), conditions_improved=wins,
                                     pass_gate=bool(improvement >= .10 and wins >= 4))
    result = dict(source_sha256=SOURCE_SHA, protocol='docs/MEASURED_IFRAP_PROTOCOL.md',
                  development_movies=15, test_movies=len(metadata[metadata.split == 'reserved_test']),
                  quarantine_movies=8, eligible_test_conditions=len(balanced),
                  condition_balanced_mse=balanced.mean().to_dict(),
                  worst_movie_mse=table.groupby('model').mse.max().to_dict(),
                  comparisons=comparisons, gate_passed=all(v['pass_gate'] for v in comparisons.values()),
                  weighted_loss_identity_max_error=maximum_identity_error,
                  fits=fits, development_drift_sensitivity=sensitivity,
                  independent_biological_validation=False, latent_rate_truth_available=False,
                  novelty='Raw-contrast objective is conventional weighted least squares; no architecture novelty.')
    output.mkdir(exist_ok=True)
    table.to_csv(output/'measured_ifrap_movies.csv', index=False)
    pd.DataFrame(predictions).to_csv(output/'measured_ifrap_predictions.csv', index=False)
    (output/'measured_ifrap.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['fits','development_drift_sensitivity']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset-root', type=Path, required=True)
    args = parser.parse_args()
    run(args.dataset_root, Path(__file__).resolve().parent/'results')
