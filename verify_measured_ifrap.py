"""Reaggregate original source and replay saved predictions without the fit code."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def verify(root, output):
    result = json.loads((output/'measured_ifrap.json').read_text())
    source_path = root/'source/combined_dataframe.csv'
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == result['source_sha256']
    source = pd.read_csv(source_path)
    metadata = pd.read_csv(root/'extrusion/movies.csv', dtype={'date_raw': str})
    recorded = pd.read_csv(output/'measured_ifrap_predictions.csv')
    metrics = pd.read_csv(output/'measured_ifrap_movies.csv')
    assert set(recorded.movie_id) == set(metadata[metadata.split == 'reserved_test'].movie_id)
    assert not set(metadata[metadata.split == 'development'].date_raw)&set(metadata[metadata.split == 'reserved_test'].date_raw)
    normalized = {}
    for row in metadata.itertuples(index=False):
        g = source[source.filename == row.filename]
        unbleached = g.pivot(index='nucleus', columns='frame', values='mean_list_interp_unbleached').to_numpy().mean(axis=0)
        bleached = g.pivot(index='nucleus', columns='frame', values='mean_list_interp_bleached').to_numpy().mean(axis=0)
        first = g[g.nucleus == g.nucleus.iloc[0]].sort_values('frame')
        denominator = (first.unfrap_cell-first.background_int).to_numpy()
        corrected = (unbleached-bleached)/(denominator/np.mean(denominator[:4]))
        normalized[row.movie_id] = corrected[5:]/corrected[5]
    error = 0.; t = np.arange(245)*.5
    replay_rows = []
    for (movie_id, model), group in recorded.groupby(['movie_id', 'model']):
        group = group.sort_values('frame'); condition = group.condition.iloc[0]
        assert group.frame.tolist() == list(range(5, 250))
        fit = result['fits'][condition][model]; p = fit['parameters']
        if model == 'empirical_mean':
            ids = metadata[(metadata.condition == condition)&(metadata.split == 'development')].movie_id
            predicted = sum(normalized[m] for m in ids)/len(ids)
        elif model == 'single':
            predicted = np.exp(-t*p[0])
        elif model == 'plateau':
            predicted = p[1]+np.exp(-t*p[0])*(1-p[1])
        else:
            predicted = np.exp(-np.outer(t, [p[0], p[0]+p[1]])) @ np.array([p[2], 1-p[2]])
        target = normalized[movie_id]
        error = max(error, float(np.max(abs(group.observed-target))), float(np.max(abs(group.prediction-predicted))))
        mse, mae = float(np.mean((target-predicted)**2)), float(np.mean(abs(target-predicted)))
        record = metrics[(metrics.movie_id == movie_id)&(metrics.model == model)].iloc[0]
        assert abs(mse-record.mse) < 1e-12 and abs(mae-record.mae) < 1e-12
        replay_rows.append(dict(condition=condition, model=model, mse=mse))
    balanced = pd.DataFrame(replay_rows).groupby(['condition', 'model']).mse.mean().groupby('model').mean()
    assert all(abs(v-result['condition_balanced_mse'][k]) < 1e-12 for k,v in balanced.items())
    assert error < 1e-12
    report = dict(status='passed', source_reaggregation=True, saved_prediction_rows=len(recorded),
                  maximum_replay_error=error, model_refitting=False, independent_biology=False)
    (output/'measured_ifrap_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--dataset-root', type=Path, required=True)
    verify(parser.parse_args().dataset_root, Path(__file__).resolve().parent/'results')
