"""Reproducible encounter-rule screen on the authors' measured CTCF workbook."""
import hashlib
import io
import re
import ssl
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import openpyxl
from scipy.optimize import minimize
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ARCHIVE_URL = ('https://zenodo.org/api/records/7409240/files/Custom%20code%20accompanying%20'
               'Davidson,%20Barth,%20Zaczek%20-%20CTCF%20is%20a%20DNA-tension-dependent%20'
               'barrier%20to%20cohesin-mediated%20DNA%20loop%20extrusion.zip/content')
ARCHIVE_SHA = '31478e8411bdbad98beb81a5c79c0c303b4a7ceffb56f809eefbccb46154ce75'
WORKBOOK_SHA = '2b1fee0cc79bc3e319187b1f6a413b0e34e03be499f142fca7dab5e58ece9d96'
MODELS = ['constant', 'published_exponential', 'force_logistic', 'force_cubic',
          'state_logistic', 'state_boosting']


def fetch_workbook(root=Path('external/measured_ctcf')):
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    archive = root / 'authors_code.zip'
    if not archive.exists():
        ca = '/etc/ssl/cert.pem'
        context = ssl.create_default_context(cafile=ca if Path(ca).exists() else None)
        archive.write_bytes(urllib.request.urlopen(ARCHIVE_URL, context=context, timeout=90).read())
    if hashlib.sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA:
        raise ValueError('Archive checksum differs from frozen release')
    with zipfile.ZipFile(io.BytesIO(archive.read_bytes())) as z:
        content = z.read('Loop extrusion data processing/CTCF_overview.xlsx')
    if hashlib.sha256(content).hexdigest() != WORKBOOK_SHA:
        raise ValueError('Workbook checksum mismatch')
    target = root / 'CTCF_overview.xlsx'; target.write_bytes(content)
    return target


def batch_id(path):
    path = str(path).replace('\\', '/')
    if '20220323_25' in path:
        return '20220323_25'
    dates = re.findall(r'(?<!\d)(2[12]\d{4})(?!\d)', path)
    if not dates:
        raise ValueError(f'Unrecognized experimental batch: {path}')
    return '20' + dates[0]


def load_rows(path, secondary=False):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = []; counts = {}
    for sheet in (['CTCF', 'iCTCF'] if secondary else ['CTCF']):
        reasons = {'worksheet_rows': 0, 'not_boolean': 0, 'wrong_side': 0,
                   'missing_or_invalid_state': 0, 'included': 0}
        for number, r in enumerate(list(wb[sheet].values)[1:], start=2):
            reasons['worksheet_rows'] += 1
            if type(r[3]) is not bool:
                reasons['not_boolean'] += 1; continue
            if r[4] != ('N' if sheet == 'CTCF' else 'C'):
                reasons['wrong_side'] += 1; continue
            if any(type(r[i]) not in (int, float) or not np.isfinite(r[i]) or r[i] < 0
                   for i in (7, 8, 9)) or r[7] == 0:
                reasons['missing_or_invalid_state'] += 1; continue
            path_key = str(r[1]).replace('\\', '/')
            rows.append({'batch': batch_id(path_key), 'molecule': path_key + '|' + str(r[2]),
                         'row': number, 'sheet': sheet, 'y': int(r[3]),
                         'force': float(r[7]), 'elapsed': float(r[8]),
                         'loop': float(r[9]), 'side': int(sheet == 'iCTCF')})
            reasons['included'] += 1
        counts[sheet] = reasons
    wb.close()
    return rows, counts


def weights(rows):
    """Equal group mass; equal molecule mass within each group."""
    out = np.empty(len(rows)); groups = sorted(set(r['batch'] for r in rows))
    for group in groups:
        ids = [i for i, r in enumerate(rows) if r['batch'] == group]
        molecules = {rows[i]['molecule'] for i in ids}
        for molecule in molecules:
            chosen = [i for i in ids if rows[i]['molecule'] == molecule]
            out[chosen] = 1 / (len(groups) * len(molecules) * len(chosen))
    return out


def design(rows, model, secondary):
    force = np.log1p(np.array([r['force'] for r in rows]) / .025)
    columns = [force]
    if model == 'force_cubic': columns += [force ** 2, force ** 3]
    if model.startswith('state_'):
        loop = np.array([r['loop'] for r in rows]) / 10
        columns += [np.log1p(np.array([r['elapsed'] for r in rows]) / 100), loop, force * loop]
    if secondary: columns += [np.array([r['side'] for r in rows])]
    return np.column_stack(columns)


def predict(train, test, model, secondary=False):
    y = np.array([r['y'] for r in train]); w = weights(train)
    neff = len(set(r['molecule'] for r in train))
    if model == 'constant':
        if not secondary:
            return np.full(len(test), (neff * (w @ y) + .5) / (neff + 1))
        probabilities = {}
        for side in (0, 1):
            ids = np.array([r['side'] == side for r in train])
            mass = w[ids].sum()
            effective = len({r['molecule'] for r in train if r['side'] == side})
            probabilities[side] = (effective * (w[ids] @ y[ids]) / mass + .5) / (effective + 1) if mass else .5
        return np.array([probabilities[r['side']] for r in test])
    if model == 'published_exponential':
        f = np.array([r['force'] for r in train]); side = np.array([r['side'] for r in train])
        def probability(theta, force, c):
            return np.clip(np.exp(-force / np.exp(theta[0] + (theta[1] * c if secondary else 0))), 1e-8, 1-1e-8)
        def objective(theta):
            p = probability(theta, f, side)
            return -w @ (y*np.log(p)+(1-y)*np.log1p(-p))
        fits = [minimize(objective, [start, 0] if secondary else [start],
                         bounds=[(-10, 3), (-5, 5)] if secondary else [(-10, 3)],
                         method='L-BFGS-B') for start in (-5., -3., -1.)]
        fit = min(fits, key=lambda result: result.fun)
        if not fit.success: raise RuntimeError(fit.message)
        return probability(fit.x, np.array([r['force'] for r in test]), np.array([r['side'] for r in test]))
    x = design(train, model, secondary); xt = design(test, model, secondary)
    if model == 'state_boosting':
        learner = HistGradientBoostingClassifier(max_iter=100, max_leaf_nodes=7,
                    learning_rate=.05, min_samples_leaf=10, l2_regularization=1,
                    random_state=21, early_stopping=False)
    else:
        scaler = StandardScaler().fit(x, sample_weight=w)
        x, xt = scaler.transform(x), scaler.transform(xt)
        learner = LogisticRegression(C=1, max_iter=2000, tol=1e-10)
    learner.fit(x, y, sample_weight=w*neff)
    return learner.predict_proba(xt)[:, 1]


def screen(rows, secondary=False):
    groups = sorted(set(r['batch'] for r in rows)); output = []
    for batch in groups:
        train = [r for r in rows if r['batch'] != batch]
        test = [r for r in rows if r['batch'] == batch]
        y = np.array([r['y'] for r in test]); w = weights(test)
        fold = {'batch': batch, 'rows': len(test), 'molecules': len(set(r['molecule'] for r in test)),
                'test_force_range': [min(r['force'] for r in test), max(r['force'] for r in test)],
                'train_force_range': [min(r['force'] for r in train), max(r['force'] for r in train)],
                'metrics': {}, 'predictions': {}}
        for model in MODELS:
            p = np.clip(predict(train, test, model, secondary), 1e-8, 1-1e-8)
            fold['metrics'][model] = {'log_loss': float(-w @ (y*np.log(p)+(1-y)*np.log1p(-p))),
                                      'brier': float(w @ (p-y)**2)}
            fold['predictions'][model] = p.tolist()
        fold['outcomes'] = y.tolist(); fold['weights'] = w.tolist()
        fold['source_rows'] = [[r['sheet'], r['row']] for r in test]
        output.append(fold)
    mean = {m: {k: float(np.mean([f['metrics'][m][k] for f in output]))
                for k in ('log_loss', 'brier')} for m in MODELS}
    controls = [m for m in MODELS if m != 'state_logistic']
    best = min(controls, key=lambda m: mean[m]['log_loss'])
    rng = np.random.default_rng(21); ids = rng.integers(len(groups), size=(10000, len(groups)))
    comparisons = {}
    candidate = np.array([f['metrics']['state_logistic']['log_loss'] for f in output])
    for control in controls:
        diff = np.array([f['metrics'][control]['log_loss'] for f in output]) - candidate
        comparisons[control] = {'improvement': float(diff.mean()),
              'bootstrap_95_interval': np.quantile(diff[ids].mean(axis=1), [.025,.975]).tolist(),
              'groups_improved': int(np.sum(diff>0))}
    gain = 1-mean['state_logistic']['log_loss']/mean[best]['log_loss']
    passes = gain>=.1 and comparisons[best]['groups_improved']>len(groups)/2 and all(
        d['bootstrap_95_interval'][0]>0 for d in comparisons.values())
    return {'rows':len(rows),'groups':len(groups),'molecules':len(set(r['molecule'] for r in rows)),
            'mean_metrics':mean,'best_control':best,'relative_log_loss_reduction':gain,
            'comparisons':comparisons,'gate_passed':bool(passes),'folds':output}
