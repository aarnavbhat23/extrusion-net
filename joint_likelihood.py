"""Split likelihood confidence sets; classical universal inference, not new physics."""
import numpy as np
from scipy.optimize import minimize
from occupancy import capture_snapshot
from continuous_occupancy import DOMAIN_LOW, DOMAIN_HIGH

STARTS = ((.15, .2, .2, .25), (.35, .5, .5, .6), (.7, .8, .8, 1.))


def multinomial_log_likelihood(counts, probabilities):
    """Parameter-dependent part only. Preserve exact zero-probability events."""
    counts = np.asarray(counts)
    p = np.asarray(probabilities, float)
    if counts.shape != p.shape or counts.ndim != 1 or np.any(counts < 0) or not np.all(counts == np.floor(counts)):
        raise ValueError('Need aligned nonnegative integer counts')
    if not np.all(np.isfinite(p)) or np.any(p < 0) or not np.isclose(p.sum(), 1., atol=1e-12, rtol=0):
        raise ValueError('Need a probability vector')
    observed = counts > 0
    if np.any(p[observed] == 0):
        return -np.inf
    return float(counts[observed] @ np.log(p[observed]))


def split_counts(counts, rng):
    counts = np.asarray(counts)
    if counts.ndim != 1 or np.any(counts < 0) or not np.all(counts == np.floor(counts)):
        raise ValueError('Need nonnegative integer counts')
    counts = counts.astype(np.int64)
    train = rng.multivariate_hypergeometric(counts, int(counts.sum()) // 2)
    return train, counts - train


def log_likelihood(parameters, observations):
    return sum(multinomial_log_likelihood(counts, capture_snapshot(parameters, gamma))
               for gamma, counts in observations)


def fit_training(observations):
    size = sum(int(c.sum()) for _, c in observations)
    fits = []
    for start in STARTS:
        fit = minimize(lambda p: -log_likelihood(p, observations) / size, start,
                       method='L-BFGS-B', bounds=list(zip(DOMAIN_LOW, DOMAIN_HIGH)),
                       options={'maxiter': 300, 'ftol': 1e-12, 'gtol': 1e-8})
        fits.append({'parameters': fit.x.tolist(), 'log_likelihood': log_likelihood(fit.x, observations),
                     'success': bool(fit.success), 'status': int(fit.status), 'iterations': int(fit.nit)})
    finite = [f for f in fits if np.isfinite(f['log_likelihood'])]
    if not finite:
        raise RuntimeError('No finite training fit; stop instead of replacing the split')
    best = max(finite, key=lambda f: f['log_likelihood'])
    return {**best, 'starts': fits}


def search_witnesses(observations, fitted, threshold):
    """Inner witnesses only. Never return a search range as a confidence interval."""
    candidates, solvers = [], []

    def add(p):
        if np.all(p >= DOMAIN_LOW) and np.all(p <= DOMAIN_HIGH):
            slack = log_likelihood(p, observations) - threshold
            if slack > 1e-6:
                candidates.append({'parameters': np.asarray(p).tolist(), 'slack': slack,
                                   'absence': float(capture_snapshot(p, removal_factor=2.)[-1])})

    starts = (tuple(fitted),) + STARTS
    for start in starts:
        add(np.asarray(start))
    for sign in (1., -1.):
        for start in starts:
            fit = minimize(lambda p: sign * capture_snapshot(p, removal_factor=2.)[-1], start,
                           method='SLSQP', bounds=list(zip(DOMAIN_LOW, DOMAIN_HIGH)),
                           constraints={'type': 'ineq', 'fun': lambda p: log_likelihood(p, observations) - threshold - 1e-5},
                           options={'maxiter': 200, 'ftol': 1e-12})
            solvers.append({'direction': sign, 'start': list(start), 'success': bool(fit.success),
                            'status': int(fit.status), 'iterations': int(fit.nit), 'parameters': fit.x.tolist(),
                            'slack': float(log_likelihood(fit.x, observations) - threshold)})
            add(fit.x)
    if not candidates:
        return {'lower': None, 'upper': None, 'attained_span': None, 'solvers': solvers}
    low = min(candidates, key=lambda p: p['absence'])
    high = max(candidates, key=lambda p: p['absence'])
    return {'lower': low, 'upper': high, 'attained_span': high['absence'] - low['absence'], 'solvers': solvers}
