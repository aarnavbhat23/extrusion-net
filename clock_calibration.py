"""Unknown-speed snapshot inference and exact finite-count clock diagnostics."""
from functools import lru_cache
from itertools import product
import numpy as np
from scipy.linalg import solve
from scipy.optimize import minimize
from scipy.stats import beta

LOW = np.array([.2, .2, .03, .001, .001, .1])
HIGH = np.array([3., 3., 2., .999, .999, 1.5])
STARTS = ((1, 1, .3, .4, .6, .6), (.6, 1.4, .6, .7, .3, 1), (1.4, .6, .15, .2, .8, .3))


def flow_snapshot(parameters, gamma, factor=1.):
    vl, vr, mu, pl, pr, h = parameters
    mu *= factor
    flow = np.zeros((9, 9, 2, 2)); flow[0, 0, 0, 0] = 1.
    result = np.zeros((9, 9))
    for a in range(9):
        for b in range(9):
            for cl in range(2):
                for cr in range(2):
                    if not flow[a, b, cl, cr]:
                        continue
                    left = vl if a < 8 and not cl else 0.
                    right = vr if b < 8 and not cr else 0.
                    mass = flow[a, b, cl, cr]/(left+right+mu*(h if cl or cr else 1.)+gamma)
                    result[a, b] += gamma*mass
                    if left:
                        p = pl if a+1 == 3 else 0.
                        flow[a+1, b, 0, cr] += left*mass*(1-p)
                        flow[a+1, b, 1, cr] += left*mass*p
                    if right:
                        p = pr if b+1 == 3 else 0.
                        flow[a, b+1, cl, 0] += right*mass*(1-p)
                        flow[a, b+1, cl, 1] += right*mass*p
    return np.r_[result.ravel(), max(0., 1-result.sum())]


@lru_cache(None)
def generator_snapshot(parameters, gamma, factor=1.):
    """Separate state-generator solve used for data generation, never fitting."""
    vl, vr, mu, pl, pr, h = parameters
    states = list(product(range(9), range(9), range(2), range(2)))
    index = {state: i for i, state in enumerate(states)}
    matrix = np.zeros((len(states), len(states)))
    for i, (a, b, cl, cr) in enumerate(states):
        if a < 8 and not cl:
            p = pl if a+1 == 3 else 0.
            matrix[i, index[a+1, b, 0, cr]] += vl*(1-p)
            matrix[i, index[a+1, b, 1, cr]] += vl*p
        if b < 8 and not cr:
            p = pr if b+1 == 3 else 0.
            matrix[i, index[a, b+1, cl, 0]] += vr*(1-p)
            matrix[i, index[a, b+1, cl, 1]] += vr*p
        matrix[i, i] = -matrix[i].sum()-factor*mu*(h if cl or cr else 1.)
    source = np.zeros(len(states)); source[index[0, 0, 0, 0]] = gamma
    latent = solve((gamma*np.eye(len(states))-matrix).T, source)
    present = np.zeros(81)
    for probability, (a, b, _, _) in zip(latent, states):
        present[9*a+b] += probability
    return np.r_[present, 1-present.sum()]


def log_likelihood(parameters, observations, conditional):
    value = 0.
    for gamma, counts in observations:
        probabilities = flow_snapshot(parameters, gamma)
        if conditional:
            probabilities = probabilities[:-1]/probabilities[:-1].sum()
        counts = np.asarray(counts)
        observed = counts > 0
        if np.any(probabilities[observed] <= 0):
            return -np.inf
        value += counts[observed] @ np.log(probabilities[observed])
    return float(value)


def fit(observations, conditional):
    size = sum(int(np.sum(c)) for _, c in observations)
    if size == 0:
        raise ValueError('No observations')
    fits = []
    for start in STARTS:
        result = minimize(lambda p: -log_likelihood(p, observations, conditional)/size,
                          start, method='L-BFGS-B', bounds=list(zip(LOW, HIGH)),
                          options={'maxiter': 300, 'ftol': 1e-11})
        fits.append({'parameters': result.x.tolist(), 'log_likelihood': log_likelihood(result.x, observations, conditional),
                     'success': bool(result.success), 'iterations': int(result.nit), 'status': int(result.status)})
    return {**max(fits, key=lambda x: x['log_likelihood']), 'starts': fits}


def odds_interval(zero_count, step_count, alpha):
    """q00/q10 via a conditional binomial CP interval; infinity encoded by None."""
    if min(zero_count, step_count) < 0 or int(zero_count) != zero_count or int(step_count) != step_count:
        raise ValueError('Need nonnegative integer counts')
    if not 0 < alpha < 1:
        raise ValueError('Need alpha in (0,1)')
    n = zero_count+step_count
    lo = beta.ppf(alpha/2, zero_count, step_count+1) if zero_count else 0.
    hi = beta.ppf(1-alpha/2, zero_count+1, step_count) if step_count else 1.
    return [float(lo/(1-lo)), float(hi/(1-hi)) if hi < 1 else None]


def speed_interval(observations, alpha=.05):
    """At least 1-alpha coverage for vL under two known distinct clocks and the family."""
    if len(observations) != 2:
        raise ValueError('Need two clocks')
    (g1, c1), (g2, c2) = sorted(observations, key=lambda x: x[0])
    if g1 <= 0 or g2 <= g1:
        raise ValueError('Need distinct positive clocks')
    first = odds_interval(int(c1[0]), int(c1[9]), alpha/2)
    second = odds_interval(int(c2[0]), int(c2[9]), alpha/2)
    lo_difference = second[0]-(first[1] if first[1] is not None else np.inf)
    hi_difference = (second[1] if second[1] is not None else np.inf)-first[0]
    if hi_difference <= 0:
        return {'bounds': None, 'empty': True, 'odds_intervals': [first, second]}
    bounds = [float((g2-g1)/hi_difference), float((g2-g1)/lo_difference) if lo_difference > 0 else None]
    return {'bounds': bounds, 'empty': False, 'odds_intervals': [first, second]}
