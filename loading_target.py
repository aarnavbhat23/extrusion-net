"""Exact interval images for a common-loading, state-independent-loss target."""
import numpy as np
from scipy.optimize import linprog
from scipy.stats import beta


def binomial_interval(k, n, error):
    low = 0. if k == 0 else float(beta.ppf(error/2, k, n-k+1))
    high = 1. if k == n else float(beta.ppf(1-error/2, k+1, n-k))
    return low, high


def linear_extreme(c, matrix, rhs, bounds):
    result = linprog(c, A_ub=np.asarray(matrix), b_ub=np.asarray(rhs), bounds=bounds, method='highs')
    if result.status == 2:
        return None
    if result.status == 3:
        return dict(status='unbounded')
    if not result.success:
        raise RuntimeError(result.message)
    return dict(value=float(result.fun), primal=result.x.tolist(),
                inequality_dual=result.ineqlin.marginals.tolist(),
                lower_dual=result.lower.marginals.tolist(), upper_dual=result.upper.marginals.tolist())


def snapshot_bounds(intervals, factor=(1.8, 2.2), loading_max=1., target=3.):
    if any(lo <= 0 for lo, _ in intervals):
        return dict(interval=[0., 1.], status='zero_probability_lower', certificates=[])
    ul, uu = 1/intervals[0][1], 1/intervals[0][0]
    vl, vu = 1/intervals[1][1], 1/intervals[1][0]
    tl, tu = (target-1)/(factor[1]-1), (target-1)/(factor[0]-1)
    # Variables U,V,W; W/(V-U) is the target/calibration multiplier ratio.
    matrix = [[1., -1., 0.], [-tl, tl, -1.], [tu, -tu, 1.], [-(target-1), 0., 1.]]
    rhs = [0., 0., 0., -(target-1)/loading_max]
    bounds = [(max(ul, 1/loading_max), uu), (vl, vu), (0., None)]
    if any(lo is not None and hi is not None and lo > hi for lo, hi in bounds):
        return dict(interval=[0., 1.], status='incompatible', certificates=[])
    c = np.array([1., 0., 1.])
    lower = linear_extreme(c, matrix, rhs, bounds)
    upper = linear_extreme(-c, matrix, rhs, bounds)
    if lower is None or upper is None:
        return dict(interval=[0., 1.], status='incompatible', certificates=[])
    zl, zu = lower['value'], -upper['value']
    return dict(interval=[max(0., 1-1/zl), min(1., 1-1/zu)], status='ok',
                problem=dict(matrix=matrix, rhs=rhs, bounds=bounds),
                certificates=[dict(objective=c.tolist(), **lower), dict(objective=(-c).tolist(), **upper)])


def timestamp_bounds(intervals, factor=(1.8, 2.2), loading=(.6, 1.), target=3.):
    # L=log(loading), A=baseline hazard/gamma, D=treated increment/gamma,
    # W=target increment/gamma. Bounds on exp(L)/(1+A+W) are conservative.
    tl, tu = (target-1)/(factor[1]-1), (target-1)/(factor[0]-1)
    matrix = [[0., 0., tl, -1.], [0., 0., -tu, 1.], [0., -(target-1), 0., 1.]]
    rhs = [0., 0., 0.]
    rows = [[1., -1., 0., 0.], [1., -2., 0., 0.], [1., -1., -1., 0.], [1., -2., -2., 0.]]
    for (lo, hi), row in zip(intervals, rows):
        matrix.append(row); rhs.append(float(np.log(hi)))
        if lo > 0:
            matrix.append((-np.array(row)).tolist()); rhs.append(float(-np.log(lo)))
    bounds = [(float(np.log(loading[0])), float(np.log(loading[1]))), (0., None), (0., None), (0., None)]
    certificates = []
    for c in [[1., 0., 0., 0.], [-1., 0., 0., 0.], [0., 1., 0., 1.], [0., -1., 0., -1.]]:
        solution = linear_extreme(c, matrix, rhs, bounds)
        if solution is None:
            return dict(interval=[0., 1.], status='incompatible', certificates=[])
        if solution.get('status') == 'unbounded':
            return dict(interval=[0., 1.], status='unbounded', certificates=[])
        certificates.append(dict(objective=c, **solution))
    lmin, lmax = certificates[0]['value'], -certificates[1]['value']
    smin, smax = certificates[2]['value'], -certificates[3]['value']
    presence_low = np.exp(lmin)/(1+smax)
    presence_high = np.exp(lmax)/(1+smin)
    return dict(interval=[max(0., 1-presence_high), min(1., 1-presence_low)], status='ok',
                problem=dict(matrix=matrix, rhs=rhs, bounds=bounds), certificates=certificates,
                sharp_timestamp_interval=False)


def witness(parameters, loading_max=1., gamma=.3, target=3.):
    u, v, w = parameters
    a = gamma*(loading_max*u-1)
    mu = loading_max*gamma*w/(target-1)
    factor = None if mu == 0 else 1+(v-u)*loading_max*gamma/mu
    return dict(loading=loading_max, mu=mu, delta=a-mu, factor=factor,
                target_absence=1-loading_max*gamma/(gamma+a+(target-1)*mu))


def drift_bounds(intervals, ratio=(.8, 1.), factor=(1.8, 2.2), loading_max=1., target=3.):
    if not 0 < ratio[0] <= ratio[1] <= 1:
        raise ValueError('The endpoint construction requires loading ratio <=1')
    if any(lo <= 0 for lo, _ in intervals):
        return dict(interval=[0., 1.], status='zero_probability_lower', certificates=[])
    ul, uu = 1/intervals[0][1], 1/intervals[0][0]
    vl, vu = 1/intervals[1][1], 1/intervals[1][0]
    tl, tu = (target-1)/(factor[1]-1), (target-1)/(factor[0]-1)
    # U,V,J,W; J=kappa*V, with linear ratio and hazard constraints.
    matrix = [[0., ratio[0], -1., 0.], [0., -ratio[1], 1., 0.],
              [1., 0., -1., 0.], [-tl, 0., tl, -1.],
              [tu, 0., -tu, 1.], [-(target-1), 0., 0., 1.]]
    rhs = [0., 0., 0., 0., 0., -(target-1)/loading_max]
    bounds = [(max(ul, 1/loading_max), uu), (vl, vu), (0., None), (0., None)]
    if bounds[0][0] > bounds[0][1]:
        return dict(interval=[0., 1.], status='incompatible', certificates=[])
    c = np.array([1., 0., 0., 1.])
    lower = linear_extreme(c, matrix, rhs, bounds)
    upper = linear_extreme(-c, matrix, rhs, bounds)
    if lower is None or upper is None:
        return dict(interval=[0., 1.], status='incompatible', certificates=[])
    zl, zu = lower['value'], -upper['value']
    return dict(interval=[max(0., 1-1/zl), min(1., 1-1/zu)], status='ok',
                problem=dict(matrix=matrix, rhs=rhs, bounds=bounds),
                certificates=[dict(objective=c.tolist(), **lower), dict(objective=(-c).tolist(), **upper)])


def drift_witness(parameters, loading_max=1., gamma=.3, target=3.):
    u, v, j, w = parameters
    result = witness([u, j, w], loading_max, gamma, target)
    result['loading_ratio'] = j/v
    return result
