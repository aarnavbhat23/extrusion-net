"""Propagate finite reference counts into a loading-ratio target interval."""
from loading_target import binomial_interval, drift_bounds


def reference_ratio(intervals, sensitivity_ratio=(1., 1.)):
    """A conservative ratio image, under zero false positives and kappa<=1."""
    ql, qu = sensitivity_ratio
    if not 0 < ql <= qu:
        raise ValueError('Sensitivity-ratio bounds must be positive and ordered')
    (a, b), (c, d) = intervals
    if a <= 0 or c <= 0:
        return None, 'zero_reference_lower'
    low, high = c/b/qu, min(1., d/a/ql)
    if low > high:
        return None, 'incompatible_reference'
    return (low, high), 'ok'


def counted_bounds(retained, initiated, reference_positive, reference_initiated,
                   sensitivity_ratio=(1., 1.), error=.05):
    retention_intervals = [binomial_interval(k, n, error/4) for k, n in zip(retained, initiated)]
    reference_intervals = [binomial_interval(k, n, error/4) for k, n in zip(reference_positive, reference_initiated)]
    ratio, status = reference_ratio(reference_intervals, sensitivity_ratio)
    result = dict(interval=[0., 1.], status=status, certificates=[]) if ratio is None else drift_bounds(retention_intervals, ratio=ratio)
    result.update(retention_intervals=retention_intervals, reference_intervals=reference_intervals,
                  ratio_contract=ratio, sensitivity_ratio=sensitivity_ratio, joint_error_bound=error)
    return result
