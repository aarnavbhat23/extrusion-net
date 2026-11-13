"""Noisy-map stress test of identifiable routing and intervention selection."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from study_identification import (admissible_initial_interval, contact_array,
                                  factor_at_initial, observable_bands)

ROOT = Path(__file__).resolve().parent


def contacts_and_jacobian(parameters, size):
    """Forward sensitivity recurrence, with independent load/left/right columns."""
    q = parameters[:size]
    left = parameters[size:2*size-1]
    right = parameters[2*size-1:]
    p = len(parameters)
    c = np.diag(q.copy())
    jac = np.zeros((size, size, p))
    diagonal = np.arange(size)
    jac[diagonal, diagonal, diagonal] = 1
    products = np.ones(size)
    product_jac = np.zeros((size, p))
    for gap in range(1, size):
        i = np.arange(size-gap)
        j = i+gap
        old_product = products[:len(i)].copy()
        products[:len(i)] *= left[j-1]
        product_jac[:len(i)] *= left[j-1, None]
        product_jac[i, size+j-1] += old_product
        c[i, j] = c[i, j-1]*right[j-1] + q[j]*products[:len(i)]
        jac[i, j] = jac[i, j-1]*right[j-1, None] + q[j, None]*product_jac[:len(i)]
        jac[i, j, 2*size+j-2] += c[i, j-1]
        jac[i, j, j] += products[:len(i)]
    ids = np.triu_indices(size)
    return c[ids], jac[ids]


def fit_map(observed, initial, intervention=None):
    size = len(observed)
    target = observed[np.triu_indices(size)]
    cache = {}
    def evaluate(parameters):
        if 'x' in cache and np.array_equal(parameters, cache['x']):
            return cache['value']
        predicted, jac = contacts_and_jacobian(parameters, size)
        residual = predicted-target
        if intervention is not None:
            site, value = intervention
            response = parameters[site+1]*parameters[size+site]
            gradient = np.zeros(len(parameters))
            gradient[site+1] = parameters[size+site]
            gradient[size+site] = parameters[site+1]
            residual = np.r_[residual, response-value]
            jac = np.vstack([jac, gradient])
        cache['x'] = parameters.copy()
        cache['value'] = residual, jac
        return residual, jac
    lower = np.r_[np.full(size, .05), np.full(2*(size-1), .005)]
    upper = np.r_[np.full(size, 2.), np.full(2*(size-1), .995)]
    return least_squares(lambda x:evaluate(x)[0], np.clip(initial,lower+1e-8,upper-1e-8),
                         jac=lambda x:evaluate(x)[1], bounds=(lower,upper),
                         max_nfev=200, ftol=1e-9, xtol=1e-9, gtol=1e-9)


def physical_family(c):
    q, sums, products = observable_bands(c)
    if np.min(products) <= 0:
        raise ValueError('Positive interior products required')
    lo, hi = admissible_initial_interval(q,sums,products)
    # Avoid poles at exact endpoints in finite precision.
    initials = [lo+1e-7*(hi-lo), (lo+hi)/2, hi-1e-7*(hi-lo)]
    factors = [factor_at_initial(q,sums,products,x) for x in initials]
    if any(not np.isfinite(np.r_[l,r]).all() or min(l.min(),r.min()) < -1e-5 or
           max(l.max(),r.max()) > 1+1e-5 for l,r in factors):
        raise ValueError('Numerically unresolved feasible family')
    return q, factors


def local_design_scores(parameters, size):
    """Local optimal design with baseline noise and an unconstrained null mode."""
    _, jac = contacts_and_jacobian(parameters,size)
    _, singular, vh = np.linalg.svd(jac,full_matrices=False)
    if singular[-2] < singular[0]*1e-10:
        raise ValueError('More than one unresolved local direction')
    null = vh[-1]
    covariance = (vh[:-1].T / singular[:-1]**2) @ vh[:-1]
    response_jac = np.zeros((size-1,len(parameters)))
    sites = np.arange(size-1)
    response_jac[sites,sites+1] = parameters[size:2*size-1]
    response_jac[sites,size+sites] = parameters[1:size]
    projected = response_jac @ covariance @ response_jac.T
    null_response = response_jac @ null
    scores=[]
    for site in sites:
        if abs(null_response[site])<1e-12:
            scores.append(np.inf)
            continue
        ratio=null_response/null_response[site]
        variance=np.diag(projected)-2*ratio*projected[:,site]+ratio**2*(1+projected[site,site])
        scores.append(float(np.mean(variance[sites!=site])))
    return np.asarray(scores)


def noise_aware_site(parameters, size):
    return int(np.argmin(local_design_scores(parameters,size)))


def robust_design_site(c):
    q,sums,products=observable_bands(c)
    low,high=admissible_initial_interval(q,sums,products)
    scores=[]
    for fraction in np.linspace(.05,.95,9):
        left,right=factor_at_initial(q,sums,products,low+fraction*(high-low))
        scores.append(local_design_scores(np.r_[q,left,right],len(q)))
    return int(np.argmin(np.max(scores,axis=0)))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--noise-aware',action='store_true')
    parser.add_argument('--robust-design',action='store_true')
    args=parser.parse_args()
    if args.robust_design:args.noise_aware=True
    seed=20261030 if args.robust_design else (20261028 if args.noise_aware else 20261026)
    rng = np.random.default_rng(seed)
    rows = []
    for size in (8,16,32):
        for case in range(10):
            q = rng.uniform(.4,.9,size)
            left, right = rng.uniform(.4,.7,(2,size-1))
            true_c = contact_array(q,left,right)
            truth = q[1:]*left
            baseline_noise = rng.normal(size=true_c.shape)
            calibration_noise = float(rng.normal())
            random_site = int(rng.integers(size-1))
            for sigma in (.001,.01,.05):
                observed = np.triu(true_c+sigma*baseline_noise)
                try:
                    physical_family(observed)
                    raw_feasible = True
                except ValueError:
                    raw_feasible = False
                starts = [np.r_[np.clip(np.diag(observed),.051,1.999),
                                np.full(size-1,a),np.full(size-1,b)] for a,b in ((.45,.65),(.65,.45))]
                fits = [fit_map(observed,start) for start in starts]
                fit = min(fits,key=lambda x:x.cost)
                fitted_c = contact_array(fit.x[:size],fit.x[size:2*size-1],fit.x[2*size-1:])
                row = {'size':size,'case':case,'sigma':sigma,'raw_band_family_feasible':raw_feasible,
                       'baseline_converged':bool(fit.success),'baseline_nfev':int(fit.nfev),
                       'baseline_map_rmse':float(np.sqrt(np.mean((fitted_c[np.triu_indices(size)]-true_c[np.triu_indices(size)])**2))),
                       'baseline_cost':float(fit.cost),'strategies':{}}
                try:
                    fq, family = physical_family(fitted_c)
                except ValueError as error:
                    row['failure'] = str(error)
                    rows.append(row)
                    continue
                responses = np.array([fq[1:]*l for l,r in family])
                span = np.abs(responses[-1]-responses[0])
                sites={'boundary':0,'largest_response_span':int(np.argmax(span)),'uniform_random_site':random_site}
                if args.noise_aware:sites['local_noise_aware']=noise_aware_site(fit.x,size)
                if args.robust_design:sites['grid_robust_local']=robust_design_site(fitted_c)
                for strategy,site in sites.items():
                    measurement = truth[site]+sigma*calibration_noise
                    joint_starts=[fit.x] if not args.robust_design else [np.r_[fq,l,r] for l,r in family]
                    joint_fits=[fit_map(observed,start,(site,measurement)) for start in joint_starts]
                    joint=min(joint_fits,key=lambda x:x.cost)
                    prediction = joint.x[1:size]*joint.x[size:2*size-1]
                    heldout = np.arange(size-1)!=site
                    row['strategies'][strategy] = {
                        'site':site,'converged':bool(joint.success),'nfev':int(joint.nfev),
                        'joint_cost_range':float(max(f.cost for f in joint_fits)-min(f.cost for f in joint_fits)),
                        'all_site_rmse':float(np.sqrt(np.mean((prediction-truth)**2))),
                        'heldout_rmse':float(np.sqrt(np.mean((prediction[heldout]-truth[heldout])**2))),
                        'uncalibrated_heldout_rmse':float(np.sqrt(np.mean((responses[1,heldout]-truth[heldout])**2)))}
                rows.append(row)
        print('completed matrix size',size,flush=True)
    summary=[]
    for size in (8,16,32):
        for sigma in (.001,.01,.05):
            group=[r for r in rows if r['size']==size and r['sigma']==sigma]
            valid=[r for r in group if 'failure' not in r]
            item={'size':size,'sigma':sigma,'cases':len(group),'family_failures':len(group)-len(valid),
                  'raw_band_feasible_fraction':float(np.mean([r['raw_band_family_feasible'] for r in group])),
                  'baseline_nonconvergence':sum(not r['baseline_converged'] for r in group),'strategies':{}}
            for strategy in ('boundary','largest_response_span','uniform_random_site') + (('local_noise_aware',) if args.noise_aware else ()) + (('grid_robust_local',) if args.robust_design else ()):
                values=[r['strategies'][strategy] for r in valid]
                item['strategies'][strategy]={
                    'mean_all_site_rmse':float(np.mean([v['all_site_rmse'] for v in values])) if values else None,
                    'mean_heldout_rmse':float(np.mean([v['heldout_rmse'] for v in values])) if values else None,
                    'mean_uncalibrated_rmse':float(np.mean([v['uncalibrated_heldout_rmse'] for v in values])) if values else None,
                    'nonconvergence':sum(not v['converged'] for v in values)}
            paired=[r['strategies']['largest_response_span']['heldout_rmse']-
                    r['strategies']['boundary']['heldout_rmse'] for r in valid]
            item['span_minus_boundary_mean']=float(np.mean(paired)) if paired else None
            item['span_beats_boundary_cases']=sum(v<0 for v in paired)
            summary.append(item)
    report={'seed':seed,'scope':'Exploratory noisy baseline and calibration; constrained least squares, synthetic mechanisms only.',
            'rows':rows,'summary':summary}
    (ROOT/'results'/('robust_identification_study.json' if args.robust_design else ('noise_aware_identification_study.json' if args.noise_aware else 'noisy_identification_study.json'))).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
