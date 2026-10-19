"""Intervention selection along the exact free-gate equivalence family."""
import json
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent


def observable_bands(c):
    """Recover identifiable coefficients from the first three contact bands."""
    q=np.diag(c).copy()
    sums=np.diag(c,1)/(q[:-1]*q[1:])
    products=q[1:-1]*sums[:-1]*sums[1:]-np.diag(c,2)/(q[:-2]*q[2:])
    return q,sums,products


def admissible_initial_interval(q,sums,products):
    """Closed feasible interval for r[0] by backward monotone constraints.

    Strictly positive interior inputs are assumed. Boundary points describe the
    closure; exact zero gates require separate degeneracy handling.
    """
    low,high=0.,1.
    for k in range(len(sums)-1,0,-1):
        a=q[k+1]*sums[k];b=q[k+1]*products[k-1]
        if a<=low:raise ValueError('No positive feasible gates')
        new_low=max(0.,products[k-1]*q[k],b/(a-low))
        new_high=min(1.,b/(a-high)) if a>high else 1.
        low,high=new_low,new_high
    low=max(low,q[1]*(sums[0]-1/q[0]))
    high=min(high,q[1]*sums[0])
    if low>high:raise ValueError('No feasible gates')
    return low,high


def factor_at_initial(q,sums,products,initial):
    r=np.empty_like(sums);l=np.empty_like(sums)
    r[0]=initial;l[0]=q[0]*(sums[0]-r[0]/q[1])
    for k in range(1,len(sums)):
        l[k]=products[k-1]*q[k]/r[k-1]
        r[k]=q[k+1]*(sums[k]-l[k]/q[k])
    return l,r


def contact_array(q,left,right):
    size=len(q)
    c=np.diag(q.copy())
    for gap in range(1,size):
        for i in range(size-gap):
            j=i+gap
            c[i,j]=c[i,j-1]*right[j-1]+q[j]*np.prod(left[i:j])
    return c


def gate_family(q,sums,products,points=4001):
    initial=np.linspace(.0001,.9999,points)
    r=np.empty((points,len(sums)));l=np.empty_like(r)
    r[:,0]=initial;l[:,0]=q[0]*(sums[0]-r[:,0]/q[1])
    with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
        for k in range(1,len(sums)):
            l[:,k]=products[k-1]*q[k]/r[:,k-1]
            r[:,k]=q[k+1]*(sums[k]-l[:,k]/q[k])
    valid=np.isfinite(l).all(1)&np.isfinite(r).all(1)&(l>0).all(1)&(l<1).all(1)&(r>0).all(1)&(r<1).all(1)
    return l[valid],r[valid]


def main():
    rng=np.random.default_rng(20261021)
    cases=[];rows=[]
    for length in (8,16,32,64):
        for seed in range(20):
            q=rng.uniform(.4,.9,length-1)
            left=rng.uniform(.4,.7,length-2);right=rng.uniform(.4,.7,length-2)
            c=contact_array(q,left,right)
            # Selection sees only contact-derived observables, not true gates.
            observed_q,sums,products=observable_bands(c)
            l,r=gate_family(observed_q,sums,products)
            if len(l)<2:raise RuntimeError('Grid did not resolve the admissible family')
            # Intervention r[k] <- 0 leaves C[k,k+1] = q[k+1]*l[k].
            responses=l*q[None,1:]
            true_responses=left*q[1:]
            spans=np.ptp(responses,axis=0)
            interval=admissible_initial_interval(observed_q,sums,products)
            endpoint_error=max(abs(interval[0]-r[0,0]),abs(interval[1]-r[-1,0]))
            # A common multiplier on both directional gates at each site has a
            # decomposition-independent effect. Use independent seeded factors.
            factors=np.random.default_rng(999+length*100+seed).uniform(.2,.9,len(right))
            common_a=contact_array(q,l[0]*factors,r[0]*factors)
            common_b=contact_array(q,l[-1]*factors,r[-1]*factors)
            local_sensitivity=np.abs(np.gradient(responses,r[:,0],axis=0)[len(r)//2])
            site_rng=np.random.default_rng(1234+length*100+seed)
            sites={'boundary':0,'middle':len(right)//2,'largest_response_span':int(np.argmax(spans)),
                   'midpoint_sensitivity':int(np.argmax(local_sensitivity)),
                   'uniform_random_site':int(site_rng.integers(len(right)))}
            equivalence_error=max(np.max(np.abs(contact_array(q,l[i],r[i])-c)) for i in (0,len(l)//2,len(l)-1))
            inverse=np.linalg.inv(c)
            offband=inverse.copy()
            for k in range(3):
                ids=np.arange(len(q)-k);offband[ids,ids+k]=0
            cases.append({'length':length,'seed':seed,'family_grid_points':len(l),
                          'max_contact_difference':float(equivalence_error),
                          'inverse_offband_max':float(np.max(np.abs(offband))),
                          'analytic_initial_interval':list(map(float,interval)),
                          'grid_endpoint_error':float(endpoint_error),
                          'common_barrier_counterfactual_difference':float(np.max(np.abs(common_a-common_b))),
                          'max_intervention_span':float(spans.max()),'selected_site':sites['largest_response_span']})
            # Same additive measurement noise across strategies; sigma is absolute
            # contact units, not a biological instrument model. Baseline map exact.
            for sigma in (.001,.01,.05):
                noise=rng.normal(0,sigma,100)
                for strategy,k in sites.items():
                    observations=true_responses[k]+noise
                    nearest=np.abs(responses[:,k,None]-observations[None,:]).argmin(axis=0)
                    prediction=responses[nearest]
                    rmse=np.sqrt(np.mean((prediction-true_responses[None,:])**2,axis=1))
                    rows.append({'length':length,'seed':seed,'sigma':sigma,'strategy':strategy,'site':k,
                                 'heldout_intervention_rmse':float(np.sqrt(np.mean((prediction[:,np.arange(len(right))!=k]-true_responses[np.arange(len(right))!=k])**2))),
                                 'all_intervention_rmse':float(np.sqrt(np.mean(rmse**2))),
                                 'parameter_rmse':float(np.sqrt(np.mean(np.concatenate([l[nearest]-left,r[nearest]-right],axis=1)**2)))})
    summaries=[]
    for length in (8,16,32,64):
        for sigma in (.001,.01,.05):
            for strategy in ('boundary','middle','largest_response_span','midpoint_sensitivity','uniform_random_site'):
                selected=[x for x in rows if x['length']==length and x['sigma']==sigma and x['strategy']==strategy]
                summaries.append({'length':length,'sigma':sigma,'strategy':strategy,
                                  'mean_heldout_intervention_rmse':float(np.mean([x['heldout_intervention_rmse'] for x in selected])),
                                  'median_parameter_rmse':float(np.median([x['parameter_rmse'] for x in selected]))})
    report={'seed':20261021,'cases':cases,'rows':rows,'summary':summaries,
            'scope':'80 synthetic mechanisms; exact baseline contacts; one noisy single-gate intervention. No actual-network or biological claim.'}
    (ROOT/'results/identification_study.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'cases':len(cases),'max_equivalence_error':max(x['max_contact_difference'] for x in cases),
                      'max_inverse_offband':max(x['inverse_offband_max'] for x in cases),
                      'summary_sigma_001':[x for x in summaries if x['sigma']==.01]},indent=2))


if __name__=='__main__':main()
