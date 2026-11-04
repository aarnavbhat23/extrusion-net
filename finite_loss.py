"""Competing-loss likelihoods with scheduled timestamps and independent reference."""
from functools import lru_cache
from itertools import product
import numpy as np
from scipy.sparse import coo_matrix,eye
from scipy.sparse.linalg import spsolve,expm_multiply
from scipy.optimize import minimize
from scipy.stats import poisson
from study_loss_perturbations import snapshot, effective_rates

TIMES=(.5,1.,2.,4.)
RATE=20.
TERMS=int(poisson.ppf(1-1e-14,RATE*max(TIMES)))+1
WEIGHTS=poisson.pmf(np.arange(TERMS)[:,None],RATE*np.array(TIMES)[None,:])
LOW=np.array([.2,.2,.03,.001,.001,.1,.001])
HIGH=np.array([2.5,2.5,1.5,.999,.999,1.5,1.])
STARTS=((1,1,.3,.4,.6,.6,.15),(.6,1.4,.6,.7,.3,1,.3),(1.4,.6,.15,.2,.8,.3,.05))


def recover_unknown_factor(baseline, perturbed, gamma=.3):
    """Exact probabilities only: two absolute snapshots, unknown perturbation.

    Requires positive capture, positive regulated removal, factor != 1,
    and captured susceptibility h != 1. Near these boundaries the algebra
    is ill-conditioned. This is not a finite-sample estimator or certificate.
    """
    rate = gamma / baseline[0]
    vl, vr = rate*baseline[9]/baseline[0], rate*baseline[1]/baseline[0]
    first = effective_rates(baseline, gamma, vl, vr)
    second = effective_rates(perturbed, gamma, vl, vr)
    if first[3] is None or second[3] is None:
        raise ValueError('Captured-state rate unavailable')
    difference = second[0]-first[0]
    if abs(difference) < 1e-12:
        raise ValueError('No resolved removal perturbation')
    h = (second[3]-first[3])/difference
    if abs(1-h) < 1e-10:
        raise ValueError('Equal susceptibilities do not identify the factor')
    mu = (first[0]-first[3])/(1-h)
    delta = first[0]-mu
    factor = 1+difference/mu
    return np.array([vl, vr, mu, first[1], first[2], h, delta, factor])


def scheduled(parameters,factor):
    """100 reachable states; nonnegative uniformization with bounded tail."""
    vl,vr,mu,pl,pr,h,delta=parameters
    rows=[]; cols=[]; data=[]; visible=[]
    for a,b in product(range(10),repeat=2):
        source=10*a+b; visible.append(9*(3 if a==9 else a)+(3 if b==9 else b))
        outgoing=[]
        if a<8:
            p=pl if a==2 else 0.
            outgoing.extend(((10*(a+1)+b,vl*(1-p)),(90+b,vl*p)))
        if b<8:
            p=pr if b==2 else 0.
            outgoing.extend(((10*a+b+1,vr*(1-p)),(10*a+9,vr*p)))
        exit_rate=sum(v for _,v in outgoing)+factor*mu*(h if a==9 or b==9 else 1.)+delta
        if exit_rate>RATE+1e-12:raise ValueError('Uniformization rate does not dominate generator')
        rows.append(source);cols.append(source);data.append(1-exit_rate/RATE)
        for target,rate in outgoing:
            rows.append(target);cols.append(source);data.append(rate/RATE)
    transition=coo_matrix((data,(rows,cols)),shape=(100,100)).tocsr()
    power=np.zeros(100);power[0]=1.;latent=WEIGHTS[0,:,None]*power[None,:]
    for j in range(1,TERMS):
        power=transition@power
        latent+=WEIGHTS[j,:,None]*power[None,:]
    present=np.array([np.bincount(visible,weights=p,minlength=81) for p in latent])
    return np.c_[present,1-present.sum(axis=1)]


@lru_cache(maxsize=256)
def reference(parameters,factor,gamma=None,extra_loss=0.):
    """Independent 324-state generator, never used for fitting."""
    vl,vr,mu,pl,pr,h,delta=parameters
    states=list(product(range(9),range(9),range(2),range(2)))
    index={s:i for i,s in enumerate(states)}; rows=[];cols=[];data=[]
    for i,(a,b,cl,cr) in enumerate(states):
        outgoing=[]
        if a<8 and not cl:
            p=pl if a+1==3 else 0.
            outgoing.extend(((index[a+1,b,0,cr],vl*(1-p)),(index[a+1,b,1,cr],vl*p)))
        if b<8 and not cr:
            p=pr if b+1==3 else 0.
            outgoing.extend(((index[a,b+1,cl,0],vr*(1-p)),(index[a,b+1,cl,1],vr*p)))
        removal=factor*mu*(h if cl or cr else 1.)+delta+extra_loss*(a>=4 or b>=4)
        rows.append(i);cols.append(i);data.append(-sum(v for _,v in outgoing)-removal)
        for j,value in outgoing:rows.append(j);cols.append(i);data.append(value)
    operator=coo_matrix((data,(rows,cols)),shape=(324,324)).tocsc()
    initial=np.zeros(324);initial[0]=1.
    if gamma is not None:
        latent=np.array([spsolve(gamma*eye(324,format='csc')-operator,gamma*initial)])
    else:
        latent=np.array([expm_multiply(operator*t,initial) for t in TIMES])
    mapping=[9*a+b for a,b,_,_ in states]
    present=np.array([np.bincount(mapping,weights=p,minlength=81) for p in latent])
    result=np.c_[present,1-present.sum(axis=1)]
    if result.min() < -1e-12:raise AssertionError('Invalid reference probability')
    result=np.maximum(result,0.);result/=result.sum(axis=1,keepdims=True)
    return result[0] if gamma is not None else result


def probabilities(parameters,method):
    if method=='three_normalized':
        return [snapshot(parameters,(1,f,g),True) for f,g in ((1,.1),(1,.3),(2,.3))]
    if method=='absolute_retention':
        return [snapshot(parameters,(1,f,.3),False) for f in (1,2)]
    if method=='scheduled_timestamps':
        return list(np.concatenate([scheduled(parameters,f) for f in (1,2)]))
    raise ValueError(method)


def likelihood(parameters,method,counts):
    value=0.
    for p,c in zip(probabilities(parameters,method),counts):
        c=np.asarray(c);mask=c>0
        if np.any(p[mask]<=0):return -np.inf
        value+=c[mask]@np.log(p[mask])
    return float(value)


def fit(method,counts):
    size=sum(np.sum(c) for c in counts); outcomes=[]
    for start in STARTS:
        result=minimize(lambda p:-likelihood(p,method,counts)/size,start,method='L-BFGS-B',
                        bounds=list(zip(LOW,HIGH)),options={'maxiter':250,'ftol':1e-10})
        outcomes.append({'parameters':result.x.tolist(),'log_likelihood':likelihood(result.x,method,counts),
                         'success':bool(result.success),'iterations':int(result.nit),'status':int(result.status)})
    return {**max(outcomes,key=lambda x:x['log_likelihood']),'starts':outcomes}
