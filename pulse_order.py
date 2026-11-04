"""Endpoint-only history diagnostics with explicit age/frailty nulls."""
import math
import numpy as np
from scipy.linalg import expm
from scipy.stats import hypergeom


def schedule(degree):
    if int(degree)!=degree or not 0<=degree<=10:raise ValueError('Need degree 0..10')
    return np.array([i.bit_count()%2 for i in range(2**(int(degree)+1))],int)


def moments(word,degree,total=4.):
    word=np.asarray(word,int);step=total/len(word);starts=np.arange(len(word))*step
    return np.array([[np.sum(((starts+step)**(j+1)-starts**(j+1))[word==u])/(j+1)
                      for j in range(degree+1)] for u in (0,1)])


def mixture_survival(word,coefficients,weights,total=4.):
    coefficients=np.asarray(coefficients,float);weights=np.asarray(weights,float)
    if coefficients.ndim!=3 or coefficients.shape[1]!=2 or len(weights)!=len(coefficients):
        raise ValueError('Need molecule by condition by polynomial coefficients')
    if np.any(coefficients<0) or np.any(weights<0) or not np.isclose(weights.sum(),1):
        raise ValueError('This numerical generator uses nonnegative coefficients and unit mixture mass')
    exposure=moments(word,coefficients.shape[2]-1,total)
    return float(weights@np.exp(-np.sum(coefficients*exposure[None,:,:],axis=(1,2))))


def generators(capture,protection,equal_hazards=False):
    rates=[(.2,capture,.05),(1.2,.05,capture)];answer=[]
    for base,on,off in rates:
        hazard=np.array([base,base if equal_hazards else base*protection])
        answer.append(np.array([[-on-hazard[0],on],[off,-off-hazard[1]]]))
    return np.array(answer)


def endpoint_survival(word,q,total=4.):
    p=np.array([1.,0.]);step=total/len(word)
    for u in word:p=p@expm(q[u]*step)
    return float(p.sum())


def event_survival(word,q,n,rng,total=4.):
    """Separate Gillespie implementation with an explicit absorbing outcome."""
    step=total/len(word);survived=0
    for _ in range(n):
        state=0;alive=True
        for u in word:
            elapsed=0.
            while alive:
                rate=-q[u,state,state];wait=rng.exponential(1/rate)
                if elapsed+wait>=step:break
                elapsed+=wait;other=1-state
                if rng.random()<q[u,state,other]/rate:state=other
                else:alive=False
            if not alive:break
        survived+=alive
    return survived/n


def equal_arm_fisher(x,y,n):
    """Exact two-sided Fisher p for equal arm sizes (hypergeometric symmetry)."""
    return np.minimum(1.,2*hypergeom.cdf(np.minimum(x,y),2*n,x+y,n))


def robust_threshold(n_a,n_b,errors_a=0,errors_b=0,alpha=.05):
    if min(n_a,n_b)<=0 or min(errors_a,errors_b)<0 or not 0<alpha<1:
        raise ValueError('Positive sample sizes and valid error/alpha budgets required')
    return math.sqrt(.5*(1/n_a+1/n_b)*math.log(2/alpha))+errors_a/n_a+errors_b/n_b


def rejection_rates(pa,pb,n,r,repetitions,rng,null=False):
    x=rng.binomial(n,pa,repetitions);y=rng.binomial(n,pb,repetitions)
    ideal=np.mean(equal_arm_fisher(x,y,n)<=.05)
    if null:
        high=np.maximum(x,y);low=np.minimum(x,y)
        observed=(np.minimum(n,high+r)-np.maximum(0,low-r))/n
    else:observed=np.maximum(0,np.abs(x-y)-2*r)/n
    robust=np.mean(observed>robust_threshold(n,n,r,r))
    return {'ideal_fisher':float(ideal),'robust_known_corruption':float(robust)}


def robust_fisher_p(x,y,n,errors):
    """Supremum two-sided equal-arm Fisher p over bounded count corrections."""
    high=np.maximum(x,y);low=np.minimum(x,y)
    closest_high=np.maximum(0,high-errors);closest_low=np.minimum(n,low+errors)
    return np.where(closest_high<=closest_low,1.,equal_arm_fisher(closest_high,closest_low,n))
