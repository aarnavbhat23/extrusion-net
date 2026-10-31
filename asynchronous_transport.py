"""Continuous-time reach and local-band timing identification."""
import numpy as np
from study_load_probes import recover_local


def reach_kernel(size,removal,left_rate=1.,right_rate=1.):
    """P(two independent Poisson arms reach a,b before shared exponential death)."""
    if removal<0 or min(left_rate,right_rate)<=0:raise ValueError('Invalid rates')
    h=np.ones((size,size))
    h[:,0]=(left_rate/(left_rate+removal))**np.arange(size)
    h[0,:]=(right_rate/(right_rate+removal))**np.arange(size)
    for a in range(1,size):
        for b in range(1,size):
            h[a,b]=(left_rate*h[a-1,b]+right_rate*h[a,b-1])/(left_rate+right_rate+removal)
    return h


def contacts(q,left,right,kernel):
    """Reach probability readout, not a polymer contact or a final-anchor map."""
    m=len(q);result=np.zeros((m+1,m+1))
    for bond,source in enumerate(q):
        pl=np.r_[1.,np.cumprod(left[:bond][::-1])];pr=np.r_[1.,np.cumprod(right[bond:])]
        cross=source*np.outer(pl,pr)*kernel[:len(pl),:len(pr)]
        result[np.ix_(bond-np.arange(len(pl)),bond+1+np.arange(len(pr)))]+=cross
    return result+result.T


def family_kernel(size,survival,family):
    if family=='asynchronous':return reach_kernel(size,1/survival-1)
    if family=='synchronized':return survival**np.maximum.outer(np.arange(size),np.arange(size))
    if family=='independent':return survival**np.add.outer(np.arange(size),np.arange(size))
    raise ValueError('Unknown timing family')


def fit_local(loads,second,third,family):
    """Point inversion only; no noisy confidence guarantee."""
    a,b=recover_local(loads,second);a=np.clip(a,0,1);b=np.clip(b,0,1);ratios=[]
    for q,c in zip(loads,third):
        static=q[:-2]*b[:-1]*b[1:]+q[1:-1]*a[:-1]*b[1:]+q[2:]*a[:-1]*a[1:]
        weight=q[1:-1]*a[:-1]*b[1:];valid=weight>1e-10
        ratios.extend(((c-static)[valid]/weight[valid]).tolist())
    theta=max(0.,float(np.median(ratios))) if ratios else 0.
    minimum=max(float(max(a.max(),b.max())),1e-4)
    survival=1/(1+theta) if family=='synchronized' else (1-theta)/(1+theta)
    survival=float(np.clip(survival,minimum,1.))
    return {'survival':survival,'left':a/survival,'right':b/survival,'theta_unclipped':theta}


def sampled_contacts(q,left,right,removal,samples,seed):
    """Independent lifetime/event-count sampler; does not call reach_kernel."""
    rng=np.random.default_rng(seed);m=len(q);result=np.zeros((m+1,m+1))
    for bond,source in enumerate(q):
        l=left[:bond][::-1];r=right[bond:]
        allowed_l=np.cumprod(rng.random((samples,len(l)))<l,axis=1).sum(axis=1)
        allowed_r=np.cumprod(rng.random((samples,len(r)))<r,axis=1).sum(axis=1)
        if removal:
            lifetime=rng.exponential(1/removal,samples)
            allowed_l=np.minimum(allowed_l,rng.poisson(lifetime));allowed_r=np.minimum(allowed_r,rng.poisson(lifetime))
        counts=np.zeros((bond+1,m-bond))
        np.add.at(counts,(allowed_l,allowed_r),1)
        reach=counts[::-1,::-1].cumsum(axis=0).cumsum(axis=1)[::-1,::-1]/samples
        result[np.ix_(bond-np.arange(bond+1),bond+1+np.arange(m-bond))]+=source*reach
    return result+result.T


def fit_asymmetric(loads,second,third,fourth,equal=False):
    """Pooled point fit; physical constraints retained, no confidence certificate."""
    from scipy.optimize import least_squares
    a,b=recover_local(loads,second);a=np.clip(a,1e-8,1-1e-8);b=np.clip(b,1e-8,1-1e-8)
    weights=[];residuals=[];design=[];target=[]
    for q,c3,c4 in zip(loads,third,fourth):
        pure3=q[:-2]*b[:-1]*b[1:]+q[1:-1]*a[:-1]*b[1:]+q[2:]*a[:-1]*a[1:]
        weights.extend((q[1:-1]*a[:-1]*b[1:]).tolist());residuals.extend((c3-pure3).tolist())
        pure4=q[:-3]*b[:-2]*b[1:-1]*b[2:]+q[3:]*a[:-2]*a[1:-1]*a[2:]
        design.extend(np.stack([q[1:-2]*a[:-2]*b[1:-1]*b[2:],q[2:-1]*a[:-2]*a[1:-1]*b[2:]],axis=1).tolist())
        target.extend((c4-pure4).tolist())
    weights=np.array(weights);residuals=np.array(residuals);design=np.array(design);target=np.array(target)
    theta_raw=float(weights@residuals/(weights@weights))
    rho_raw=np.linalg.lstsq(design,target,rcond=None)[0]
    umin=float(a.max()/(1-a.max()));vmin=float(b.max()/(1-b.max()))
    maximum=1/(1+(2*max(umin,vmin) if equal else umin+vmin))
    theta=float(np.clip(theta_raw,1e-8,maximum));total=1/theta-1
    u_raw=.5*(total+(rho_raw[1]-rho_raw[0])/theta**2)
    fraction=float(np.clip((u_raw-umin)/max(total-umin-vmin,1e-15),0,1))
    def parameters(x):
        t=x[0];total=1/t-1
        u=total/2 if equal else umin+x[1]*max(0.,total-umin-vmin)
        return t,u,total-u
    def errors(x):
        t,u,v=parameters(x);rho=np.array([1+t+(v+1)*t*t,1+t+(u+1)*t*t])
        return np.r_[weights*t-residuals,design@rho-target]
    initial=[theta] if equal else [theta,fraction]
    fit=least_squares(errors,initial,bounds=([1e-8] if equal else [1e-8,0.],[maximum] if equal else [maximum,1.]),
                      ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=200)
    t,u,v=parameters(fit.x);sl=u/(u+1);sr=v/(v+1)
    return {'left_rate_over_removal':u,'right_rate_over_removal':v,'theta':t,
            'left':a/sl,'right':b/sr,'raw_theta':theta_raw,'raw_directional_coefficients':rho_raw.tolist(),
            'initial_projection':bool(theta!=theta_raw or not 0<=((u_raw-umin)/max(total-umin-vmin,1e-15))<=1),
            'physical_boundary_fit':bool(np.any(fit.active_mask!=0)),'fit_success':bool(fit.success),
            'band_residual_sum_squares':float(np.sum(fit.fun**2))}
