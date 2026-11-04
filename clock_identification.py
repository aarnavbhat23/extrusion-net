"""Conservative shared-clock inversion from two local contact bands."""
import numpy as np
from study_load_probes import local_design,recover_local


def identify_clock(loads,second,third,eta):
    """All observed band errors must be bounded by eta simultaneously."""
    eta=float(eta)+1e-14  # outward numerical padding for bounded inputs
    a,b=recover_local(loads,second);estimate=np.stack([a,b])
    inverse=np.linalg.inv(local_design(loads));radius=eta*np.sum(np.abs(inverse),axis=2).T
    low=np.maximum(0.,estimate-radius);high=np.minimum(1.,estimate+radius)
    if np.any(low>high):return {'empty':True}
    theta_low=0.;theta_high=np.inf;ratios=[]
    for view,q in enumerate(loads):
        def static(g):
            x,y=g
            return q[:-2]*y[:-1]*y[1:]+q[1:-1]*x[:-1]*y[1:]+q[2:]*x[:-1]*x[1:]
        residual_low=third[view]-eta-static(high);residual_high=third[view]+eta-static(low)
        wlow=q[1:-1]*low[0,:-1]*low[1,1:];whigh=q[1:-1]*high[0,:-1]*high[1,1:]
        if np.any(residual_high< -1e-12):return {'empty':True}
        lower=np.divide(np.maximum(0.,residual_low),whigh,out=np.zeros_like(whigh),where=whigh>0)
        upper=np.divide(np.maximum(0.,residual_high),wlow,out=np.full_like(wlow,np.inf),where=wlow>0)
        theta_low=max(theta_low,float(np.max(lower)));theta_high=min(theta_high,float(np.min(upper)))
        denom=q[1:-1]*a[:-1]*b[1:];valid=denom>1e-10
        ratios.extend(((third[view]-static(estimate))[valid]/denom[valid]).tolist())
    slo=max(0. if np.isinf(theta_high) else 1/(1+theta_high),float(np.max(low)))
    shi=1/(1+theta_low)
    if slo>shi+1e-12:return {'empty':True}
    point=float(np.clip(1/(1+max(0.,np.median(ratios))),slo,shi)) if ratios else (slo+shi)/2
    gate_low=low/shi;gate_high=np.minimum(1.,high/max(slo,1e-15))
    gates=np.clip(estimate/max(point,1e-15),0,1)
    return {'empty':False,'survival':point,'survival_interval':[slo,shi],
            'gate_estimate':gates,'gate_low':gate_low,'gate_high':gate_high}
