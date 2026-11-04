"""Sharp real-arithmetic output bounds for one unsaturated routing intervention.

These numerical evaluations are not outward-rounded interval certificates.
"""
import numpy as np

from study_identification import (admissible_initial_interval, contact_array,
                                  factor_at_initial, observable_bands)


def symmetric_contacts(reduced):
    size=len(reduced)
    upper=np.zeros((size+1,size+1),dtype=reduced.dtype)
    i,j=np.triu_indices(size)
    upper[i,j+1]=reduced[i,j]
    return upper+upper.T


def cut_endpoints(baseline,site,multiplier=0.):
    """Two contact endpoints for a fixed right-gate multiplier in [0,1]."""
    if not 0<=multiplier<=1:raise ValueError('Multiplier must be in [0,1]')
    q,sums,products=observable_bands(baseline)
    if np.any(q<=0) or np.any(products<=0):
        raise ValueError('Positive loads and interior family products required')
    low,high=admissible_initial_interval(q,sums,products)
    matrices=[];prefix=[]
    for initial in (low,high):
        left,right=factor_at_initial(q,sums,products,initial)
        if not np.allclose(contact_array(q,left,right),baseline,atol=1e-9,rtol=1e-9):
            raise ValueError('Baseline is not an exact feasible routing map')
        prefix.append(float(np.prod(right[:site])))
        right[site]*=multiplier
        matrices.append(symmetric_contacts(contact_array(q,left,right)))
    if not np.all(np.isfinite(prefix)) or np.any(np.asarray(prefix)<=0):
        raise FloatingPointError('Unresolved prefix product; bounds not returned')
    return np.asarray(matrices),(low,high),np.asarray(prefix)


def initial_for_mix(theta,interval,prefix):
    """Invert the common matrix-segment coordinate to the initial right gate."""
    fraction=theta*prefix[0]/((1-theta)*prefix[1]+theta*prefix[0])
    return interval[0]+fraction*(interval[1]-interval[0])


def mixer_outputs(contacts,values,weight,bias=None):
    """One-head production formula: linear messages, log1p mass, affine output."""
    features=np.c_[contacts@values,np.log1p(contacts.sum(axis=1))]
    return features@weight.T+(0 if bias is None else bias)


def sharp_output_bounds(endpoints,values,weight,bias=None):
    """Min/max per output, including the possible interior stationary point.

    Input endpoints already include the fixed positive contact scale beta.
    Each output is a+b*theta+c*log(d+e*theta), theta in [0,1].
    """
    low,high=endpoints
    messages=low@values;difference=(high-low)@values
    a=messages@weight[:,:-1].T+(0 if bias is None else bias)
    b=difference@weight[:,:-1].T
    c=np.broadcast_to(weight[:,-1],a.shape)
    d=np.broadcast_to(1+low.sum(axis=1,keepdims=True),a.shape)
    e=np.broadcast_to((high-low).sum(axis=1,keepdims=True),a.shape)
    stationary=np.zeros_like(a)
    valid=(b!=0)&(e!=0)
    with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
        stationary[valid]=-c[valid]/b[valid]-d[valid]/e[valid]
    stationary=np.where(np.isfinite(stationary)&(stationary>0)&(stationary<1),stationary,0.)
    candidates=np.stack([np.zeros_like(a),np.ones_like(a),stationary],axis=-1)
    evaluated=a[...,None]+b[...,None]*candidates+c[...,None]*np.log(d[...,None]+e[...,None]*candidates)
    lower_index=evaluated.argmin(axis=-1);upper_index=evaluated.argmax(axis=-1)
    return {'lower':evaluated.min(axis=-1),'upper':evaluated.max(axis=-1),
            'lower_theta':np.take_along_axis(candidates,lower_index[...,None],axis=-1)[...,0],
            'upper_theta':np.take_along_axis(candidates,upper_index[...,None],axis=-1)[...,0]}
