"""Experimental symmetric-power scans for a bounded saturation approximation.

Tensor lifting and polynomial kernel approximations are established techniques.
This file specializes them to the existing two-state extrusion transfer matrix.
It does not implement the saturated kernel exactly or establish a new model class.
"""
import math

import torch


def symmetric_power(a, b, c, degree):
    """Symmetric power of [[a,0],[b,c]] in a monomial basis of size degree+1."""
    rows = []
    zero = torch.zeros_like(a)
    for row in range(degree+1):
        rows.append(torch.stack([
            math.comb(row,column)*a**(degree-row)*b**(row-column)*c**column
            if column<=row else zero for column in range(degree+1)
        ],dim=-1))
    return torch.stack(rows,dim=-2)


def contact_power_scan(load, pass_l, pass_r, values, degree):
    """Compute (D ** degree) @ values; ** means entrywise power, not D @ D.

    Time is linear in sequence length for fixed degree; the dense lifted transfer
    multiply is quadratic in degree. This readable implementation is not a kernel.
    """
    if degree<1:
        raise ValueError('degree must be positive')
    n=load.shape[-1]
    if n==1:
        return values*0
    transfers=symmetric_power(pass_r[..., :-1],load[...,1:],pass_l[...,1:],degree)
    state=values.new_zeros(*values.shape[:-2],degree+1,values.shape[-1])
    initial=torch.nn.functional.one_hot(torch.tensor(0,device=values.device),degree+1).to(values.dtype)
    final=initial.flip(0)
    upper=[None]*n
    for i in range(n-1,-1,-1):
        if i<n-1:
            state=transfers[...,i,:,:]@state
        upper[i]=state[...,-1,:]
        state=state+initial[:,None]*values[...,i,None,:]
    state=torch.zeros_like(state)
    lower=[]
    for i in range(n):
        if i>0:
            state=transfers[...,i-1,:,:].transpose(-1,-2)@state
        lower.append(state[...,0,:])
        state=state+final[:,None]*values[...,i,None,:]
    return torch.stack(upper,dim=-2)+torch.stack(lower,dim=-2)


def saturated_polynomial_scan(load, pass_l, pass_r, values, beta, degree):
    """Taylor approximation of (1-exp(-beta*D)) @ values.

    beta is a scalar or has the shape of load's batch/head dimensions. No error
    claim is made without checking a bounded argument and floating-point error.
    """
    if degree<1:
        raise ValueError('degree must be positive')
    beta=torch.as_tensor(beta,dtype=values.dtype,device=values.device)
    result=torch.zeros_like(values)
    for power in range(1,degree+1):
        term=contact_power_scan(load,pass_l,pass_r,values,power)
        coefficient=(-1)**(power+1)*beta**power/math.factorial(power)
        result=result+coefficient[...,None,None]*term
    return result


def contact_argument_bound(load, pass_l, pass_r, beta):
    """Conservative bound on every nonnegative beta*D entry, for gates in [0,1]."""
    beta=torch.as_tensor(beta,dtype=load.dtype,device=load.device)
    if torch.any((load<0)|(load>1)|(pass_l<0)|(pass_l>1)|(pass_r<0)|(pass_r>1)) or torch.any(beta<0):
        raise ValueError('Need nonnegative beta and gates in [0,1]')
    if load.shape[-1]==1:
        return torch.zeros_like(load[...,0])
    p=torch.maximum(pass_l.amax(-1),pass_r.amax(-1))
    distance=torch.arange(1,load.shape[-1],dtype=load.dtype,device=load.device)
    return beta*load.amax(-1)*(distance*p[...,None]**(distance-1)).amax(-1)


def taylor_entry_error(argument_bound, degree):
    """Exact-arithmetic Taylor remainder bound, excluding numerical roundoff."""
    return argument_bound**(degree+1)/math.factorial(degree+1)
