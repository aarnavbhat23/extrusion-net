"""Four vector responses and a global strong-barrier recovery bound."""
import numpy as np
from study_load_probes import patterns
from study_message_probes import message_response


def coded_values(size):
    """Physical positions: constant and repeated ++++---- values."""
    return np.stack([np.ones(size+1),(-1.)**(np.arange(size+1)//4)])


def probe_response(gates,size,values=None):
    if values is None:values=coded_values(size)
    loads=patterns(size,'alternating',None)
    responses=[];jacobians=[]
    for q in loads:
        y,j=message_response(q,gates[:size-1],gates[size-1:],values)
        responses.append(y);jacobians.append(j)
    return np.concatenate(responses),np.concatenate(jacobians)


def recovery_constants(rho,q0=.5,contrast=.8):
    """Sufficient global secant bounds on the convex gate box [0,rho]."""
    if not (0<=rho<1 and q0>0 and 0<contrast<1):
        raise ValueError('Require 0 <= rho < 1, positive q0 and contrast in (0,1)')
    tail=4*np.sqrt(2)*q0*(1+contrast)*((1-rho)**-3-1)
    return {'lower':4*q0*contrast-tail,'upper':4*q0+tail,'tail':tail}


def recovery_error_bound(residual_norm,noise_norm_bound,rho):
    lower=recovery_constants(rho)['lower']
    if lower<=0:raise ValueError('The sufficient recovery theorem does not apply')
    return (residual_norm+noise_norm_bound)/lower
