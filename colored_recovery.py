"""Spaced value probes and globally contractive directional-gate recovery."""
import numpy as np
from study_load_probes import patterns


def scan_messages(q,left,right,values):
    """NumPy reference of the production four-recurrence unsaturated scan."""
    n=len(q)+1;load=np.r_[0.,q];pl=np.r_[0.,left,0.];pr=np.r_[0.,right,0.]
    output=np.zeros_like(values,dtype=float);Q=np.zeros(values.shape[1]);U=Q.copy()
    for k in range(n-1,-1,-1):
        if k+1<n:U=load[k+1]*Q+pl[k+1]*U
        output[k]=U;Q=values[k]+pr[k]*Q
    G=np.zeros(values.shape[1]);H=G.copy()
    for i in range(1,n):
        G=values[i-1]+pl[i-1]*G;H=load[i]*G+pr[i-1]*H;output[i]+=H
    return output


def alias_contraction(rho,colors,contrast=.8):
    if not (0<=rho<1 and colors>=5 and 0<contrast<1):
        raise ValueError('Require rho in [0,1), at least five colors and positive contrast below one')
    if rho==0:return 0.
    r=rho**colors
    moment1=r/(1-r);moment2=r/(1-r)**2;moment3=r*(1+r)/(1-r)**3
    plus=colors**2*moment3+3*colors*moment2+2*moment1
    minus=(colors**2*moment3-5*colors*moment2+6*moment1)/rho**4
    return float((1+contrast)/(2*contrast)*(plus+minus))


def choose_colors(rho,size,target=.5):
    if not (0<target<1):raise ValueError('Contraction target must be in (0,1)')
    for colors in range(5,size+1):
        rate=alias_contraction(rho,colors)
        if rate<=target:return colors,rate
    return size+1,0.  # Unique colors for all physical positions: no aliases.


def colored_response(gates,size,colors):
    values=np.eye(colors)[np.arange(size+1)%colors]
    return np.stack([scan_messages(q,gates[:size-1],gates[size-1:],values)
                     for q in patterns(size,'alternating',None)])


def local_inverse(responses,size,colors):
    rows=np.arange(size-1);selected=responses[:,rows,(rows+2)%colors]
    sums=(selected[0]+selected[1])  # 2*q0 = 1
    differences=(selected[0]-selected[1])/(.8*(-1.)**rows)
    return np.r_[(sums-differences)/2,(sums+differences)/2]


def recover(observed,size,rho,colors,rate,tolerance=1e-10,max_iterations=100):
    """Project a preconditioned residual step; rate must bound alias contraction."""
    if colors==size+1:expected=0.
    else:expected=alias_contraction(rho,colors)
    if not (0<=expected<=rate<1):raise ValueError('A valid contraction bound below one is required')
    estimate=np.zeros(2*(size-1));updates=[]
    for iteration in range(1,max_iterations+1):
        residual=local_inverse(observed-colored_response(estimate,size,colors),size,colors)
        updated=np.clip(estimate+residual,0,rho);step=float(np.max(abs(updated-estimate)))
        updates.append(step);estimate=updated
        if step<=tolerance:break
    residual=local_inverse(observed-colored_response(estimate,size,colors),size,colors)
    fixed_residual=float(np.max(abs(np.clip(estimate+residual,0,rho)-estimate)))
    return estimate,{'iterations':iteration,'converged':step<=tolerance,'step_history':updates,
                     'fixed_point_residual':fixed_residual,'fixed_point_error_bound':fixed_residual/(1-rate)}
