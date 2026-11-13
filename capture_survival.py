"""Competing-risk identity and rectangle extrema for aggregate absence."""
import numpy as np


def _down(x):return np.nextafter(x,-np.inf)
def _up(x):return np.nextafter(x,np.inf)
def _point(x):return np.asarray(x),np.asarray(x)
def _add(x,y):return _down(x[0]+y[0]),_up(x[1]+y[1])
def _multiply(x,y):return np.maximum(0.,_down(x[0]*y[0])),_up(x[1]*y[1])
def _divide(x,y):return np.maximum(0.,_down(x[0]/y[1])),_up(x[1]/y[0])
def _complement(x):return np.maximum(0.,_down(1-x[1])),np.minimum(1.,_up(1-x[0]))


def _present_interval(mu,pl,pr,h,gamma,factor):
    mu=_multiply(_point(mu),_point(factor));s=_add(mu,_point(gamma))
    one_step=_divide(_point(1.),_add(s,_point(1.)))
    single=_multiply(_multiply(one_step,one_step),one_step)
    two_step=_divide(_point(1.),_add(s,_point(2.)));power=_point(np.ones_like(np.asarray(pl)))
    both=_point(np.zeros_like(np.asarray(pl)))
    for exponent in range(1,6):
        power=_multiply(power,two_step)
        if exponent>=3:both=_add(both,_multiply(_point({3:2.,4:6.,5:12.}[exponent]),power))
    pleft,pright=_point(pl),_point(pr)
    alone=_add(_multiply(pleft,_complement(pright)),_multiply(pright,_complement(pleft)))
    laplace=_add(_multiply(alone,single),_multiply(_multiply(pleft,pright),both))
    laplace=(np.maximum(0.,laplace[0]),np.minimum(1.,laplace[1]))
    before=_divide(_point(gamma),s)
    after=_divide(_point(gamma),_add(_multiply(mu,_point(h)),_point(gamma)))
    return _add(_multiply(before,_complement(laplace)),_multiply(after,laplace))


def absence_bounds(low,high,gamma=.3,removal_factor=1.):
    """Exact mathematical extrema at coupled-order corners; rounded outward."""
    low=np.atleast_2d(low);high=np.atleast_2d(high)
    most_capture=np.where(low[:,3,None]<=1.,high[:,1:3],low[:,1:3])
    least_capture=np.where(high[:,3,None]<=1.,low[:,1:3],high[:,1:3])
    present_max=_present_interval(low[:,0],most_capture[:,0],most_capture[:,1],low[:,3],gamma,removal_factor)
    present_min=_present_interval(high[:,0],least_capture[:,0],least_capture[:,1],high[:,3],gamma,removal_factor)
    return _complement(present_max)[0],_complement(present_min)[1]
