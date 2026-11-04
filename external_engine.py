"""Observable bridge to unmodified lefsim episodes and conditional prediction boxes."""
import numpy as np
from scipy.stats import beta
from asynchronous_transport import family_kernel

SOURCE_COMMIT='db4f5dadd28d16d0c76911202ff9d1b6205bd70e'


def make_model(mu,p_left,p_right,mode='matched'):
    import lefsim
    lattice=lefsim.Lattice(17,{'site':[5,11],'captures_from':['right','right' if mode=='eviction' else 'left'],
                              'strength':[p_left,p_right]})
    world=lefsim.World(lattice,{'left_pos':[8],'right_pos':[8]})
    rates=lefsim.Rates(step=1.,respawn=mu,respawn_blocked=mu,respawn_captured=mu/4 if mode=='protected' else mu)
    return lefsim.Model(world,rates)


def episodes(model,count,seed):
    """Independent seeded trials, stopped before the first engine removal/reset."""
    import lefsim
    terminal=np.empty((count,2),np.int16);actions={'Respawn':0,'Evict':0}
    for i in range(count):
        simulation=lefsim.Sim(model,seed=seed+i);left=right=8
        while True:
            event=simulation.step()
            if event is None:raise RuntimeError('Engine exhausted before removal')
            _,action,leg=event
            if action!='Step':actions[action]+=1;break
            if leg==0:left-=1
            else:right+=1
        terminal[i]=[8-left,right-8]
    histogram=np.zeros((9,9),int);np.add.at(histogram,(terminal[:,0],terminal[:,1]),1)
    reach=histogram[::-1,::-1].cumsum(0).cumsum(1)[::-1,::-1]/count
    return reach,actions


def binomial_interval(successes,count,alpha):
    lo=0. if successes==0 else float(beta.ppf(alpha/2,successes,count-successes+1))
    hi=1. if successes==count else float(beta.ppf(1-alpha/2,successes+1,count-successes))
    return lo,hi


def calibrate(reach,count,alpha=.025):
    # Source location, step rate, gate locations and nonbarrier gates are supplied.
    observed=np.array([reach[1,0],reach[4,0],reach[0,4]])
    intervals=np.array([binomial_interval(int(round(x*count)),count,alpha/3) for x in observed])
    slo,shi=intervals[0];s=float(observed[0]);point=np.r_[s,np.clip(observed[1:]/s**4,0,1) if s else np.ones(2)]
    low=np.r_[slo,np.maximum(0.,intervals[1:,0]/shi**4)]
    high=np.r_[shi,np.minimum(1.,intervals[1:,1]/max(slo**4,1e-300))]
    if np.any(low>high):return {'empty':True,'observed':observed.tolist(),'probability_intervals':intervals.tolist()}
    return {'empty':False,'point':point.tolist(),'low':low.tolist(),'high':high.tolist(),
            'observed':observed.tolist(),'probability_intervals':intervals.tolist()}


def predict(parameters,family='asynchronous',intervention='none'):
    s,left,right=parameters
    if intervention in ('left','both'):left=0.
    if intervention=='both':right=0.
    axes=np.arange(9)
    kernel=family_kernel(9,s,family) if s else np.pad(np.ones((1,1)),((0,8),(0,8)))
    return kernel*np.where(axes>=4,left,1.)[:,None]*np.where(axes>=4,right,1.)[None,:]
