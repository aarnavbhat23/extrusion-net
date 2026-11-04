"""Finite supplied mechanism bank: adaptive observations and conditional envelopes."""
import numpy as np
from scipy.stats import beta
from occupancy import external_snapshots

GAMMAS=np.array([.1,.3,1.,3.])


def confidence_intervals(counts,alpha=.05/(4*8*82)):
    n=int(counts.sum());counts=np.asarray(counts)
    low=np.zeros(len(counts));high=np.ones(len(counts));positive=counts>0;not_all=counts<n
    low[positive]=beta.ppf(alpha/2,counts[positive],n-counts[positive]+1)
    high[not_all]=beta.ppf(1-alpha/2,counts[not_all]+1,n-counts[not_all])
    return low,high


def choose_action(predictions,targets,alive,policy,round_index):
    if policy=='round_robin':return round_index%4
    ids=np.flatnonzero(alive);target=targets[ids]
    target_gap=np.max(np.abs(target[:,None,:]-target[None,:,:]),axis=2)
    pairs=np.triu(np.ones((len(ids),len(ids)),bool),1)
    if policy=='target_disagreement':pairs&=target_gap>.02
    scores=[]
    for i,gamma in enumerate(GAMMAS):
        roots=np.sqrt(predictions[ids,i]);distance=np.maximum(0.,1-roots@roots.T)
        selected=distance[pairs]
        value=(float(np.min(selected)) if policy=='target_disagreement' else float(np.mean(selected))) if len(selected) else 0.
        scores.append(value/(1+1/gamma))
    return int(np.argmax(scores))


def run(parameters,predictions,targets,policy,seed):
    alive=np.ones(len(targets),bool);counts=np.zeros((4,82),int);trace=[];episodes=0;cost=0.
    for index in range(8):
        lower=targets[alive].min(axis=0);upper=targets[alive].max(axis=0)
        if np.max(upper-lower)<=.02:break
        action=choose_action(predictions,targets,alive,policy,index);gamma=GAMMAS[action]
        counts[action]+=external_snapshots(parameters,gamma,1000,seed+index*10000)
        lo,hi=confidence_intervals(counts[action]);alive&=np.all((predictions[:,action]>=lo-1e-12)&(predictions[:,action]<=hi+1e-12),axis=1)
        episodes+=1000;cost+=1000*(1+1/gamma)
        trace.append({'action':int(action),'gamma':float(gamma),'episodes':episodes,'quoted_cost':cost,
                      'retained':np.flatnonzero(alive).tolist(),'counts':counts[action].tolist()})
        if not alive.any():return {'rejected':True,'attained':False,'episodes':episodes,'quoted_cost':cost,'trace':trace}
    lower=targets[alive].min(axis=0);upper=targets[alive].max(axis=0)
    return {'rejected':False,'attained':bool(np.max(upper-lower)<=.02),'episodes':episodes,'quoted_cost':cost,
            'max_target_width':float(np.max(upper-lower)),'lower':lower.tolist(),'upper':upper.tolist(),
            'retained':np.flatnonzero(alive).tolist(),'trace':trace}
