"""Grouped time-reversal diagnostic; no motor-state identification."""
import itertools
import numpy as np


def track_scores(frames, distances, thresholds, lags=(1,2,4)):
    frames=np.asarray(frames);distances=np.asarray(distances,float)
    if frames.ndim!=1 or distances.shape!=frames.shape or not np.isfinite(frames).all():
        raise ValueError('Need finite frame indices and matching distances')
    if not np.equal(frames,np.floor(frames)).all() or len(np.unique(frames))!=len(frames):
        raise ValueError('Need unique integer frames')
    thresholds=np.asarray(thresholds,float)
    if thresholds.shape!=(2,) or not np.isfinite(thresholds).all() or thresholds[0]>=thresholds[1]:
        raise ValueError('Need two distinct finite thresholds')
    valid=np.isfinite(distances);f=frames[valid].astype(int);b=np.digitize(distances[valid],thresholds)
    states=dict(zip(f.tolist(),b.tolist()));scores=[];counts=[]
    for lag in lags:
        if int(lag)!=lag or lag<1:raise ValueError('Need positive integer lag')
        pairs=[(state,states[frame+lag]) for frame,state in states.items() if frame+lag in states]
        counts.append(len(pairs))
        values=[0 if a==c else (1 if (c-a)%3==1 else -1) for a,c in pairs]
        scores.append(float(np.mean(values)) if values else 0.)
    return np.array(scores),counts


def reversal_test(group_scores, seed=2117, draws=65535):
    """Exact sign orbit for <=20 independent reversal-invariant clusters."""
    x=np.asarray(group_scores,float)
    if x.ndim!=2 or len(x)==0 or not np.isfinite(x).all():raise ValueError('Need finite cluster x feature scores')
    norm=np.sqrt((x*x).sum(axis=0));norm[norm==0]=1.
    observed=float(np.max(np.abs(x.sum(axis=0))/norm));n=len(x);hits=0;total=0
    if n<=20:
        signs=itertools.product((-1,1),repeat=n);exact=True;number=2**n
    else:
        rng=np.random.default_rng(seed);signs=(rng.choice((-1,1),size=n) for _ in range(draws));exact=False;number=draws
    while total<number:
        block=np.asarray(list(itertools.islice(signs,min(2048,number-total))))
        stats=np.max(np.abs(block@x)/norm,axis=1)
        hits+=int(np.sum(stats>=observed-1e-12));total+=len(block)
    return {'statistic':observed,'adjusted_p':hits/total if exact else (hits+1)/(total+1),
            'exact':exact,'randomizations':total,'cluster_count':n,
            'mean_scores':x.mean(axis=0).tolist(),'sufficient_clusters':n>=8,
            'gate_passed':bool(n>=8 and (hits/total if exact else (hits+1)/(total+1))<=.05)}
