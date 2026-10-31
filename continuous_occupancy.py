"""Continuous-domain interval enclosure; retained boxes are never replaced by points."""
import numpy as np
from scipy.stats import beta

DOMAIN_LOW=np.array([.03,0.,0.,.15])
DOMAIN_HIGH=np.array([1.5,1.,1.,1.2])


def down(x):return np.nextafter(x,-np.inf)
def up(x):return np.nextafter(x,np.inf)


def interval_snapshot(low,high,gamma=.3,removal_factor=1.,aggregate=False):
    """Outward elementary operations, with nonnegativity and simplex constraints."""
    low=np.atleast_2d(low);high=np.atleast_2d(high);count=len(low)
    if np.any(low>high) or np.any(low[:,[0,3]]<=0) or np.any(low[:,1:3]<0) or np.any(high[:,1:3]>1):
        raise ValueError('Invalid parameter box')
    fl=np.zeros((count,9,9,2,2));fu=fl.copy();fl[:,0,0,0,0]=fu[:,0,0,0,0]=1.
    rl=np.zeros((count,9,9));ru=rl.copy()
    ml=np.maximum(0.,down(low[:,0]*removal_factor));mu=up(high[:,0]*removal_factor)
    for a in range(9):
        for b in range(9):
            for cl in range(2):
                if cl and a!=3:continue
                for cr in range(2):
                    if cr and b!=3:continue
                    left=a<8 and not cl;right=b<8 and not cr
                    dl=np.maximum(0.,down(ml*low[:,3])) if cl or cr else ml
                    du=up(mu*high[:,3]) if cl or cr else mu
                    steps=int(left)+int(right)
                    dl=down(down(dl+steps)+gamma);du=up(up(du+steps)+gamma)
                    tl=np.maximum(0.,down(fl[:,a,b,cl,cr]/du));tu=up(fu[:,a,b,cl,cr]/dl)
                    rl[:,a,b]=np.maximum(0.,down(rl[:,a,b]+down(gamma*tl)))
                    ru[:,a,b]=up(ru[:,a,b]+up(gamma*tu))
                    for moving,axis in ((left,0),(right,1)):
                        if not moving:continue
                        na,nb=(a+1,b) if axis==0 else (a,b+1);arrival=na if axis==0 else nb
                        pl=low[:,1+axis] if arrival==3 else np.zeros(count)
                        pu=high[:,1+axis] if arrival==3 else np.zeros(count)
                        for captured in (0,1):
                            ll=pl if captured else np.maximum(0.,down(1-pu))
                            uu=pu if captured else up(1-pl)
                            nc,nr=(captured,cr) if axis==0 else (cl,captured)
                            fl[:,na,nb,nc,nr]=np.maximum(0.,down(fl[:,na,nb,nc,nr]+down(tl*ll)))
                            fu[:,na,nb,nc,nr]=up(fu[:,na,nb,nc,nr]+up(tu*uu))
    rl=rl.reshape(count,81);ru=ru.reshape(count,81);sl=np.zeros(count);su=np.zeros(count)
    for i in range(81):sl=np.maximum(0.,down(sl+rl[:,i]));su=up(su+ru[:,i])
    lower=np.c_[rl,np.maximum(0.,down(1-su))];upper=np.c_[ru,np.minimum(1.,up(1-sl))]
    if aggregate:
        from capture_survival import absence_bounds
        al,au=absence_bounds(low,high,gamma,removal_factor)
        lower[:,-1]=np.maximum(lower[:,-1],al);upper[:,-1]=np.minimum(upper[:,-1],au)
    return np.maximum(0.,lower),np.minimum(1.,upper)


def calibration_intervals(counts,assays=2):
    counts=np.asarray(counts);n=int(counts.sum());alpha=.05/(assays*82)
    low=np.zeros(82);high=np.ones(82);positive=counts>0;not_all=counts<n
    low[positive]=beta.ppf(alpha/2,counts[positive],n-counts[positive]+1)
    high[not_all]=beta.ppf(1-alpha/2,counts[not_all]+1,n-counts[not_all])
    return np.maximum(0.,low-1e-12),np.minimum(1.,high+1e-12)


def enclose(observations,max_boxes=4096,max_refinements=32,precision=.02,aggregate=False):
    low=DOMAIN_LOW[None].copy();high=DOMAIN_HIGH[None].copy();trace=[]
    intervals=[(gamma,*calibration_intervals(counts,len(observations))) for gamma,counts in observations]
    for iteration in range(max_refinements+1):
        keep=np.ones(len(low),bool)
        for gamma,observed_low,observed_high in intervals:
            predicted_low,predicted_high=interval_snapshot(low,high,gamma,aggregate=aggregate)
            keep&=np.all((predicted_low<=observed_high)&(predicted_high>=observed_low),axis=1)
        low=low[keep];high=high[keep]
        if not len(low):return {'rejected':True,'attained':False,'trace':trace},low,high
        tl,tu=interval_snapshot(low,high,.3,2.,aggregate=aggregate);lower=tl.min(axis=0);upper=tu.max(axis=0)
        width=float(np.max(upper-lower));trace.append({'iteration':iteration,'retained_boxes':len(low),'max_target_width':width})
        if width<=precision or len(low)>=max_boxes or iteration==max_refinements:break
        split_count=min(len(low),512,max_boxes-len(low))
        selected=np.argsort(-np.max(tu-tl,axis=1),kind='stable')[:split_count]
        dimensions=np.argmax((high[selected]-low[selected])/(DOMAIN_HIGH-DOMAIN_LOW),axis=1)
        mid=low[selected].copy();new_high=high[selected].copy();new_low=low[selected].copy()
        for index,(box,dimension) in enumerate(zip(selected,dimensions)):
            value=low[box,dimension]+(high[box,dimension]-low[box,dimension])/2
            new_high[index,dimension]=value;new_low[index,dimension]=value
        retain=np.ones(len(low),bool);retain[selected]=False
        low=np.vstack([low[retain],mid,new_low]);high=np.vstack([high[retain],new_high,high[selected]])
    return {'rejected':False,'attained':width<=precision,'max_target_width':width,'lower':lower.tolist(),'upper':upper.tolist(),
            'retained_boxes':len(low),'trace':trace},low,high
