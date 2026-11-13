"""Random-delay anchor observations: Markov bridge, hidden capture, external trials."""
import numpy as np


def capture_snapshot(parameters,gamma=.3,removal_factor=1.):
    mu,pl,pr,h=parameters;mu*=removal_factor
    flow=np.zeros((9,9,2,2));flow[0,0,0,0]=1.;result=np.zeros((9,9))
    for a in range(9):
        for b in range(9):
            for cl in range(2):
                for cr in range(2):
                    if flow[a,b,cl,cr]==0:continue
                    left=float(a<8 and not cl);right=float(b<8 and not cr)
                    mass=flow[a,b,cl,cr]/(left+right+mu*(h if cl or cr else 1.)+gamma)
                    result[a,b]+=gamma*mass
                    if left:
                        probability=pl if a+1==3 else 0.
                        flow[a+1,b,0,cr]+=mass*(1-probability);flow[a+1,b,1,cr]+=mass*probability
                    if right:
                        probability=pr if b+1==3 else 0.
                        flow[a,b+1,cl,0]+=mass*(1-probability);flow[a,b+1,cl,1]+=mass*probability
    return np.r_[result.ravel(),max(0.,1-result.sum())]


def markov_snapshot(parameters,gamma=.3,removal_factor=1.):
    """Equivalent to the mapped published recurrence, not a new forward operator."""
    left=np.r_[parameters[:8],0.];right=np.r_[parameters[8:16],0.];mu=parameters[16]*removal_factor
    result=np.zeros((9,9))
    for a in range(9):
        for b in range(9):
            incoming=gamma if a==b==0 else 0.
            if a:incoming+=left[a-1]*result[a-1,b]
            if b:incoming+=right[b-1]*result[a,b-1]
            result[a,b]=incoming/(left[a]+right[b]+mu+gamma)
    return np.r_[result.ravel(),max(0.,1-result.sum())]


def external_snapshots(parameters,gamma,count,seed,removal_factor=1.):
    """One independent censor clock and engine seed per episode; removal is absent."""
    import lefsim
    mu,pl,pr,h=parameters;mu*=removal_factor
    lattice=lefsim.Lattice(17,{'site':[5,11],'captures_from':['right','left'],'strength':[pl,pr]})
    world=lefsim.World(lattice,{'left_pos':[8],'right_pos':[8]})
    model=lefsim.Model(world,lefsim.Rates(step=1.,respawn=mu,respawn_blocked=mu,respawn_captured=mu*h))
    rng=np.random.default_rng(seed+2**40);clocks=rng.exponential(1/gamma,count);counts=np.zeros(82,int)
    for trial,clock in enumerate(clocks):
        simulation=lefsim.Sim(model,seed=int(seed+trial));a=b=0
        while True:
            event=simulation.step()
            if event is None:raise RuntimeError('No removal event')
            time,action,leg=event
            if time>clock:counts[a*9+b]+=1;break
            if action!='Step':counts[-1]+=1;break
            if leg==0:a+=1
            else:b+=1
    return counts


def fit_snapshot(observed,kind):
    """Three fixed starts, identical Hellinger objective; point fits only."""
    from scipy.optimize import least_squares
    import time
    start=time.perf_counter();fits=[];forward=capture_snapshot if kind=='capture' else markov_snapshot
    lower=[.03,0.,0.,.15] if kind=='capture' else [1e-5]*16+[.03]
    upper=[1.5,1.,1.,1.2] if kind=='capture' else [4.]*16+[1.5]
    for mu,p,h in ((.15,.2,.25),(.35,.5,.6),(.7,.8,1.)):
        initial=[mu,p,p,h] if kind=='capture' else np.r_[np.full(16,1.),mu]
        if kind!='capture':initial[[3,11]]=1-p
        objective=lambda x:np.sqrt(np.maximum(forward(x),1e-16))-np.sqrt(observed)
        fit=least_squares(objective,initial,bounds=(lower,upper),max_nfev=300,ftol=1e-10,xtol=1e-10,gtol=1e-10)
        fits.append({'parameters':fit.x.tolist(),'loss':float(np.sum(fit.fun**2)),
                     'success':bool(fit.success),'evaluations':int(fit.nfev)})
    chosen=min(fits,key=lambda x:x['loss'])
    return {**chosen,'starts':fits,'seconds':time.perf_counter()-start}
