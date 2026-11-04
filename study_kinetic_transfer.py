"""Independent lifetime enumeration: matched and coupled-arm transfer controls."""
import json
from pathlib import Path
import numpy as np
from colored_recovery import choose_colors,recover,scan_messages
from study_load_probes import patterns,recover_local

ROOT=Path(__file__).resolve().parent


def kinetic_contacts(q,left,right,survival,shared):
    """Expected cross-arm reach by enumerating the geometric removal age."""
    m=len(q);d=np.zeros((m+1,m+1))
    if not 0<survival<=1:raise ValueError('Survival must lie in (0,1]')
    for bond,source in enumerate(q):
        reach_l=np.r_[1.,np.cumprod(left[:bond][::-1])]
        reach_r=np.r_[1.,np.cumprod(right[bond:])]
        if not shared:
            reach_l*=survival**np.arange(len(reach_l))
            reach_r*=survival**np.arange(len(reach_r))
            cross=np.outer(reach_l,reach_r)
        else:
            cross=np.zeros((len(reach_l),len(reach_r)));maximum=max(len(reach_l),len(reach_r))-1
            for age in range(maximum+1):
                probability=survival**age*((1-survival) if age<maximum else 1.)
                a=min(age+1,len(reach_l));b=min(age+1,len(reach_r))
                cross[:a,:b]+=probability*np.outer(reach_l[:a],reach_r[:b])
        d[np.ix_(bond-np.arange(len(reach_l)),bond+1+np.arange(len(reach_r)))]+=source*cross
    return d+d.T


def shared_clock_signature(maps,loads):
    m=len(loads[0]);observed=np.stack([np.diag(d,2) for d in maps]);left,right=recover_local(loads,observed)
    q=loads[0]
    static=q[:-2]*right[:-1]*right[1:]+q[1:-1]*left[:-1]*right[1:]+q[2:]*left[:-1]*left[1:]
    excess=(np.diag(maps[0],3)-static)/(q[1:-1]*left[:-1]*right[1:])
    return {'inferred_shared_survival':float(1/(1+np.median(excess))),
            'signature_spread':float(np.ptp(excess))}


def positive_bounds(q,estimate,radius,rho,cut):
    m=len(q);low=np.clip(estimate-radius,0,rho);high=np.clip(estimate+radius,0,rho)
    if cut is not None:low[m-1+cut]=0;high[m-1+cut]=0
    values=np.ones((m+1,1))
    return (scan_messages(q,low[:m-1],low[m-1:],values).ravel(),
            scan_messages(q,high[:m-1],high[m-1:],values).ravel())


def main():
    rng=np.random.default_rng(20261113);rows=[];rho=.6
    for size in (16,32):
        colors,rate=choose_colors(rho,size);values=np.eye(colors)[np.arange(size+1)%colors]
        loads=patterns(size,'alternating',None)
        for repeat in range(5):
            left=rng.uniform(.2,.6,size-1);right=rng.uniform(.2,.6,size-1)
            noise=rng.normal(size=(2,size+1,colors));hold_noise=rng.normal(size=(2,size+1))
            for survival in (1.,.95,.8,.5):
                for shared in (False,True):
                    maps=[kinetic_contacts(q,left,right,survival,shared) for q in loads]
                    clean=np.stack([d@values for d in maps]);signature=shared_clock_signature(maps,loads)
                    hold_q=np.full(size,.5);held=[]
                    for cut in (None,size//2):
                        r=right.copy()
                        if cut is not None:r[cut]=0
                        held.append(kinetic_contacts(hold_q,left,r,survival,shared).sum(axis=1))
                    for sigma in (0.,.0001,.001):
                        estimate,fit=recover(clean+sigma*noise,size,rho,colors,rate)
                        calibration_radius=sigma*np.sqrt(2*np.log(2*2*(size-1)/.025))
                        gate_radius=fit['fixed_point_error_bound']+calibration_radius/(.8*(1-rate))
                        held_radius=sigma*np.sqrt(2*np.log(2*2*(size+1)/.025))
                        rejected=False;errors=[]
                        for j,cut in enumerate((None,size//2)):
                            low,high=positive_bounds(hold_q,estimate,gate_radius,rho,cut)
                            observed=held[j]+sigma*hold_noise[j]
                            rejected|=bool(np.any((observed<low-held_radius-1e-10)|(observed>high+held_radius+1e-10)))
                            predicted=positive_bounds(hold_q,estimate,0.,rho,cut)[0]
                            errors.extend((predicted-held[j]).tolist())
                        effective=survival*np.r_[left,right]
                        rows.append({'size':size,'repeat':repeat,'survival':survival,'shared_clock':shared,'noise_sd':sigma,
                                     'vector_responses':2*colors,'gate_rmse_vs_marginal':float(np.sqrt(np.mean((estimate-effective)**2))),
                                     'heldout_message_rmse':float(np.sqrt(np.mean(np.array(errors)**2))),
                                     'heldout_rejected':rejected,'gate_radius':gate_radius,**signature,
                                     'recovery_converged':fit['converged']})
    # Post-screen diagnostic: a single central source supplies an exact rank witness.
    witnesses=[]
    for rank in (4,8,16,32):
        m=2*rank-1;q=np.zeros(m);q[rank-1]=1.;gates=np.ones(m-1)
        block=kinetic_contacts(q,gates,gates,.8,True)[:rank,rank:][::-1]
        singular=np.linalg.svd(block,compute_uv=False)
        witnesses.append({'block_size':rank,'survival':.8,
                          'relative_tolerance':1e-12,'numerical_rank':int(np.sum(singular>singular[0]*1e-12)),
                          'smallest_singular_value':float(singular[-1]),
                          'kernel_max_error':float(np.max(np.abs(block-.8**np.maximum.outer(np.arange(rank),np.arange(rank)))))})
    result={'scope':'Structural transfer to a separately enumerated kinetic generator, not measured biological validation.',
            'seed':20261113,'calibration_alpha':.025,'heldout_alpha':.025,'rows':rows,
            'post_screen_rank_witnesses':witnesses}
    (ROOT/'results/kinetic_transfer_study.json').write_text(json.dumps(result,indent=2)+'\n')
    for shared in (False,True):
        for sigma in (0.,.0001,.001):
            print(shared,sigma,[(s,sum(x['heldout_rejected'] for x in rows if x['shared_clock']==shared and x['noise_sd']==sigma and x['survival']==s)) for s in (1.,.95,.8,.5)],flush=True)


if __name__=='__main__':main()
