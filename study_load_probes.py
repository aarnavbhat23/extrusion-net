"""Input-design screen with free gates and a shared logistic controller."""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit

from study_identification import contact_array
from study_noisy_identification import contacts_and_jacobian

ROOT=Path(__file__).resolve().parent


def patterns(size,kind,rng,contrast=.8):
    signs=(np.ones(size) if kind=='uniform' else
           rng.choice([-1.,1.],size=size) if kind=='random' else (-1.)**np.arange(size))
    amplitude=0. if kind=='uniform' else contrast
    return .5*np.stack([1+amplitude*signs,1-amplitude*signs])


def local_design(loads):
    """One 2-by-2 design per adjacent contact; columns are left, right."""
    return np.stack([loads[:,1:],loads[:,:-1]],axis=-1).transpose(1,0,2)


def recover_local(loads,observed):
    design=local_design(loads)
    recovered=np.einsum('ijk,ki->ij',np.linalg.pinv(design),observed)
    return recovered[:,0],recovered[:,1]


def controller_output(theta,features,loads,full):
    size=len(features)+1
    gates=expit(features@theta.reshape(2,-1).T)
    gradients=gates*(1-gates)
    chain=np.zeros((2*(size-1),len(theta)))
    chain[:size-1,:features.shape[1]]=gradients[:,0,None]*features
    chain[size-1:,features.shape[1]:]=gradients[:,1,None]*features
    outputs=[];jacobians=[]
    for q in loads:
        if full:
            output,jac=contacts_and_jacobian(np.r_[q,gates[:,0],gates[:,1]],size)
            i,j=np.triu_indices(size)
            outputs.append(output[i!=j]);jacobians.append(jac[i!=j,size:]@chain)
        else:
            outputs.append(q[1:]*gates[:,0]+q[:-1]*gates[:,1])
            jacobians.append(np.c_[q[1:,None]*chain[:size-1,:features.shape[1]],
                                   q[:-1,None]*chain[size-1:,features.shape[1]:]])
    return np.concatenate(outputs),np.concatenate(jacobians)


def fit_controller(features,loads,target,full):
    starts=np.stack([np.zeros(8),np.linspace(-.4,.4,8),np.linspace(.4,-.4,8)])
    fits=[]
    for start in starts:
        cache={}
        def evaluate(theta):
            if 'theta' not in cache or not np.array_equal(theta,cache['theta']):
                cache['theta']=theta.copy()
                output,jac=controller_output(theta,features,loads,full)
                cache['result']=(output-target,jac)
            return cache['result']
        fits.append(least_squares(lambda t:evaluate(t)[0],start,jac=lambda t:evaluate(t)[1],
                                 max_nfev=200,ftol=1e-10,xtol=1e-10,gtol=1e-10))
    best=min(fits,key=lambda f:f.cost)
    return best, max(f.cost for f in fits)-best.cost


def main():
    rng=np.random.default_rng(20261105);sigma=.01
    free=[];shared=[];max_sum_error=0.
    for size in (8,32,128):
        for repetition in range(200):
            left=rng.uniform(.2,.8,size-1);right=rng.uniform(.2,.8,size-1)
            noise=rng.normal(0,sigma,(2,size-1))
            for kind in ('uniform','random','alternating'):
                loads=patterns(size,kind,rng)
                exact=loads[:,1:]*left+loads[:,:-1]*right
                if repetition==0:
                    summed=np.stack([np.diag(contact_array(q,left,right),1) for q in loads])
                    max_sum_error=max(max_sum_error,float(np.max(abs(summed-exact))))
                l,r=recover_local(loads,exact+noise)
                ranks=np.linalg.matrix_rank(local_design(loads))
                free.append({'size':size,'repeat':repetition,'pattern':kind,
                             'gate_mse':float(np.mean(np.r_[l-left,r-right]**2)),
                             'identified_fraction':float(np.mean(ranks==2))})
    for repetition in range(10):
        features=rng.normal(size=(15,4));features[:,0]=1
        test=rng.normal(size=(63,4));test[:,0]=1
        truth=rng.normal(0,.7,8)
        true_test=expit(test@truth.reshape(2,4).T)
        load_sets={k:patterns(16,k,rng) for k in ('uniform','random','alternating')}
        for full in (False,True):
            noise=rng.normal(0,sigma,2*(120 if full else 15))
            for kind,loads in load_sets.items():
                clean,jac=controller_output(truth,features,loads,full)
                fit,spread=fit_controller(features,loads,clean+noise,full)
                predicted=expit(test@fit.x.reshape(2,4).T)
                shared.append({'repeat':repetition,'observations':'full' if full else 'nearest',
                               'pattern':kind,'training_local_rank':int(np.linalg.matrix_rank(jac)),
                               'test_gate_rmse':float(np.sqrt(np.mean((predicted-true_test)**2))),
                               'test_right_cut_rmse':float(.5*np.sqrt(np.mean((predicted[:,0]-true_test[:,0])**2))),
                               'fit_cost':float(fit.cost),'start_cost_spread':float(spread),
                               'optimizer_success':bool(fit.success)})
    drift=[]
    for amplitude in (0.,.005,.02,.05):
        rows=[]
        for repetition in range(200):
            left=rng.uniform(.3,.7,31);right=rng.uniform(.3,.7,31)
            loads=patterns(32,'alternating',rng)
            dl=rng.uniform(-amplitude,amplitude,(2,31))
            dr=rng.uniform(-amplitude,amplitude,(2,31))
            observed=loads[:,1:]*(left+dl)+loads[:,:-1]*(right+dr)
            l,r=recover_local(loads,observed)
            rows.append(float(np.max(np.abs(np.r_[l-left,r-right]))))
        drift.append({'max_gate_drift':amplitude,'largest_recovery_error':max(rows),
                      'deterministic_error_bound':amplitude/.8})
    report={'scope':'Exploratory known-load, fixed-transport input design; not task-label learning.',
            'seed':20261105,'noise_sd':sigma,'free_gate':free,'shared_controller':shared,
            'gate_drift':drift,'direct_sum_max_error':max_sum_error,
            'alternating_predicted_gate_mse':sigma**2/(8*.5**2)*(1+.8**-2)}
    (ROOT/'results/load_probe_study.json').write_text(json.dumps(report,indent=2)+'\n')
    for full in ('nearest','full'):
        for kind in ('uniform','random','alternating'):
            rows=[r for r in shared if r['observations']==full and r['pattern']==kind]
            print(full,kind,'mean test gate RMSE',np.mean([r['test_gate_rmse'] for r in rows]),flush=True)


if __name__=='__main__':main()
