"""Identify transport from vector responses instead of contact-map labels."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit
from study_load_probes import patterns
from study_noisy_identification import contacts_and_jacobian

ROOT=Path(__file__).resolve().parent


def message_response(q,left,right,values):
    """Return D(q) v and its free-gate Jacobian; dense reference for this screen."""
    m=len(q);n=m+1
    contacts,jac=contacts_and_jacobian(np.r_[q,left,right],m)
    jac=jac[:,m:];a,b=np.triu_indices(m)
    responses=[];derivatives=[]
    for v in values:
        output=np.zeros(n);gradient=np.zeros((n,2*(m-1)))
        np.add.at(output,a,contacts*v[b+1]);np.add.at(output,b+1,contacts*v[a])
        np.add.at(gradient,a,jac*v[b+1,None]);np.add.at(gradient,b+1,jac*v[a,None])
        responses.append(output);derivatives.append(gradient)
    return np.concatenate(responses),np.concatenate(derivatives)


def shared_response(theta,features,loads,values):
    gates=expit(features@theta.reshape(2,4).T);p=len(features)
    derivatives=gates*(1-gates);chain=np.zeros((2*p,8))
    chain[:p,:4]=features*derivatives[:,0,None]
    chain[p:,4:]=features*derivatives[:,1,None]
    outputs=[];jacobians=[]
    for q in loads:
        y,j=message_response(q,gates[:,0],gates[:,1],values)
        outputs.append(y);jacobians.append(j@chain)
    return np.concatenate(outputs),np.concatenate(jacobians)


def main():
    rng=np.random.default_rng(20261106);local=[];shared=[]
    for m in (8,16,32,64):
        for repeat in range(10):
            left=rng.uniform(.2,.8,m-1);right=rng.uniform(.2,.8,m-1)
            values=rng.choice([-1.,1.],size=(2,m+1))
            for kind in ('uniform','alternating','random'):
                # Separate stream preserves the original matched mechanisms/probes.
                design_rng=np.random.default_rng(20261107+100*m+repeat)
                loads=patterns(m,kind,design_rng)
                for probes in (1,2):
                    jac=np.concatenate([message_response(q,left,right,values[:probes])[1] for q in loads])
                    s=np.linalg.svd(jac,compute_uv=False);rank=int(sum(s>s[0]*1e-9))
                    local.append({'size':m,'repeat':repeat,'pattern':kind,'values_per_load':probes,
                                  'response_scalars':len(jac),'parameters':len(s),'rank':rank,
                                  'smallest_singular_value':float(s[-1]),
                                  'mean_variance_over_noise_variance':float(np.mean(s**-2)) if rank==len(s) else None})
    # This is a second exploratory screen, not an independent confirmation.
    for repeat in range(10):
        features=rng.normal(size=(15,4));features[:,0]=1
        test=rng.normal(size=(63,4));test[:,0]=1
        truth=rng.normal(0,.7,8);expected=expit(test@truth.reshape(2,4).T)
        values=rng.choice([-1.,1.],size=(2,17));noise=rng.normal(0,.01,68)
        for kind in ('uniform','random','alternating'):
            loads=patterns(16,kind,rng)
            clean,j=shared_response(truth,features,loads,values);target=clean+noise;fits=[]
            for start in (np.zeros(8),np.linspace(-.4,.4,8),np.linspace(.4,-.4,8)):
                cache={}
                def evaluate(theta):
                    if 'theta' not in cache or not np.array_equal(theta,cache['theta']):
                        cache['theta']=theta.copy();y,jac=shared_response(theta,features,loads,values)
                        cache['result']=(y-target,jac)
                    return cache['result']
                fits.append(least_squares(lambda t:evaluate(t)[0],start,jac=lambda t:evaluate(t)[1],
                                         max_nfev=200,ftol=1e-10,xtol=1e-10,gtol=1e-10))
            fit=min(fits,key=lambda f:f.cost);predicted=expit(test@fit.x.reshape(2,4).T)
            shared.append({'repeat':repeat,'pattern':kind,'observed_response_scalars':68,
                           'training_local_rank':int(np.linalg.matrix_rank(j)),
                           'test_gate_rmse':float(np.sqrt(np.mean((predicted-expected)**2))),
                           'test_right_cut_rmse':float(.5*np.sqrt(np.mean((predicted[:,0]-expected[:,0])**2))),
                           'fit_cost':float(fit.cost),'optimizer_success':bool(fit.success),
                           'start_cost_spread':float(max(f.cost for f in fits)-fit.cost)})
    report={'scope':'Exploratory vector-response supervision, known inputs and invariant gates; no task-label claim.',
            'seed':20261106,'noise_sd':.01,'local_information':local,'shared_controller':shared}
    (ROOT/'results/message_probe_study.json').write_text(json.dumps(report,indent=2)+'\n')
    for kind in ('uniform','random','alternating'):
        rows=[r for r in shared if r['pattern']==kind]
        print(kind,'test gate RMSE',np.mean([r['test_gate_rmse'] for r in rows]),flush=True)


if __name__=='__main__':main()
