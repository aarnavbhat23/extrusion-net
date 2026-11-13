"""Stress local/global bounds and fit free gates from four noisy responses."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
from deterministic_probes import coded_values,probe_response,recovery_constants,recovery_error_bound
from study_message_probes import message_response

ROOT=Path(__file__).resolve().parent


def main():
    rng=np.random.default_rng(20261108);screens=[];fits=[]
    for size in (8,16,32,64):
        p=2*(size-1);zero=np.zeros(p)
        _,j0=probe_response(zero,size)
        for rho in (0.,.02,.05,.08,.2,.8):
            constants=recovery_constants(rho)
            for repeat in range(10):
                gate=rng.uniform(0,rho,p);other=rng.uniform(0,rho,p)
                random=rng.choice([-1.,1.],size=(2,size+1))
                for design,values in (('coded',coded_values(size)),('random',random)):
                    y,j=probe_response(gate,size,values);y2,_=probe_response(other,size,values)
                    s=np.linalg.svd(j,compute_uv=False);full=bool(s[-1]>s[0]*1e-9)
                    screens.append({'size':size,'rho':rho,'repeat':repeat,'design':design,
                                    'smallest_singular_value':float(s[-1]),'full_rank':full,
                                    'mean_variance_over_noise_variance':float(np.mean(s**-2)) if full else None,
                                    'secant_ratio':float(np.linalg.norm(y-y2)/np.linalg.norm(gate-other)) if rho else None,
                                    'observed_tail_norm':float(np.linalg.norm(j-j0,2)) if design=='coded' else None,
                                    'proved_lower':constants['lower'],'proved_tail':constants['tail']})
    rho=.05
    for size in (8,16,32):
        p=2*(size-1)
        for sigma in (.001,.01):
            for repeat in range(5):
                gates=rng.uniform(0,rho,p);noise=rng.normal(0,sigma,4*(size+1))
                random=rng.choice([-1.,1.],size=(2,size+1))
                hold_q=rng.uniform(.1,.9,size);hold_v=rng.choice([-1.,1.],size=(1,size+1))
                cut=size//2
                true_cut=gates.copy();true_cut[size-1+cut]=0
                held,_=message_response(hold_q,true_cut[:size-1],true_cut[size-1:],hold_v)
                for design,values in (('coded',coded_values(size)),('random',random)):
                    clean,_=probe_response(gates,size,values);target=clean+noise;candidates=[]
                    for fraction in (.1,.5,.9):
                        cache={}
                        def evaluate(g):
                            if 'x' not in cache or not np.array_equal(g,cache['x']):
                                cache['x']=g.copy();y,j=probe_response(g,size,values);cache['value']=(y-target,j)
                            return cache['value']
                        candidates.append(least_squares(lambda g:evaluate(g)[0],np.full(p,fraction*rho),
                                                       jac=lambda g:evaluate(g)[1],bounds=(0,rho),
                                                       max_nfev=200,ftol=1e-10,xtol=1e-10,gtol=1e-10))
                    fit=min(candidates,key=lambda f:f.cost);estimate=fit.x.copy();estimate[size-1+cut]=0
                    predicted,_=message_response(hold_q,estimate[:size-1],estimate[size-1:],hold_v)
                    norm_error=float(np.linalg.norm(fit.x-gates));residual=float(np.linalg.norm(fit.fun))
                    radius=float(sigma*np.sqrt(chi2.ppf(.95,len(noise))))
                    bound=recovery_error_bound(residual,radius,rho) if design=='coded' else None
                    exact_bound=recovery_error_bound(residual,float(np.linalg.norm(noise)),rho) if design=='coded' else None
                    fits.append({'size':size,'rho':rho,'noise_sd':sigma,'repeat':repeat,'design':design,
                                 'gate_rmse':norm_error/np.sqrt(p),'heldout_cut_message_rmse':float(np.sqrt(np.mean((held-predicted)**2))),
                                 'error_l2':norm_error,'residual_l2':residual,'noise_l2':float(np.linalg.norm(noise)),
                                 'gaussian_noise_radius_95':radius,'recovery_radius_95':bound,
                                 'exact_noise_error_bound':exact_bound,
                                 'error_bound_covers':bool(norm_error<=bound+1e-10) if bound is not None else None,
                                 'optimizer_success':bool(fit.success),
                                 'start_cost_spread':float(max(f.cost for f in candidates)-fit.cost)})
                    print(size,sigma,repeat,design,'gate RMSE',round(norm_error/np.sqrt(p),6),flush=True)
    report={'scope':'Restricted-domain global recovery theorem and exploratory matched-generator screen.',
            'seed':20261108,'local_and_secant_checks':screens,'free_gate_fits':fits,
            'strong_barrier_constants':recovery_constants(.05)}
    (ROOT/'results/deterministic_probe_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
