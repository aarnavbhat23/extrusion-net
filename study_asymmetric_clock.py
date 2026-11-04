"""Theory-driven unequal-speed follow-up; models use identical observed bands."""
import json
from pathlib import Path
import numpy as np
from asynchronous_transport import reach_kernel,contacts,fit_asymmetric
from study_load_probes import patterns

ROOT=Path(__file__).resolve().parent


def main():
    rows=[]
    for m in (16,64,128):
        loads=patterns(m,'alternating',None)
        for seed in range(3):
            rng=np.random.default_rng(12000+seed+m);gates=rng.uniform(.2,.8,(2,m-1))
            for right_rate in (.5,1.,2.):
                for removal in (.1,.5,1.):
                    h=reach_kernel(m,removal,1.,right_rate);maps=[contacts(q,*gates,h) for q in loads]
                    clean=[np.stack([np.diag(c,k) for c in maps]) for k in (2,3,4)];noise=[rng.normal(size=x.shape) for x in clean]
                    for sigma in (0.,.0001,.001):
                        observed=[x+sigma*z for x,z in zip(clean,noise)]
                        fits={'unequal':fit_asymmetric(loads,*observed),'equal':fit_asymmetric(loads,*observed,equal=True)}
                        errors={k:[] for k in fits}
                        for intervention in ('none','left','both'):
                            actual=gates.copy();estimated={k:np.stack([f['left'],f['right']]) for k,f in fits.items()}
                            for g in [actual,*estimated.values()]:
                                if intervention in ('left','both'):g[0,m//3]=0
                                if intervention=='both':g[1,2*m//3]=0
                            q=np.full(m,.5);truth=contacts(q,*actual,h).sum(axis=1)
                            for k,f in fits.items():
                                kernel=reach_kernel(m,1.,f['left_rate_over_removal'],f['right_rate_over_removal'])
                                prediction=contacts(q,*estimated[k],kernel).sum(axis=1)
                                errors[k].extend((prediction-truth).tolist())
                        fit=fits['unequal'];point=np.stack([fit['left'],fit['right']])
                        row={'size':m,'seed':seed,'right_rate':right_rate,'removal':removal,'noise_sd':sigma,'observed_scalars':sum(x.size for x in observed),
                             'estimated_rate_ratios':[fit['left_rate_over_removal'],fit['right_rate_over_removal']],
                             'true_rate_ratios':[1/removal,right_rate/removal],
                             'rate_ratio_relative_error':float(np.mean(np.abs(np.array([fit['left_rate_over_removal'],fit['right_rate_over_removal']])/[1/removal,right_rate/removal]-1))),
                             'maximum_gate_error':float(np.max(abs(point-gates))),
                             'heldout_rmse':{k:float(np.sqrt(np.mean(np.square(e)))) for k,e in errors.items()},
                             'initial_projection':fit['initial_projection'],'physical_boundary_fit':fit['physical_boundary_fit'],
                             'fit_success':fit['fit_success'],'band_residual_sum_squares':{k:f['band_residual_sum_squares'] for k,f in fits.items()}}
                        rows.append(row)
        (ROOT/'results/asymmetric_clock_study.json').write_text(json.dumps({'scope':'Exploratory unequal-speed follow-up; point fits only.','protocol_commit':'4a910b9','rows':rows},indent=2)+'\n')
    for sigma in (0.,.0001,.001):
        r=[x for x in rows if x['noise_sd']==sigma and x['right_rate']!=1.]
        print(sigma,'unequal',np.mean([x['heldout_rmse']['unequal'] for x in r]),'equal',np.mean([x['heldout_rmse']['equal'] for x in r]),
              'relative rate error',np.mean([x['rate_ratio_relative_error'] for x in r]),'boundary',sum(x['physical_boundary_fit'] for x in r),'of',len(r),flush=True)


if __name__=='__main__':main()
