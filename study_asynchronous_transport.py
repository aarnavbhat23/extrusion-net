"""Frozen timing-family discrimination and held-out intervention study."""
import json
from pathlib import Path
import numpy as np
from asynchronous_transport import contacts,family_kernel,fit_local,sampled_contacts
from study_load_probes import patterns

ROOT=Path(__file__).resolve().parent


def main():
    report={'scope':'Conditional reach model, not polymer contacts or measured biological transfer.',
            'protocol_commit':'0c2d447','rows':[],'sampler_checks':[]}
    for m in (16,64):
        loads=patterns(m,'alternating',None)
        for seed in range(3):
            rng=np.random.default_rng(10200+seed+m);gates=rng.uniform(.2,.6,(2,m-1))
            for s in (.5,.8,.95):
                for family in ('asynchronous','synchronized'):
                    kernel=family_kernel(m,s,family);maps=[contacts(q,*gates,kernel) for q in loads]
                    clean=[np.stack([np.diag(c,k) for c in maps]) for k in (2,3,4)]
                    noises=[rng.normal(size=c.shape) for c in clean]
                    for sigma in (0.,.0001,.001,.01):
                        observed=[c+sigma*z for c,z in zip(clean,noises)]
                        fits={k:fit_local(loads,*observed[:2],k) for k in ('asynchronous','synchronized')}
                        scores={}
                        for k,f in fits.items():
                            predicted=[contacts(q,f['left'],f['right'],family_kernel(m,f['survival'],k)) for q in loads]
                            scores[k]=float(np.mean((np.stack([np.diag(c,4) for c in predicted])-observed[2])**2))
                        selected=min(scores,key=scores.get);errors={k:[] for k in fits};errors['static_oracle']=[]
                        for intervention in ('none','left','both'):
                            actual=gates.copy();estimated={k:np.stack([f['left'],f['right']]) for k,f in fits.items()};effective=s*gates.copy()
                            for g in [actual,effective,*estimated.values()]:
                                if intervention in ('left','both'):g[0,m//3]=0
                                if intervention=='both':g[1,2*m//3]=0
                            q=np.full(m,.5);truth=contacts(q,*actual,kernel).sum(axis=1)
                            for k,f in fits.items():
                                prediction=contacts(q,*estimated[k],family_kernel(m,f['survival'],k)).sum(axis=1)
                                errors[k].extend((prediction-truth).tolist())
                            errors['static_oracle'].extend((contacts(q,*effective,np.ones((m,m))).sum(axis=1)-truth).tolist())
                        rmse={k:float(np.sqrt(np.mean(np.square(e)))) for k,e in errors.items()}
                        row={'size':m,'seed':seed,'true_family':family,'true_survival':s,'noise_sd':sigma,
                             'observed_scalars':sum(x.size for x in observed),'fitting_scalars':sum(x.size for x in observed[:2]),
                             'selection_scalars':observed[2].size,'selected_family':selected,'selected_correct':selected==family,
                             'fourth_band_mean_squared_residual':scores,'heldout_rmse':rmse,'selected_heldout_rmse':rmse[selected],
                             'estimated_survivals':{k:f['survival'] for k,f in fits.items()}}
                        report['rows'].append(row)
        (ROOT/'results/asynchronous_transport_study.json').write_text(json.dumps(report,indent=2)+'\n')
    m=8;gates=np.random.default_rng(10400).uniform(.3,.8,(2,m-1));q=np.full(m,.5)
    for seed,mu in enumerate((.1,.5,2.)):
        expected=contacts(q,*gates,family_kernel(m,1/(1+mu),'asynchronous'))
        sampled=sampled_contacts(q,*gates,mu,200000,10401+seed)
        # Each edge is a weighted sum of independent per-source Bernoulli means.
        eta=np.sqrt(np.sum(q*q)*np.log(2*(m+1)**2/.01)/(2*200000))
        report['sampler_checks'].append({'removal':mu,'samples_per_source':200000,'max_absolute_error':float(np.max(abs(expected-sampled))),
                                         'simultaneous_hoeffding_radius':float(eta),'all_within_bound':bool(np.all(abs(expected-sampled)<=eta))})
    (ROOT/'results/asynchronous_transport_study.json').write_text(json.dumps(report,indent=2)+'\n')
    for sigma in (0.,.0001,.001,.01):
        rows=[r for r in report['rows'] if r['noise_sd']==sigma];a=[r for r in rows if r['true_family']=='asynchronous']
        print(sigma,'selection',sum(r['selected_correct'] for r in rows),'/',len(rows),
              'async', {k:np.mean([r['heldout_rmse'][k] for r in a]) for k in ('asynchronous','synchronized','static_oracle')},
              'selected',np.mean([r['selected_heldout_rmse'] for r in a]),flush=True)
    print(report['sampler_checks'])


if __name__=='__main__':main()
