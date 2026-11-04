"""Run the frozen endpoint-order design screen, retaining all schedule families."""
import json
from pathlib import Path
import numpy as np
from pulse_order import (schedule,mixture_survival,generators,endpoint_survival,event_survival,
                         rejection_rates,robust_threshold)


def main():
    rng=np.random.default_rng(2101);event_rng=np.random.default_rng(2102)
    report={'protocol_commit':'2dd76b6','scope':'Synthetic experiment-design screen, no measured biological outcomes.',
            'n_each_arm':1000,'errors_each_arm':5,'repetitions':2000,'total_duration':4.,
            'robust_threshold':robust_threshold(1000,1000,5,5),'null':[],'memory':[],'counterexamples':[]}
    for case in range(100):
        coeff=np.stack([rng.uniform(.01,1,(32,2)),rng.uniform(0,.2,(32,2))],axis=2)
        weights=rng.dirichlet(np.ones(32))
        for degree in (0,1,2):
            word=schedule(degree);pa=mixture_survival(word,coeff,weights);pb=mixture_survival(1-word,coeff,weights)
            report['null'].append({'case':case,'degree_balanced':degree,'survivals':[pa,pb],
                'contrast':pa-pb,**rejection_rates(pa,pb,1000,5,2000,rng,null=True)})
    for k in (.15,.5,1.5):
        for h in (.1,.4):
            q=generators(k,h)
            for degree in (0,1,2):
                word=schedule(degree);pa=endpoint_survival(word,q);pb=endpoint_survival(1-word,q)
                estimates=[event_survival(w,q,20000,event_rng) for w in (word,1-word)]
                row={'capture':k,'protection':h,'degree_balanced':degree,'word':word.tolist(),
                     'switches_per_arm':int(np.sum(np.diff(word)!=0)),'survivals':[pa,pb],
                     'contrast':pa-pb,'event_estimates':estimates,
                     **rejection_rates(pa,pb,1000,5,2000,rng)}
                report['memory'].append(row)
                print('memory',row,flush=True)
    coeff=np.array([[[.02,0,.05],[.02,0,0]]]);weights=np.ones(1)
    for degree in (0,1,2):
        w=schedule(degree)
        report['counterexamples'].append({'type':'quadratic_age_without_hidden_state','degree_balanced':degree,
            'survivals':[mixture_survival(w,coeff,weights),mixture_survival(1-w,coeff,weights)]})
        q=generators(.5,.1,equal_hazards=True)
        report['counterexamples'].append({'type':'hidden_state_with_equal_event_hazards','degree_balanced':degree,
            'survivals':[endpoint_survival(w,q),endpoint_survival(1-w,q)]})
    primary_null=[r for r in report['null'] if r['degree_balanced']==1]
    primary_memory=[r for r in report['memory'] if r['degree_balanced']==1]
    # 36 engine comparisons, union-bound sampling check at .01.
    epsilon=float(np.sqrt(np.log(2*36/.01)/(2*20000)))
    errors=[abs(p-e) for row in report['memory'] for p,e in zip(row['survivals'],row['event_estimates'])]
    report['summary']={'primary_max_null_contrast':max(abs(r['contrast']) for r in primary_null),
        'primary_max_robust_null_rejection':max(r['robust_known_corruption'] for r in primary_null),
        'primary_cases_at_least_80pct_power':sum(r['robust_known_corruption']>=.8 for r in primary_memory),
        'primary_robust_power':[r['robust_known_corruption'] for r in primary_memory],
        'engine_max_error':max(errors),'engine_simultaneous_tolerance':epsilon}
    s=report['summary'];s['gate_passed']=bool(s['primary_max_null_contrast']<1e-12 and
        s['primary_max_robust_null_rejection']<=.065 and s['primary_cases_at_least_80pct_power']>=4)
    assert max(errors)<=epsilon
    Path('results/pulse_order.json').write_text(json.dumps(report,indent=2)+'\n')
    print('SUMMARY',json.dumps(s,indent=2),flush=True)

if __name__=='__main__':main()
