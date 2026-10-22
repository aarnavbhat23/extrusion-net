"""Replay original count draws and diagnose conservatism with a supremum exact p."""
import json
from pathlib import Path
import numpy as np
from pulse_order import equal_arm_fisher,robust_fisher_p,robust_threshold

if __name__=='__main__':
    original=json.loads(Path('results/pulse_order.json').read_text());rng=np.random.default_rng(2101)
    report={'refinement_protocol_commit':'437c8b0','fresh_confirmation':False,'null':[],'memory':[]}
    def compare(row,null):
        pa,pb=row['survivals'];n=1000;r=5
        x=rng.binomial(n,pa,2000);y=rng.binomial(n,pb,2000)
        assert float(np.mean(equal_arm_fisher(x,y,n)<=.05))==row['ideal_fisher']
        high=np.maximum(x,y);low=np.minimum(x,y)
        if null:
            ox=np.minimum(n,high+r);oy=np.maximum(0,low-r)
        else:
            overlap=high-low<=2*r;mid=(high+low)//2
            ox=np.where(overlap,mid,high-r);oy=np.where(overlap,mid,low+r)
        rate=float(np.mean(np.abs(ox-oy)/n>robust_threshold(n,n,r,r)))
        assert rate==row['robust_known_corruption']
        return {'degree_balanced':row['degree_balanced'],'robust_fisher_rejection':float(np.mean(robust_fisher_p(ox,oy,n,r)<=.05))}
    for case in range(100):
        # Consume exactly the same frailty-generation draws as the frozen study.
        rng.uniform(.01,1,(32,2));rng.uniform(0,.2,(32,2));rng.dirichlet(np.ones(32))
        for degree in (0,1,2):
            row=original['null'][case*3+degree];report['null'].append({'case':case,**compare(row,True)})
    for row in original['memory']:
        report['memory'].append({'capture':row['capture'],'protection':row['protection'],**compare(row,False)})
    primary=[r for r in report['memory'] if r['degree_balanced']==1]
    report['summary']={'original_count_draws_reproduced':True,
        'primary_exact_robust_power':[r['robust_fisher_rejection'] for r in primary],
        'primary_cases_at_least_80pct_power':sum(r['robust_fisher_rejection']>=.8 for r in primary),
        'primary_max_robust_null_rejection':max(r['robust_fisher_rejection'] for r in report['null'] if r['degree_balanced']==1)}
    Path('results/pulse_test_refinement.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['summary'],indent=2))
