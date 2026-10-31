"""Frozen finite-bank action comparison with retained off-grid failures."""
import itertools,json,time
from pathlib import Path
import numpy as np
from occupancy import capture_snapshot
from bank_acquisition import run,GAMMAS

ROOT=Path(__file__).resolve().parent


def main():
    bank=np.array(list(itertools.product((.15,.35,.7),(.2,.5,.8),(.2,.5,.8),(.25,1.))))
    predictions=np.array([[capture_snapshot(x,gamma) for gamma in GAMMAS] for x in bank])
    targets=np.array([capture_snapshot(x,removal_factor=2.) for x in bank])
    cases=[(mu,pl,pr,h) for mu in (.15,.35,.7) for pl,pr in ((.2,.8),(.5,.5),(.8,.2)) for h in (.25,1.)]
    cases.extend(((.25,.35,.65,.6),(.55,.7,.4,.4),(.9,.4,.7,.8)))
    report={'protocol_commit':'75e4b4a','scope':'Conditional on membership in a supplied finite mechanism bank; not continuous confidence.',
            'bank':bank.tolist(),'gammas':GAMMAS.tolist(),'rows':[]}
    for case_index,parameters in enumerate(cases):
        for repeat in (range(3) if case_index<18 else range(1)):
            for policy_index,policy in enumerate(('round_robin','global_disagreement','target_disagreement')):
                seed=2000000000+case_index*10000000+repeat*1000000+policy_index*100000
                start=time.perf_counter();result=run(parameters,predictions,targets,policy,seed)
                truth=capture_snapshot(parameters,removal_factor=2.)
                covered=False if result['rejected'] else bool(np.all((truth>=np.array(result['lower'])-1e-12)&(truth<=np.array(result['upper'])+1e-12)))
                report['rows'].append({'parameters':parameters,'on_bank':case_index<18,'repeat':repeat,'policy':policy,
                                       'seed':seed,'result':result,'covers_generating_target':covered,'seconds':time.perf_counter()-start})
                (ROOT/'results/bank_acquisition_study.json').write_text(json.dumps(report,indent=2)+'\n')
                print(case_index,repeat,policy,'attained',result['attained'],'rejected',result['rejected'],'cost',result['quoted_cost'],'covered',covered,flush=True)


if __name__=='__main__':main()
