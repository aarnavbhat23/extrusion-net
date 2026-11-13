"""Prespecified fresh off-grid mechanisms; continuous enclosures versus point bank."""
import itertools,json,time
from pathlib import Path
import numpy as np
from occupancy import capture_snapshot,external_snapshots
from continuous_occupancy import enclose,calibration_intervals

ROOT=Path(__file__).resolve().parent


def main():
    cases=((.21,.3,.7,.4),(.45,.65,.35,.8),(.8,.4,.6,.55))
    bank=np.array(list(itertools.product((.15,.35,.7),(.2,.5,.8),(.2,.5,.8),(.25,1.))))
    predictions={gamma:np.array([capture_snapshot(x,gamma) for x in bank]) for gamma in (.3,.1)}
    targets=np.array([capture_snapshot(x,removal_factor=2.) for x in bank]);arrays={}
    report={'protocol_commit':'833a651','scope':'Continuous capture-family/domain enclosure; not arbitrary model robustness or clinical validity.','rows':[]}
    for index,parameters in enumerate(cases):
        for repeat in range(2):
            seed=4000000000+index*1000000+repeat*100000
            observations=[(gamma,external_snapshots(parameters,gamma,10000,seed+j*20000)) for j,gamma in enumerate((.3,.1))]
            target_counts=external_snapshots(parameters,.3,20000,seed+50000,2.)
            start=time.perf_counter();result,lo,hi=enclose(observations);seconds=time.perf_counter()-start
            truth=capture_snapshot(parameters,removal_factor=2.);covered=not result['rejected'] and bool(np.all((truth>=np.array(result['lower'])-1e-12)&(truth<=np.array(result['upper'])+1e-12)))
            contains=bool(np.any(np.all((np.array(parameters)>=lo)&(np.array(parameters)<=hi),axis=1)))
            alive=np.ones(len(bank),bool)
            for gamma,counts in observations:
                ol,oh=calibration_intervals(counts);alive&=np.all((predictions[gamma]>=ol)&(predictions[gamma]<=oh),axis=1)
            finite={'rejected':not bool(alive.any()),'retained':np.flatnonzero(alive).tolist()}
            if alive.any():
                bl=targets[alive].min(axis=0);bh=targets[alive].max(axis=0)
                finite.update({'lower':bl.tolist(),'upper':bh.tolist(),'covers_target':bool(np.all((truth>=bl-1e-12)&(truth<=bh+1e-12)))})
            key=f'case_{index}_{repeat}';arrays[key+'_low']=lo;arrays[key+'_high']=hi
            report['rows'].append({'parameters':parameters,'repeat':repeat,'seed':seed,'observations':[{'gamma':g,'counts':c.tolist()} for g,c in observations],
                                   'target_counts':target_counts.tolist(),'exact_target':truth.tolist(),'continuous':result,'seconds':seconds,
                                   'covers_generating_target':covered,'contains_generating_parameters':contains,'finite_bank':finite,
                                   'engine_target_rmse':float(np.sqrt(np.mean((target_counts/20000-truth)**2)))})
            (ROOT/'results/continuous_occupancy_study.json').write_text(json.dumps(report,indent=2)+'\n')
            np.savez_compressed(ROOT/'results/continuous_occupancy_boxes.npz',**arrays)
            print(index,repeat,'cover',covered,'parameters',contains,'attained',result['attained'],'width',result.get('max_target_width'),
                  'boxes',result.get('retained_boxes'),'seconds',seconds,'finite',finite['rejected'],flush=True)


if __name__=='__main__':main()
