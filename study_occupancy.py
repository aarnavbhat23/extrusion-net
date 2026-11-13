"""Frozen common-observable comparison; all fits use the same baseline counts."""
import itertools,json,time
from pathlib import Path
import numpy as np
from occupancy import capture_snapshot,markov_snapshot,external_snapshots,fit_snapshot
from dlem_comparator import load_model,snapshot,REVISION,SHA256,URL

ROOT=Path(__file__).resolve().parent


def main():
    import lefsim
    bank=np.array(list(itertools.product((.15,.35,.7),(.2,.5,.8),(.2,.5,.8),(.25,1.))))
    bank_calibration=np.array([capture_snapshot(x) for x in bank])
    cases=[(mu,pl,pr,h) for mu in (.15,.35,.7) for pl,pr in ((.2,.8),(.5,.5),(.8,.2)) for h in (.25,1.)]
    cases.extend(((.25,.35,.65,.6),(.55,.7,.4,.4),(.9,.4,.7,.8)))
    report={'protocol_commit':'75e4b4a','scope':'Random-delay single-motor snapshot, not Hi-C, polymer contacts or new physics.',
            'published_forward':{'revision':REVISION,'sha256':SHA256,'url':URL},'engine_version':lefsim.__version__,
            'calibration_episodes':20000,'target_episodes':20000,'gamma':.3,'removal_factor':2.,'rows':[]}
    model=load_model()
    for index,parameters in enumerate(cases):
        start=time.perf_counter();seed=800000000+index*1000000
        counts=external_snapshots(parameters,.3,20000,seed);observed=counts/20000
        target_counts=external_snapshots(parameters,.3,20000,seed+100000,removal_factor=2.);target=target_counts/20000
        exact=capture_snapshot(parameters,removal_factor=2.)
        fits={kind:fit_snapshot(observed,kind) for kind in ('capture','markov')}
        mirror=markov_snapshot(fits['markov']['parameters'],removal_factor=2.)
        published=snapshot(fits['markov']['parameters'],removal_factor=2.,model=model)
        discrepancy=float(np.max(np.abs(mirror-published)))
        if discrepancy>1e-10:raise AssertionError('Published forward bridge failed')
        bank_index=int(np.argmin(np.sum((np.sqrt(bank_calibration)-np.sqrt(observed))**2,axis=1)))
        predictions={'capture':capture_snapshot(fits['capture']['parameters'],removal_factor=2.),
                     'published_dlem':published,'finite_bank':capture_snapshot(bank[bank_index],removal_factor=2.)}
        calibration_predictions={'capture':capture_snapshot(fits['capture']['parameters']),
                                 'published_dlem':snapshot(fits['markov']['parameters'],model=model),
                                 'finite_bank':bank_calibration[bank_index]}
        row={'parameters':parameters,'on_finite_bank':index<18,'seed':seed,'baseline_counts':counts.tolist(),
             'target_counts':target_counts.tolist(),'fits':fits,'bank_selected':bank[bank_index].tolist(),
             'published_mirror_max_error':discrepancy,'seconds':time.perf_counter()-start,
             'target_predictions':{k:v.tolist() for k,v in predictions.items()},'exact_target':exact.tolist(),
             'target_rmse':{k:float(np.sqrt(np.mean((v-target)**2))) for k,v in predictions.items()},
             'exact_target_rmse':{k:float(np.sqrt(np.mean((v-exact)**2))) for k,v in predictions.items()},
             'baseline_rmse':{k:float(np.sqrt(np.mean((v-observed)**2))) for k,v in calibration_predictions.items()},
             'engine_target_rmse':float(np.sqrt(np.mean((target-exact)**2)))}
        report['rows'].append(row)
        (ROOT/'results/occupancy_study.json').write_text(json.dumps(report,indent=2)+'\n')
        print(index,parameters,'baseline',row['baseline_rmse'],'target',row['target_rmse'],'seconds',row['seconds'],flush=True)


if __name__=='__main__':main()
