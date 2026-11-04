"""Frozen independent-engine transfer, with retained mismatches and observation cost."""
import importlib.metadata,json,time
from pathlib import Path
import numpy as np
from external_engine import make_model,episodes,calibrate,predict,SOURCE_COMMIT

ROOT=Path(__file__).resolve().parent


def run_case(mu,p_left,p_right,mode,seed):
    count=10000;start=time.perf_counter();model=make_model(mu,p_left,p_right,mode)
    calibration,cal_actions=episodes(model,count,seed);fit=calibrate(calibration,count)
    row={'removal':mu,'capture_left':p_left,'capture_right':p_right,'mode':mode,'seed':seed,
         'calibration_episodes':count,'fit':fit,'calibration_actions':cal_actions,'targets':[]}
    if fit['empty']:row['rejected']=True;return row
    eta=float(np.sqrt(np.log(2*81*3/.025)/(2*count)))
    for index,intervention in enumerate(('none','left','both')):
        pl=1. if intervention in ('left','both') else p_left;pr=1. if intervention=='both' else p_right
        observed,actions=episodes(make_model(mu,pl,pr,mode),count,seed+100000*(index+1))
        low=predict(fit['low'],intervention=intervention);high=predict(fit['high'],intervention=intervention)
        errors={family:float(np.sqrt(np.mean((predict(fit['point'],family,intervention)-observed)**2)))
                for family in ('asynchronous','synchronized','independent')}
        rejected=bool(np.any((observed<low-eta-1e-12)|(observed>high+eta+1e-12)))
        row['targets'].append({'intervention':intervention,'episodes':count,'rmse':errors,'rejected':rejected,
                               'mean_model_interval_width':float(np.mean(high-low)),'test_noise_radius':eta,'actions':actions,
                               'observed_reach':observed.tolist(),'predicted_reach':predict(fit['point'],intervention=intervention).tolist(),
                               'lower':low.tolist(),'upper':high.tolist()})
    row['rejected']=any(x['rejected'] for x in row['targets']);row['seconds']=time.perf_counter()-start
    return row


def main():
    import lefsim
    if lefsim.__version__!='0.1.0':raise ValueError('Pinned engine version required')
    report={'scope':'Independent alpha event simulator, restricted reach bridge; not measured biology or polymer contacts.',
            'protocol_commit':'e71840a','engine_version':lefsim.__version__,'engine_source_commit':SOURCE_COMMIT,
            'engine_wheel_sha256_macos_arm64':'afc3013a418fd6da62f3016a52fedf5039f82f33b39d96ca1635eefe9e3cc47c',
            'pyarrow_version':importlib.metadata.version('pyarrow'),'rows':[]}
    for mode in ('matched','protected','eviction'):
        for mu_index,mu in enumerate((.15,.5,1.)):
            strengths=[(.2,.8),(.8,.2)] if mode=='matched' else [(.8,.8)]
            for strength_index,(pl,pr) in enumerate(strengths):
                for repeat in (range(3) if mode=='matched' else range(1)):
                    seed=13000000+1000000*('matched','protected','eviction').index(mode)+100000*mu_index+10000*strength_index+repeat*10000000
                    row=run_case(mu,pl,pr,mode,seed);report['rows'].append(row)
                    (ROOT/'results/external_engine_study.json').write_text(json.dumps(report,indent=2)+'\n')
                    print(mode,mu,pl,pr,repeat,'reject',row['rejected'],'mean_rmse',
                          {k:float(np.mean([x['rmse'][k] for x in row['targets']])) for k in ('asynchronous','synchronized','independent')},flush=True)


if __name__=='__main__':main()
