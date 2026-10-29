"""Frozen finite-data comparison and model-mismatch stress screen."""
import json
from pathlib import Path
import time
import numpy as np
from finite_loss import reference,probabilities,likelihood,fit,TIMES

ROOT=Path(__file__).resolve().parent
MECHANISMS=((1.05,.95,.32,.4,.7,.55,.12),(.85,1.15,.55,.65,.25,.8,.2),(1.25,.75,.65,.5,.6,.4,.08))
METHODS=('three_normalized','absolute_retention','scheduled_timestamps')


def conditions(method):
    if method=='three_normalized':return [(1.,.1,None),(1.,.3,None),(2.,.3,None)]
    if method=='absolute_retention':return [(1.,.3,None),(2.,.3,None)]
    return [(f,None,t) for f in (1.,2.) for t in TIMES]


def main():
    rng=np.random.default_rng(20261030);rows=[]
    for stress in ('matched','turnover_error','state_dropout'):
        for case,parameters in enumerate(MECHANISMS):
            extra=.15 if stress=='state_dropout' else 0.
            target=reference(parameters,3.,.3,extra)
            for n in ((6000,24000) if stress=='matched' else (24000,)):
                for repeat in range(3 if stress=='matched' else 2):
                    for method in METHODS:
                        settings=conditions(method);saved=[];training=[];validation=[]
                        for factor,gamma,t in settings:
                            actual=1.8 if stress=='turnover_error' and factor==2 else factor
                            p=reference(parameters,actual,gamma,extra)
                            if gamma is None:p=p[TIMES.index(t)]
                            if stress=='matched':
                                expected=probabilities(parameters,method)[len(saved)]
                                if method=='three_normalized':
                                    np.testing.assert_allclose(p[:-1]/p[:-1].sum(),expected,atol=2e-12,rtol=0)
                                else:np.testing.assert_allclose(p,expected,atol=2e-12,rtol=0)
                            count=n//(2*len(settings))
                            train=rng.multinomial(count,p);holdout=rng.multinomial(count,p)
                            saved.append({'factor':factor,'actual_factor':actual,'gamma':gamma,'time':t,
                                          'training_counts':train.tolist(),'validation_counts':holdout.tolist()})
                            training.append(train[:-1] if method=='three_normalized' else train)
                            validation.append(holdout[:-1] if method=='three_normalized' else holdout)
                        started=time.perf_counter();fitted=fit(method,training);seconds=time.perf_counter()-started
                        prediction=reference(tuple(fitted['parameters']),3.,.3)
                        fitted_validation=likelihood(fitted['parameters'],method,validation)
                        base_validation=likelihood(parameters,method,validation)
                        threshold=fitted_validation+np.log(.05)
                        rows.append({'stress':stress,'case':case,'parameters':parameters,'initiated_episodes':n,
                                     'repeat':repeat,'method':method,'observations':saved,'fit':fitted,'fit_seconds':seconds,
                                     'training_retained':sum(sum(x['training_counts'][:-1]) for x in saved),
                                     'expected_sum_observation_duration':sum(n/len(settings)*(1/g if g else t) for _,g,t in settings),
                                     'holdout_log_likelihood_at_fit':fitted_validation,'holdout_threshold':threshold,
                                     'holdout_log_likelihood_at_base_parameter':base_validation,
                                     'base_parameter_in_split_set':bool(base_validation>=threshold),
                                     'prediction':prediction.tolist(),'target':target.tolist(),
                                     'absence_error':float(abs(prediction[-1]-target[-1])),
                                     'target_total_variation':float(np.abs(prediction-target).sum()/2)})
                        (ROOT/'results/finite_loss_progress.json').write_text(json.dumps({'complete':False,'rows':rows},indent=2)+'\n')
                        print(len(rows),stress,case,n,repeat,method,'error',round(rows[-1]['absence_error'],5),'seconds',round(seconds,1),flush=True)
    report={'protocol_commit':'e4ebb15','seed':20261030,'complete':True,'rows':rows,
            'scope':'Development simulations; split-set membership is not predictive interval width or empirical coverage validation.'}
    (ROOT/'results/finite_loss.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (ROOT/'results/finite_loss_progress.json').unlink()


if __name__=='__main__':main()
