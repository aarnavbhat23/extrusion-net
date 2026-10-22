"""Post-outcome numerical witnesses; lower bounds on ambiguity, not global optima."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,minimize
from occupancy import capture_snapshot
from continuous_occupancy import DOMAIN_LOW,DOMAIN_HIGH,calibration_intervals

ROOT=Path(__file__).resolve().parent


def main():
    source=json.loads((ROOT/'results/continuous_occupancy_study.json').read_text());rows=[]
    for row in source['rows']:
        observations=[(x['gamma'],np.array(x['counts'])) for x in row['observations']]
        intervals=[calibration_intervals(counts) for _,counts in observations]
        def slack(parameters):
            return np.concatenate([np.r_[capture_snapshot(parameters,g)-lo,hi-capture_snapshot(parameters,g)]
                                   for (g,_),(lo,hi) in zip(observations,intervals)])
        def constraints(parameters):
            result=slack(parameters)
            inset=np.concatenate([np.r_[(lo>0)*1e-8,(hi<1)*1e-8] for lo,hi in intervals])
            return result-inset
        def residual(parameters):
            return np.concatenate([np.sqrt(capture_snapshot(parameters,g))-np.sqrt(counts/counts.sum()) for g,counts in observations])
        fitted=least_squares(residual,[.35,.5,.5,.6],bounds=(DOMAIN_LOW,DOMAIN_HIGH),max_nfev=300)
        known=np.array(row['parameters']);candidates=[];solvers=[]
        for sign in (1.,-1.):
            for initial in (known,fitted.x,(DOMAIN_LOW+DOMAIN_HIGH)/2):
                fit=minimize(lambda x:sign*capture_snapshot(x,removal_factor=2.)[-1],initial,
                             method='SLSQP',bounds=list(zip(DOMAIN_LOW,DOMAIN_HIGH)),constraints={'type':'ineq','fun':constraints},
                             options={'maxiter':200,'ftol':1e-12})
                solvers.append({'success':bool(fit.success),'status':int(fit.status)})
                for fraction in (1.,.999,.99,.9,0.):
                    candidate=known+fraction*(fit.x-known);minimum=float(slack(candidate).min())
                    if minimum>=0 and np.all(candidate>=DOMAIN_LOW) and np.all(candidate<=DOMAIN_HIGH):
                        candidates.append({'parameters':candidate.tolist(),'minimum_calibration_slack':minimum,
                                           'target':capture_snapshot(candidate,removal_factor=2.).tolist()});break
        if not candidates:raise AssertionError('No feasible diagnostic witness')
        low=min(candidates,key=lambda x:x['target'][-1]);high=max(candidates,key=lambda x:x['target'][-1])
        gap=high['target'][-1]-low['target'][-1]
        rows.append({'parameters':row['parameters'],'repeat':row['repeat'],'lower_absence_witness':low,'upper_absence_witness':high,
                     'attained_absence_span':gap,'span_exceeds_precision':gap>.02,'solvers':solvers})
        print(row['parameters'],row['repeat'],'attained span',gap,'solver successes',sum(x['success'] for x in solvers),flush=True)
    report={'scope':'Exploratory numerical feasibility witnesses after outcomes. Same data; no new confirmation or globally optimal extrema.',
            'interpretation':'A span >.02 prevents a .02 envelope covering every parameter consistent with the chosen rectangular calibration constraints; not an impossibility for every statistical procedure.',
            'rows':rows}
    (ROOT/'results/occupancy_ambiguity_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
