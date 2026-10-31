"""Noisy band identification and held-out intervention predictions."""
import json
from pathlib import Path
import numpy as np
from clock_identification import identify_clock
from study_kinetic_transfer import kinetic_contacts
from study_load_probes import patterns
from counterfactual_bounds import symmetric_contacts
from study_identification import contact_array

ROOT=Path(__file__).resolve().parent


def main():
    rng=np.random.default_rng(20261116);rows=[]
    for m in (16,64):
        loads=patterns(m,'alternating',None)
        for repeat in range(3):
            l=rng.uniform(.2,.6,m-1);r=rng.uniform(.2,.6,m-1)
            for s in (.95,.8,.5,'mixture'):
                def generate(q,left,right):
                    if s=='mixture':return .5*(kinetic_contacts(q,left,right,.3,True)+kinetic_contacts(q,left,right,.95,True))
                    return kinetic_contacts(q,left,right,s,True)
                maps=[generate(q,l,r) for q in loads]
                second=np.stack([np.diag(d,2) for d in maps]);third=np.stack([np.diag(d,3) for d in maps])
                noise2=rng.normal(size=second.shape);noise3=rng.normal(size=third.shape)
                for sigma in (0.,.0001,.001,.01):
                    eta=sigma*np.sqrt(2*np.log(2*(second.size+third.size)/.05))
                    fit=identify_clock(loads,second+sigma*noise2,third+sigma*noise3,eta)
                    row={'size':m,'repeat':repeat,'true_survival':s,'noise_sd':sigma,'observed_scalars':second.size+third.size,'empty':fit['empty']}
                    if not fit['empty']:
                        slo,shi=fit['survival_interval'];sref=.625 if s=='mixture' else s
                        errors=[];static_errors=[];widths=[];covered=True;effective=sref*np.stack([l,r])
                        for cuts in ('none','left','both'):
                            raw=np.stack([l,r]).copy();point=fit['gate_estimate'].copy();low=fit['gate_low'].copy();high=fit['gate_high'].copy();static=effective.copy()
                            if cuts in ('left','both'):
                                for g in (raw,point,low,high,static):g[0,m//3]=0
                            if cuts=='both':
                                for g in (raw,point,low,high,static):g[1,2*m//3]=0
                            q=np.full(m,.5);truth=generate(q,*raw).sum(axis=1)
                            predicted=kinetic_contacts(q,*point,fit['survival'],True).sum(axis=1)
                            lower=kinetic_contacts(q,*low,slo,True).sum(axis=1) if slo>0 else np.r_[q,0]+np.r_[0,q]
                            upper=kinetic_contacts(q,*high,shi,True).sum(axis=1)
                            baseline=symmetric_contacts(contact_array(q,*static)).sum(axis=1)
                            errors.extend((predicted-truth).tolist());static_errors.extend((baseline-truth).tolist());widths.extend((upper-lower).tolist())
                            covered&=bool(np.all((lower-1e-10<=truth)&(truth<=upper+1e-10)))
                        row.update({'estimated_survival':fit['survival'],'survival_interval':fit['survival_interval'],
                                    'clock_absolute_error':abs(fit['survival']-sref),'clock_covered':bool(slo-1e-12<=sref<=shi+1e-12) if s!='mixture' else None,
                                    'heldout_all_covered':covered,'message_rmse':float(np.sqrt(np.mean(np.array(errors)**2))),
                                    'static_message_rmse':float(np.sqrt(np.mean(np.array(static_errors)**2))),'mean_prediction_width':float(np.mean(widths))})
                    rows.append(row)
    (ROOT/'results/clock_identification_study.json').write_text(json.dumps({'scope':'Exploratory simulated band recovery; mixture rows violate the assumed family.','rows':rows},indent=2)+'\n')
    for sigma in (0.,.0001,.001,.01):
        selected=[x for x in rows if x['noise_sd']==sigma and x['true_survival']!='mixture'];valid=[x for x in selected if not x['empty']]
        print(sigma,'nonempty',len(valid),'of',len(selected),'clock MAE',np.mean([x['clock_absolute_error'] for x in valid]),'RMSE',np.mean([x['message_rmse'] for x in valid]),'static',np.mean([x['static_message_rmse'] for x in valid]),'all coverage',sum(x['heldout_all_covered'] for x in valid),'width',np.mean([x['mean_prediction_width'] for x in valid]),'mixture rejected',sum(x['empty'] for x in rows if x['noise_sd']==sigma and x['true_survival']=='mixture'),flush=True)


if __name__=='__main__':main()
