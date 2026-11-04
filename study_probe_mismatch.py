"""An explicit source-invariance alternative invisible to calibration probes."""
import json
from pathlib import Path
import numpy as np
from colored_recovery import choose_colors,colored_response,recover,scan_messages

ROOT=Path(__file__).resolve().parent


def positive_message_bounds(q,estimate,size,radius,rho):
    low=np.clip(estimate-radius,0,rho);high=np.clip(estimate+radius,0,rho)
    value=np.ones((size+1,1))
    return (scan_messages(q,low[:size-1],low[size-1:],value).ravel(),
            scan_messages(q,high[:size-1],high[size-1:],value).ravel())


def main():
    rng=np.random.default_rng(20261112);rows=[];stress=[];size=128;sigma=.001
    for rho in (.2,.5,.8):
        colors,rate=choose_colors(rho,size)
        for repeat in range(10):
            truth=rng.uniform(.1*rho,.7*rho,2*(size-1));clean=colored_response(truth,size,colors)
            observed=clean+rng.normal(0,sigma,clean.shape);estimate,fit=recover(observed,size,rho,colors,rate)
            eta=sigma*np.sqrt(2*np.log(2*2*(size-1)/.025))
            radius=fit['fixed_point_error_bound']+eta/(.8*(1-rate))
            q=np.full(size,.5);low,high=positive_message_bounds(q,estimate,size,radius,rho)
            hold_noise=rng.normal(0,sigma,size+1);hold_radius=sigma*np.sqrt(2*np.log(2*(size+1)/.025))
            for gamma in (0.,.005,.02,.05):
                # h(q)=0 on both calibration loads; h(q)=1 on this held-out load.
                source_factor=float(np.mean(1-((q-.5)/.4)**2))
                shifted=np.clip(truth+gamma*source_factor,0,rho)
                actual=scan_messages(q,shifted[:size-1],shifted[size-1:],np.ones((size+1,1))).ravel()
                response=actual+hold_noise
                rejected=bool(np.any((response<low-hold_radius)|(response>high+hold_radius)))
                rows.append({'rho':rho,'repeat':repeat,'gate_shift':gamma,'calibration_gate_radius':radius,
                             'actual_calibration_gate_max_error':float(np.max(abs(estimate-truth))),
                             'rejected':rejected,'heldout_noise_radius':hold_radius,
                             'observed_response_min':float(response.min()),'observed_response_max':float(response.max()),
                             'mean_prediction_interval_width':float(np.mean(high-low+2*hold_radius))})
    for size in (128,512):
        for rho in (.8,.95):
            colors,rate=choose_colors(rho,size);truth=np.full(2*(size-1),.999*rho)
            clean=colored_response(truth,size,colors);noise=rng.normal(size=clean.shape)
            for sigma in (0.,.001):
                estimate,fit=recover(clean+sigma*noise,size,rho,colors,rate)
                stress.append({'size':size,'rho':rho,'noise_sd':sigma,'colors':colors,'contraction_bound':rate,
                               'gate_rmse':float(np.sqrt(np.mean((estimate-truth)**2))),
                               'gate_max_error':float(np.max(abs(estimate-truth))),**fit})
    report={'scope':'Adversarial source-dependent alternative and near-transparent stress, not biological validation.',
            'seed':20261112,'invariant_null_joint_alpha':.05,'mismatch_rows':rows,'near_transparent_stress':stress}
    (ROOT/'results/probe_mismatch_study.json').write_text(json.dumps(report,indent=2)+'\n')
    for rho in (.2,.5,.8):
        print(rho,[(gamma,sum(x['rejected'] for x in rows if x['rho']==rho and x['gate_shift']==gamma)) for gamma in (0.,.005,.02,.05)],flush=True)


if __name__=='__main__':main()
