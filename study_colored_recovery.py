"""Test the response-budget and accuracy tradeoff of spaced probing."""
import json
from pathlib import Path
import numpy as np
from colored_recovery import choose_colors,colored_response,recover,scan_messages

ROOT=Path(__file__).resolve().parent


def main():
    rng=np.random.default_rng(20261111);rows=[];alpha=.05
    for size in (32,128,512):
        for rho in (.2,.5,.8,.95):
            colors,rate=choose_colors(rho,size)
            for repeat in range(5):
                truth=rng.uniform(0,rho,2*(size-1));clean=colored_response(truth,size,colors)
                # Only these selected response coordinates enter the recovery step.
                noise=rng.normal(size=clean.shape);indices=np.arange(size-1)
                used_noise=noise[:,indices,(indices+2)%colors]
                q=rng.uniform(.1,.9,size);v=rng.choice([-1.,1.],size=(size+1,1));cut=size//2
                true_cut=truth.copy();true_cut[size-1+cut]=0
                target=scan_messages(q,true_cut[:size-1],true_cut[size-1:],v)
                for sigma in (0.,.001):
                    observed=clean+sigma*noise;estimate,fit=recover(observed,size,rho,colors,rate)
                    eta=float(sigma*np.sqrt(2*np.log(2*2*(size-1)/alpha)))
                    bound=fit['fixed_point_error_bound']+eta/(.8*(1-rate))
                    exact_bound=fit['fixed_point_error_bound']+float(sigma*np.max(abs(used_noise)))/(.8*(1-rate))
                    error=estimate-truth;estimated_cut=estimate.copy();estimated_cut[size-1+cut]=0
                    predicted=scan_messages(q,estimated_cut[:size-1],estimated_cut[size-1:],v)
                    rows.append({'size':size,'rho':rho,'repeat':repeat,'noise_sd':sigma,'colors':colors,
                                 'vector_responses':2*colors,'response_scalars':2*colors*(size+1),
                                 'selected_response_scalars':2*(size-1),'total_value_energy_both_loads':2*(size+1),
                                 'contraction_bound':rate,'gate_rmse':float(np.sqrt(np.mean(error**2))),
                                 'gate_max_error':float(np.max(abs(error))),'noise_bound_95':eta,
                                 'error_bound_95':bound,'actual_noise_error_bound':exact_bound,
                                 'bound_covers':bool(np.max(abs(error))<=bound+1e-10),
                                 'heldout_cut_message_rmse':float(np.sqrt(np.mean((predicted-target)**2))),**fit})
                print(size,rho,repeat,'responses',2*colors,'rate',round(rate,4),'iterations',fit['iterations'],flush=True)
    report={'scope':'Global recovery with spaced probing under a known uniform gate bound and exact known loading.',
            'seed':20261111,'joint_gate_error_alpha':alpha,'rows':rows}
    (ROOT/'results/colored_recovery_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
