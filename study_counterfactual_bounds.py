"""Check exact matrix segments and signed production-mixer output bounds."""
import json
from pathlib import Path

import numpy as np

from counterfactual_bounds import (cut_endpoints, initial_for_mix, mixer_outputs,
                                  sharp_output_bounds, symmetric_contacts)
from study_identification import contact_array, factor_at_initial, observable_bands

ROOT=Path(__file__).resolve().parent


def main():
    rng=np.random.default_rng(20261103);rows=[];witness=None
    for size in (8,16,32):
        for case in range(20):
            q=rng.uniform(.4,.9,size)
            left,right=rng.uniform(.4,.7,(2,size-1))
            baseline=contact_array(q,left,right)
            site=int(rng.integers(size-1));multiplier=float(rng.choice([0.,.3,.7]))
            endpoints,interval,prefix=cut_endpoints(baseline,site,multiplier)
            values=rng.normal(size=(size+1,4))
            weight=rng.normal(size=(6,5));weight[:,-1]*=4
            bias=rng.normal(size=6)
            bounds=sharp_output_bounds(endpoints,values,weight,bias)
            endpoint_outputs=np.asarray([mixer_outputs(c,values,weight,bias) for c in endpoints])
            endpoint_low=endpoint_outputs.min(0);endpoint_high=endpoint_outputs.max(0)
            q,sums,products=observable_bands(baseline)
            segment_error=0.;enclosure_error=0.;signed_enclosure_error=0.
            for theta in np.linspace(0,1,101):
                initial=initial_for_mix(theta,interval,prefix)
                l,r=factor_at_initial(q,sums,products,initial);r[site]*=multiplier
                actual=symmetric_contacts(contact_array(q,l,r))
                segment=(1-theta)*endpoints[0]+theta*endpoints[1]
                segment_error=max(segment_error,float(np.max(np.abs(actual-segment))))
                output=mixer_outputs(actual,values,weight,bias)
                enclosure_error=max(enclosure_error,float(np.max(bounds['lower']-output)),float(np.max(output-bounds['upper'])))
                signed=actual@values
                signed_endpoints=endpoints@values
                signed_enclosure_error=max(signed_enclosure_error,float(np.max(signed_endpoints.min(0)-signed)),
                                           float(np.max(signed-signed_endpoints.max(0))))
            # Reconstruct physical gates at every claimed extremizer and evaluate
            # the original contact sum, independently of the scalar optimizer.
            attainment_error=0.
            for side in ('lower','upper'):
                for position in range(size+1):
                    for channel in range(6):
                        theta=bounds[side+'_theta'][position,channel]
                        initial=initial_for_mix(theta,interval,prefix)
                        l,r=factor_at_initial(q,sums,products,initial);r[site]*=multiplier
                        output=mixer_outputs(symmetric_contacts(contact_array(q,l,r)),values,weight,bias)
                        attainment_error=max(attainment_error,abs(float(output[position,channel]-bounds[side][position,channel])))
            missed=np.maximum(endpoint_low-bounds['lower'],bounds['upper']-endpoint_high)
            row={'size':size,'case':case,'site':site,'multiplier':multiplier,
                 'matrix_segment_max_error':segment_error,'signed_message_enclosure_error':signed_enclosure_error,
                 'full_output_enclosure_error':enclosure_error,'extremum_attainment_error':attainment_error,
                 'endpoint_only_missed_outputs':int(np.sum(missed>1e-8)),
                 'endpoint_only_max_miss':float(missed.max()),'output_coordinates':int(missed.size)}
            rows.append(row)
            if row['endpoint_only_max_miss']>1e-5 and (witness is None or row['endpoint_only_max_miss']>witness['miss']):
                position,channel=np.unravel_index(np.argmax(missed),missed.shape)
                witness={'size':size,'case':case,'site':site,'multiplier':multiplier,
                         'q':q.tolist(),'left':left.tolist(),'right':right.tolist(),'values':values.tolist(),
                         'weight':weight.tolist(),'bias':bias.tolist(),'position':int(position),'channel':int(channel),
                         'endpoint_interval':[float(endpoint_low[position,channel]),float(endpoint_high[position,channel])],
                         'sharp_interval':[float(bounds['lower'][position,channel]),float(bounds['upper'][position,channel])],
                         'extremizer_theta':[float(bounds['lower_theta'][position,channel]),float(bounds['upper_theta'][position,channel])],
                         'miss':float(missed.max())}
        print('completed',size,flush=True)
    report={'seed':20261103,'scope':'60 exact baseline mechanisms; fixed one-gate intervention; single unsaturated head and affine output. Not a certificate for a full nonlinear network or noisy baseline.',
            'rows':rows,'interior_extremum_witness':witness,
            'summary':{'cases':len(rows),'matrix_segment_max_error':max(r['matrix_segment_max_error'] for r in rows),
                       'output_enclosure_max_error':max(r['full_output_enclosure_error'] for r in rows),
                       'attainment_max_error':max(r['extremum_attainment_error'] for r in rows),
                       'endpoint_only_missed_outputs':sum(r['endpoint_only_missed_outputs'] for r in rows),
                       'output_coordinates':sum(r['output_coordinates'] for r in rows)}}
    (ROOT/'results/counterfactual_output_bounds.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['summary'],indent=2))


if __name__=='__main__':main()
