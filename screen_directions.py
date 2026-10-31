"""Small independent structural screens, not architecture performance claims."""
import json
from pathlib import Path

import numpy as np
import torch

from models import extrusion_log_contacts

ROOT = Path(__file__).resolve().parent


def contacts(q, left, right):
    """Reduced upper-triangular contact array C[i,j] = D[i,j+1]."""
    size = len(q)
    eye = torch.eye(size,dtype=q.dtype,device=q.device)
    return torch.linalg.solve(eye-torch.diag(left,1), torch.diag(q)) @ torch.linalg.inv(eye-torch.diag(right,1))


def equivalent(q, left, right, first_right):
    """Factor the same inverse band operator using a different initial gate."""
    sums = left/q[:-1] + right/q[1:]
    products = right[:-1]*left[1:]/q[1:-1]
    new_r = [torch.as_tensor(first_right,dtype=q.dtype,device=q.device)]
    new_l = [q[0]*(sums[0]-new_r[0]/q[1])]
    for k in range(1,len(right)):
        new_l.append(products[k-1]*q[k]/new_r[-1])
        new_r.append(q[k+1]*(sums[k]-new_l[-1]/q[k]))
    return torch.stack(new_l),torch.stack(new_r)


def svd_rank(jac):
    values = torch.linalg.svdvals(jac.reshape(-1,jac.shape[-1]))
    return {'rank':int((values > values[0]*1e-9).sum()),
            'parameters':int(jac.shape[-1]),'singular_values':values.tolist()}


def main():
    torch.set_default_dtype(torch.float64)
    torch.set_num_threads(1)
    generator = torch.Generator().manual_seed(20261020)
    free=[]; shared=[]; saturation=[]
    max_direct_error=0.
    for size in (7,15,31):
        q=.4+.5*torch.rand(size,generator=generator)
        left=.4+.3*torch.rand(size-1,generator=generator)
        right=.4+.3*torch.rand(size-1,generator=generator)
        z=torch.cat([left,right])
        raw=lambda z:contacts(q,z[:size-1],z[size-1:])
        c=raw(z)
        # Independent production sum uses positions 1..n-2 as physical gates.
        load=torch.cat([torch.ones(1),q])
        pl=torch.cat([torch.ones(1),left,torch.ones(1)])
        pr=torch.cat([torch.ones(1),right,torch.ones(1)])
        production=extrusion_log_contacts(load.log(),pl.log(),pr.log()).exp()
        check=torch.triu(production[:-1,1:])
        max_direct_error=max(max_direct_error,float((check-c).abs().max()))
        rank=svd_rank(torch.autograd.functional.jacobian(raw,z))
        saturated=svd_rank(torch.autograd.functional.jacobian(lambda z:-torch.expm1(-4*raw(z)),z))
        candidates=[]
        for r0 in np.linspace(.001,.999,2001):
            l,r=equivalent(q,left,right,r0)
            if bool(((l>0)&(l<1)&(r>0)&(r<1)).all()):
                error=float((contacts(q,l,r)-c).abs().max())
                cut_l,cut_r=l.clone(),r.clone();cut_r[0]=0
                candidates.append((float(r0),float(contacts(q,cut_l,cut_r)[0,1]),error))
        reference_r=right.clone();reference_r[0]=0
        intervened=contacts(q,left,reference_r)
        recovered_r0=(c[0,1]-intervened[0,1])/q[0]
        recovered_l,recovered_r=equivalent(q,left,right,recovered_r0)
        free.append({'positions':size+1,'raw':rank,'saturated':saturated,
                     'admissible_grid_points':len(candidates),
                     'admissible_first_right_range':[candidates[0][0],candidates[-1][0]],
                     'intervention_contact_range':[min(x[1] for x in candidates),max(x[1] for x in candidates)],
                     'max_observational_difference':max(x[2] for x in candidates),
                     'single_intervention_gate_recovery_error':float(torch.cat([recovered_l-left,recovered_r-right]).abs().max())})
    size=31
    q=.4+.5*torch.rand(size,generator=generator)
    left=.65+.3*torch.rand(size-1,generator=generator)
    right=.65+.3*torch.rand(size-1,generator=generator)
    c=contacts(q,left,right)
    # A separated block corresponds to disjoint physical endpoint sets.
    block=c[:12,17:]
    for beta in (.01,.1,1.,4.,16.):
        w=-torch.expm1(-beta*block)
        singular=torch.linalg.svdvals(w)
        saturation.append({'beta':beta,'raw_block_rank':int(torch.linalg.matrix_rank(block)),
                           'saturated_block_rank_relative_1e8':int((singular>singular[0]*1e-8).sum()),
                           'singular_values':singular.tolist()})
    for contexts in (1,2,4):
        size=15;features=4
        x=torch.randn(contexts,size-1,features,generator=generator)
        x[:,:,0]=1
        q=.4+.5*torch.rand(contexts,size,generator=generator)
        parameter=.3*torch.randn(2*features,generator=generator)
        def output(theta):
            gates=torch.sigmoid(x@theta.reshape(features,2))
            return torch.stack([contacts(q[k],gates[k,:,0],gates[k,:,1]) for k in range(contexts)])
        shared.append({'contexts':contexts,**svd_rank(torch.autograd.functional.jacobian(output,parameter))})
    depth=[]
    for size in (8,16):
        for layers in (1,2,3,4):
            theta=.3+.5*torch.rand(layers*(size-1),generator=generator)
            def cascade(parameters):
                output=torch.eye(size)
                for row in parameters.reshape(layers,size-1):
                    output=output@torch.linalg.inv(torch.eye(size)-torch.diag(row,1))
                return output
            band_coordinates=sum(size-k for k in range(1,layers+1))
            check=svd_rank(torch.autograd.functional.jacobian(cascade,theta))
            depth.append({'positions':size,'cascade_depth':layers,
                          'inverse_band_coordinates':band_coordinates,
                          'dimension_nullity_lower_bound':layers*(layers-1)//2,**check})
    report={'seed':20261020,'direct_production_max_error':max_direct_error,
            'free_gate_identification':free,'saturation_rank':saturation,'shared_logistic_controller':shared,
            'positive_cascade_depth':depth,
            'limitations':['Known positive loads and full noiseless contact observations.',
                          'Local Jacobian evidence is not global identifiability.',
                          'A free-gate equivalence family need not be realizable by a shared neural controller.',
                          'Entrywise saturation is invertible at known scale in exact arithmetic.']}
    (ROOT/'results/direction_screen.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Wrote structural screens to results/direction_screen.json')


if __name__=='__main__':
    main()
