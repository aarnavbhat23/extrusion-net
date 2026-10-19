"""Numerical scope check for the experimental polynomial contact scan."""
import json
import time
from pathlib import Path

import torch

from models import extrusion_log_contacts
from polynomial_scan import contact_argument_bound,saturated_polynomial_scan,taylor_entry_error


def main():
    torch.set_num_threads(1)
    rows=[]
    for n in (32,256):
        generator=torch.Generator().manual_seed(20261018+n)
        load=torch.rand(n,generator=generator,dtype=torch.float64)*0.2+0.2
        pl=torch.rand(n,generator=generator,dtype=torch.float64)*0.3+0.5
        pr=torch.rand(n,generator=generator,dtype=torch.float64)*0.3+0.5
        values=torch.randn(n,3,generator=generator,dtype=torch.float64)
        beta=0.8
        start=time.perf_counter()
        contact=extrusion_log_contacts(load.log(),pl.log(),pr.log()).exp()
        reference=-torch.expm1(-beta*contact)@values
        dense_seconds=time.perf_counter()-start
        bound=contact_argument_bound(load,pl,pr,beta)
        for degree in (1,2,4,6,8):
            start=time.perf_counter()
            actual=saturated_polynomial_scan(load,pl,pr,values,beta,degree)
            seconds=time.perf_counter()-start
            message_bound=taylor_entry_error(bound,degree)*values.abs().sum(0).max()
            error=(actual-reference).abs().max()
            assert error<=message_bound+1e-10
            rows.append({'length':n,'degree':degree,'max_argument':float((beta*contact).max()),
                         'argument_upper_bound':float(bound),'max_message_abs_error':float(error),
                         'message_error_bound_exact_arithmetic':float(message_bound),
                         'scan_seconds':seconds,'dense_seconds':dense_seconds})
    n=32
    ones=torch.ones(n,dtype=torch.float64)
    values=torch.ones(n,1,dtype=torch.float64)
    beta=4.0
    contact=extrusion_log_contacts(ones.log(),ones.log(),ones.log()).exp()
    approximation=saturated_polynomial_scan(ones,ones,ones,values,beta,8)
    exact=-torch.expm1(-beta*contact)@values
    failure={'length':n,'degree':8,'beta':beta,
             'argument_upper_bound':float(contact_argument_bound(ones,ones,ones,beta)),
             'max_message_abs_error':float((approximation-exact).abs().max()),
             'scope':'Deliberately outside the useful bounded regime; not a usable saturated mixer.'}
    result={'seed':20261018,'dtype':'float64','rows':rows,'large_argument_failure':failure,
            'scope':'Operator consistency and approximation only. No new trained model or task improvement.'}
    Path('results/polynomial_scan.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
