"""Post-screen resource check; saved entries are not peak allocated memory."""
import json,time
from pathlib import Path
import numpy as np
import torch
from causal_transport import CausalModel


def main():
    torch.set_num_threads(2);rows=[]
    for n in (64,256,1024,4096):
        for variant in ('static','shared','independent','attention'):
            torch.manual_seed(91);model=CausalModel(variant,25,8);tokens=torch.randint(25,(2,n));times=[]
            saved=[];parameters={p.untyped_storage().data_ptr() for p in model.parameters()}
            def pack(t):
                if t.untyped_storage().data_ptr() not in parameters:saved.append(t.numel())
                return t
            with torch.autograd.graph.saved_tensors_hooks(pack,lambda t:t):loss=model(tokens).square().mean()
            loss.backward();model.zero_grad(set_to_none=True)
            for repeat in range(7):
                start=time.perf_counter();loss=model(tokens).square().mean();loss.backward();model.zero_grad(set_to_none=True)
                if repeat>=2:times.append(time.perf_counter()-start)
            row={'length':n,'variant':variant,'batch':2,'median_forward_backward_seconds':float(np.median(times)),
                 'tokens_per_second':2*n/float(np.median(times)),'saved_nonparameter_tensor_entries':sum(saved)}
            rows.append(row);print(n,variant,round(row['tokens_per_second']),row['saved_nonparameter_tensor_entries'],flush=True)
    Path('results/causal_transport_resources.json').write_text(json.dumps({'scope':'Untrained CPU resource screen, two threads. Saved entries can count aliases and are not peak memory.','rows':rows},indent=2)+'\n')


if __name__=='__main__':main()
