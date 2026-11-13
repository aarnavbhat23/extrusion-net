"""Matched causal architecture screen; no long-length checkpoint selection."""
import json,time,copy
from pathlib import Path
import numpy as np
import torch
from torch.nn import functional as F
import task_dyck,task_recall
from causal_transport import CausalModel

ROOT=Path(__file__).resolve().parent


def tensors(split):return {n:tuple(torch.from_numpy(a) for a in pair) for n,pair in split.items()}


@torch.no_grad()
def accuracy(model,data):
    correct=total=0;model.eval()
    for x,y in data.values():
        for i in range(0,len(x),16):
            target=y[i:i+16];keep=target!=-100;pred=model(x[i:i+16]).argmax(-1)
            correct+=int(((pred==target)&keep).sum());total+=int(keep.sum())
    model.train();return correct/total


def main():
    torch.set_num_threads(2);rows=[]
    for name,data in [('recall',task_recall),('dyck',task_dyck)]:
        train=tensors(data.make_split(200,800,[24,32,48,64]));val=tensors(data.make_split(201,192,[24,32,48,64]))
        tests={n:tensors(data.make_split(202+n,128,[n])) for n in (64,128,256)}
        for seed in (0,1):
            for variant in ('static','shared','independent','attention'):
                torch.manual_seed(seed);rng=np.random.default_rng(7100+seed)
                model=CausalModel(variant,len(data.TOKENS),data.N_CLASSES);opt=torch.optim.AdamW(model.parameters(),lr=.002,weight_decay=.01)
                lengths=sorted(train);best=-1.;state=None;curve=[];start=time.perf_counter()
                for step in range(1,401):
                    n=int(rng.choice(lengths));x,y=train[n];ids=rng.choice(len(x),min(24,len(x)),replace=False)
                    logits=model(x[ids]);loss=F.cross_entropy(logits.flatten(0,1),y[ids].flatten(),ignore_index=-100)
                    opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);opt.step()
                    if step%100==0:
                        score=accuracy(model,val);curve.append({'step':step,'loss':float(loss.detach()),'validation':score})
                        if score>best:best=score;state=copy.deepcopy(model.state_dict())
                elapsed=time.perf_counter()-start;model.load_state_dict(state)
                row={'task':name,'seed':seed,'variant':variant,'params':sum(p.numel() for p in model.parameters()),'train_seconds':elapsed,'best_validation':best,'curve':curve,
                     'test_accuracy':{str(n):accuracy(model,d) for n,d in tests.items()}}
                rows.append(row);print(name,seed,variant,row['test_accuracy'],'seconds',round(elapsed,1),flush=True)
                (ROOT/'results/causal_clock_study.json').write_text(json.dumps({'scope':'Exploratory causal architecture screen, two seeds, synthetic tasks.','rows':rows},indent=2)+'\n')


if __name__=='__main__':main()
