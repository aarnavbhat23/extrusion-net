"""Frozen unseen-key and overload screen, supplied-binding operators only."""
import copy,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from torch.nn import functional as F
from moment_memory import RetrievalOperator

ROOT=Path(__file__).resolve().parent
KINDS=('atom','mean','mean_same_width','delta','rls','attention','attention_fifo')


def pool(seed,count):
    rng=np.random.default_rng(seed);x=rng.normal(size=(count,16)).astype('float32');x/=np.linalg.norm(x,axis=1,keepdims=True)
    return torch.from_numpy(x)


def episode(keys,rng,n,batch):
    ids=np.stack([rng.choice(len(keys),n,replace=False) for _ in range(batch)])
    return keys[ids],torch.from_numpy(rng.integers(8,size=(batch,n)))


@torch.no_grad()
def evaluate(model,keys,labels):
    model=model.double();raw=keys.double();start=time.perf_counter();logits,compatible=model(raw,raw,labels)
    predictions=logits.argmax(-1);truth=labels;result={'accuracy':float((predictions==truth).double().mean()),'evaluation_seconds':time.perf_counter()-start}
    if compatible is not None:
        certified=compatible.sum(-1)==1
        result.update({'certified_coverage':float(certified.double().mean()),'false_certificates':int((certified&(predictions!=truth)).sum()),
                       'accuracy_counting_abstention_as_failure':float((certified&(predictions==truth)).double().mean())})
    return result


def main():
    torch.set_num_threads(2);train=pool(9200,1024);validation=pool(9201,256);test=pool(9202,1024)
    pools=[set(map(bytes,x.numpy())) for x in (train,validation,test)]
    assert all(not pools[i].intersection(pools[j]) for i in range(3) for j in range(i+1,3))
    report={'scope':'Supplied-binding operators, not learned token-stream binding or a language model.','protocol_commit':'016e546','rows':[],
            'pool_sha256':{name:hashlib.sha256(x.numpy().tobytes()).hexdigest() for name,x in [('train',train),('validation',validation),('test',test)]}}
    vrng=np.random.default_rng(9203);val_keys,val_labels=episode(validation,vrng,8,64)
    fixed_tests={n:episode(test,np.random.default_rng(9300+n),n,64) for n in (8,16,32,64)}
    for seed in (0,1,2):
        for kind in KINDS:
            torch.manual_seed(seed);model=RetrievalOperator(kind);opt=torch.optim.Adam(model.parameters(),lr=.002);rng=np.random.default_rng(9400+seed)
            best=-1.;state=None;curve=[];start=time.perf_counter()
            for step in range(1,601):
                keys,labels=episode(train,rng,8,16);logits,_=model(keys,keys,labels)
                loss=F.cross_entropy(logits.flatten(0,1),labels.flatten())
                opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);opt.step()
                if step%100==0:
                    with torch.no_grad():score=float((model(val_keys,val_keys,val_labels)[0].argmax(-1)==val_labels).float().mean())
                    curve.append({'step':step,'loss':float(loss.detach()),'validation_accuracy':score})
                    if score>best:best=score;state=copy.deepcopy(model.state_dict())
            elapsed=time.perf_counter()-start;model.load_state_dict(state)
            tests={str(n):evaluate(model,*data) for n,data in fixed_tests.items()}
            row={'seed':seed,'kind':kind,'params':sum(p.numel() for p in model.parameters()),'projection_buffer_scalars':model.projection.numel(),
                 'state_scalars':{str(n):model.state_scalars(n) for n in fixed_tests},'train_seconds':elapsed,'validation_accuracy':best,'curve':curve,'tests':tests}
            if kind=='atom':
                keys,labels=fixed_tests[16];queries=pool(9500,1024)[:len(keys)].double()[:,None,:]
                with torch.no_grad():
                    _,compatible=model(keys.double(),queries,labels)
                    row['absent_query_singleton_rate']=float((compatible.sum(-1)==1).double().mean())
                    duplicated=keys.double().clone();duplicated[:,1]=duplicated[:,0];duplicate_labels=labels.clone();duplicate_labels[:,1]=(duplicate_labels[:,0]+1)%8
                    _,compatible=model(duplicated,duplicated[:,:1],duplicate_labels)
                    row['conflicting_key_singleton_rate']=float((compatible.sum(-1)==1).double().mean())
                torch.save({'kind':kind,'state_dict':model.float().state_dict()},ROOT/'results'/f'moment_operator_seed{seed}.pt')
            report['rows'].append(row)
            (ROOT/'results/moment_retrieval_study.json').write_text(json.dumps(report,indent=2)+'\n')
            print(seed,kind,[(n,round(v['accuracy'],4),round(v.get('certified_coverage',0),4)) for n,v in tests.items()],'seconds',round(elapsed,1),flush=True)


if __name__=='__main__':main()
