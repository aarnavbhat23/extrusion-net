"""Small causal character-model diagnostic with fixed contiguous data splits."""
import argparse,copy,hashlib,json,math,time
from pathlib import Path
import numpy as np
import torch
from torch.nn import functional as F
from moment_text import MomentLanguageModel

SOURCE='https://raw.githubusercontent.com/karpathy/char-rnn/370cbcd448eb7daf32f21a6be560b70e0b33c4e3/data/tinyshakespeare/input.txt'
EXPECTED='86c4e6aa9db7c042ec79f339dcb96d42b0075e16b8fc2e86bf0ca57e2dc565ed'
ROOT=Path(__file__).resolve().parent


@torch.no_grad()
def evaluate(model,data):
    model.eval();loss=0.;count=0;starts=list(range(0,len(data)-65,64))
    for offset in range(0,len(starts),8):
        group=starts[offset:offset+8];x=torch.stack([data[i:i+64] for i in group]);y=torch.stack([data[i+1:i+65] for i in group])
        logits=model(x);loss+=float(F.cross_entropy(logits.flatten(0,1),y.flatten(),reduction='sum'));count+=y.numel()
    model.train();return loss/count,count


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--corpus',required=True);args=parser.parse_args()
    torch.set_num_threads(2);raw=Path(args.corpus).read_bytes();assert hashlib.sha256(raw).hexdigest()==EXPECTED
    text=raw.decode();first=int(.8*len(text));second=int(.9*len(text));vocabulary=sorted(set(text[:first]));mapping={c:i for i,c in enumerate(vocabulary)}
    unknown=len(vocabulary);encoded=torch.tensor([mapping.get(c,unknown) for c in text]);train=encoded[:first];val=encoded[first:second];test=encoded[second:]
    report={'scope':'Small character-model diagnostic, reset context 64; no future-token certificate or scaling claim.',
            'protocol_commit':'016e546','source_url':SOURCE,'source_sha256':EXPECTED,'split_char_indices':[0,first,second,len(text)],
            'training_vocabulary':vocabulary,'unknown_token':unknown,'validation_unknown_characters':int((val==unknown).sum()),'test_unknown_characters':int((test==unknown).sum()),'rows':[]}
    for seed in (0,1,2):
        for kind in ('atom','mean','static','attention'):
            torch.manual_seed(seed);model=MomentLanguageModel(kind,len(vocabulary)+1);opt=torch.optim.AdamW(model.parameters(),lr=.002,weight_decay=.01)
            rng=np.random.default_rng(9800+seed);best=float('inf');state=None;curve=[];train_time=0.;all_start=time.perf_counter()
            for step in range(1,601):
                start=time.perf_counter();ids=rng.integers(0,len(train)-65,size=8)
                x=torch.stack([train[i:i+64] for i in ids]);y=torch.stack([train[i+1:i+65] for i in ids])
                logits=model(x);loss=F.cross_entropy(logits.flatten(0,1),y.flatten())
                opt.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1.);opt.step();train_time+=time.perf_counter()-start
                if step%150==0:
                    score,count=evaluate(model,val);curve.append({'step':step,'training_loss':float(loss.detach()),'validation_nats':score})
                    if score<best:best=score;state=copy.deepcopy(model.state_dict())
            model.load_state_dict(state);test_loss,count=evaluate(model,test)
            row={'seed':seed,'kind':kind,'params':sum(p.numel() for p in model.parameters()),'train_seconds_excluding_validation':train_time,
                 'total_seconds':time.perf_counter()-all_start,'train_tokens':600*8*64,'training_tokens_per_second':600*8*64/train_time,
                 'best_validation_nats':best,'test_nats':test_loss,'test_bits_per_character':test_loss/math.log(2),'test_target_characters':count,'curve':curve}
            if kind in ('atom','mean'):
                bits=model.value_dim;row['content_memory_state_scalars']=32*(bits+1)**2
                row['content_projection_buffer_scalars']=model.projection.numel()
            report['rows'].append(row);print(seed,kind,'test nats',round(test_loss,4),'tokens/s',round(row['training_tokens_per_second']),flush=True)
            (ROOT/'results/moment_text_study.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
