"""Causal local-context memory diagnostic; no certificate for future tokens."""
import math
import torch
from torch import nn
from torch.nn import functional as F
from causal_transport import CausalModel,delay
from moment_memory import positive_features,value_codes


class MomentLanguageModel(nn.Module):
    def __init__(self,kind,vocab,d=32):
        super().__init__();self.kind=kind;self.base=CausalModel('attention' if kind=='attention' else 'static',vocab,vocab,d=d)
        if kind in ('atom','mean'):
            bits=math.ceil(math.log2(vocab));self.value_dim=bits
            self.address=nn.Linear(d,16,bias=False);self.memory_out=nn.Linear(bits,d,bias=False)
            self.gate=nn.Linear(d,1);nn.init.constant_(self.gate.bias,-2.)
            self.log_temperature=nn.Parameter(torch.tensor(0.))
            width=32 if kind=='atom' else 32*(bits+1)
            generator=torch.Generator().manual_seed(9700)
            self.register_buffer('projection',torch.randn(width,16,generator=generator)/4)
            self.register_buffer('codes',value_codes(bits)[:vocab])
    def memory_read(self,local,tokens):
        features=positive_features(torch.tanh(self.address(local)),self.projection)
        keys=delay(features,1);values=self.codes[tokens]
        if self.kind=='atom':
            augmented=torch.cat([torch.ones_like(values[...,:1]),values],dim=-1)
            increments=keys[...,None,None]*augmented[:,:,None,:,None]*augmented[:,:,None,None,:]
            prefix=increments.cumsum(dim=1)
            gram=(features[...,None,None]*prefix).sum(dim=2)
            gram=gram+1e-3*torch.eye(self.value_dim+1,device=local.device,dtype=local.dtype)
            catalogue=torch.cat([torch.ones_like(self.codes[:,:1]),self.codes],dim=-1)
            solved=torch.linalg.solve(gram,catalogue.T.expand(*gram.shape[:-2],*catalogue.T.shape))
            scores=torch.einsum('ci,bnic->bnc',catalogue,solved)
            predicted=torch.softmax(-scores*self.log_temperature.exp().clamp(max=20),dim=-1)@self.codes
        else:
            augmented=torch.cat([torch.ones_like(values[...,:1]),values],dim=-1)
            prefix=(keys[...,None]*augmented[:,:,None,:]).cumsum(dim=1)
            mixed=(features[...,None]*prefix).sum(dim=2)
            predicted=mixed[...,1:]/mixed[...,:1].clamp_min(1e-8)
        return self.memory_out(predicted)*torch.sigmoid(self.gate(local))
    def forward(self,tokens):
        base=self.base;h=base.embed(tokens)
        h=h+F.gelu(base.conv(F.pad(base.norms[0](h).transpose(1,2),(4,0))).transpose(1,2))
        if self.kind in ('atom','mean'):h=h+self.memory_read(base.norms[1](h),tokens)
        h=h+base.mixer(base.norms[1](h));h=h+base.mlp(base.norms[2](h))
        return base.head(base.norms[3](h))
