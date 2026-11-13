"""Positive-feature value moments and a conditional atom-compatibility read."""
import itertools,math
import torch
from torch import nn
from torch.nn import functional as F


def positive_features(x,projection):
    values=F.relu(x@projection.T).pow(4)+1e-6
    return F.normalize(values,dim=-1)


def value_codes(bits=3,dtype=torch.float32):
    return torch.tensor(list(itertools.product((-1.,1.),repeat=bits)),dtype=dtype)/math.sqrt(bits)


def moment_state(keys,values):
    augmented=torch.cat([torch.ones_like(values[...,:1]),values],dim=-1)
    return torch.einsum('bnd,bni,bnj->bdij',keys,augmented,augmented)


def atom_scores(state,queries,codes,ridge=1e-3):
    gram=torch.einsum('bqd,bdij->bqij',queries,state)
    gram=gram+ridge*torch.eye(gram.shape[-1],device=gram.device,dtype=gram.dtype)
    augmented=torch.cat([torch.ones_like(codes[...,:1]),codes],dim=-1)
    solved=torch.linalg.solve(gram,augmented.T.expand(*gram.shape[:-2],*augmented.T.shape))
    scores=torch.einsum('ci,bqic->bqc',augmented,solved)
    return scores


class RetrievalOperator(nn.Module):
    widths={'atom':32,'mean':128,'mean_same_width':32,'delta':170,'rls':21,'attention':16,'attention_fifo':16}
    def __init__(self,kind,input_dim=16):
        super().__init__();self.kind=kind;self.width=self.widths[kind]
        self.encoder=nn.Linear(input_dim,input_dim,bias=False);self.log_temperature=nn.Parameter(torch.tensor(2.))
        with torch.no_grad():self.encoder.weight.copy_(torch.eye(input_dim)+.02*torch.randn(input_dim,input_dim))
        generator=torch.Generator().manual_seed(9101)
        self.register_buffer('projection',torch.randn(self.width,input_dim,generator=generator)/math.sqrt(input_dim))
        self.register_buffer('codes',value_codes())
    def features(self,raw):
        encoded=F.normalize(self.encoder(raw),dim=-1)
        return encoded if self.kind.startswith('attention') else positive_features(encoded,self.projection)
    def forward(self,raw_keys,raw_queries,labels):
        keys=self.features(raw_keys);queries=self.features(raw_queries);values=self.codes[labels]
        b,n,d=keys.shape;compatibility=None
        if self.kind=='atom':
            scores=atom_scores(moment_state(keys,values),queries,self.codes);logits=-scores
            compatibility=scores<=1+1e-6
        elif self.kind in ('mean','mean_same_width'):
            state=torch.einsum('bnd,bnv->bdv',keys,values)
            mass=torch.einsum('bqd,bd->bq',queries,keys.sum(dim=1))
            predicted=torch.einsum('bqd,bdv->bqv',queries,state)/mass.clamp_min(1e-8)[...,None]
            logits=-(predicted[...,None,:]-self.codes).square().sum(dim=-1)
        elif self.kind=='delta':
            state=torch.zeros(b,d,3,device=keys.device,dtype=keys.dtype)
            for i in range(n):
                prediction=torch.einsum('bd,bdv->bv',keys[:,i],state)
                state=state+keys[:,i,:,None]*(values[:,i]-prediction)[:,None,:]
            predicted=torch.einsum('bqd,bdv->bqv',queries,state)
            logits=-(predicted[...,None,:]-self.codes).square().sum(dim=-1)
        elif self.kind=='rls':
            covariance=keys.transpose(1,2)@keys+1e-3*torch.eye(d,device=keys.device,dtype=keys.dtype)
            state=torch.linalg.solve(covariance,keys.transpose(1,2)@values)
            predicted=queries@state;logits=-(predicted[...,None,:]-self.codes).square().sum(dim=-1)
        else:
            if self.kind=='attention_fifo':keys=keys[:,-26:];values=values[:,-26:]
            weights=torch.softmax(20*(queries@keys.transpose(1,2)),dim=-1)
            predicted=weights@values;logits=-(predicted[...,None,:]-self.codes).square().sum(dim=-1)
        return logits*self.log_temperature.exp().clamp(max=100),compatibility
    def state_scalars(self,bindings):
        d=self.width
        if self.kind=='atom':return d*16
        if self.kind in ('mean','mean_same_width'):return d*4
        if self.kind=='delta':return d*3
        if self.kind=='rls':return d*d+d*3
        return min(bindings,26)*19 if self.kind=='attention_fifo' else bindings*19
