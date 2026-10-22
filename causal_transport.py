"""Causal transport with a finite shared-age correction and matched controls."""
import torch
from torch import nn
from torch.nn import functional as F


def delay(x,steps):
    if steps==0:return x
    if steps>=x.shape[-2]:return torch.zeros_like(x)
    return F.pad(x[...,:-steps,:],(0,0,steps,0))


def causal_sum(values,gates,window=None):
    """Exclusive previous-position gates; associative doubling, no pair matrix."""
    limit=values.shape[-2] if window is None else min(window,values.shape[-2])
    total=torch.zeros_like(values);path=torch.ones_like(gates[...,None]);offset=0
    block=values;product=delay(gates[...,None],1);size=1
    while size<=limit:
        if limit&size:
            total=total+path*delay(block,offset)
            path=path*delay(product,offset);offset+=size
        if size*2<=limit:
            block=block+product*delay(block,size)
            product=product*delay(product,size)
        size*=2
    return total


def causal_contacts(load,left,right,values,age=None):
    inner=causal_sum(delay(values,1),left,None if age is None else age+1)
    return causal_sum(load[...,None]*inner,right,None if age is None else age+1)


def clock_messages(load,left,right,values,survival,max_age,shared=True):
    """Capped geometric clock. Independent-clock control has the same reach cap."""
    if not shared:
        return causal_contacts(load,left*survival,right*survival,values,max_age)
    result=torch.zeros_like(values)
    for age in range(max_age+1):
        weight=survival**age*((1-survival) if age<max_age else 1.)
        result=result+weight[...,None]*causal_contacts(load,left,right,values,age)
    return result


class CausalTransportMixer(nn.Module):
    def __init__(self,d,heads,variant='static',max_age=8):
        super().__init__();self.heads=heads;self.variant=variant;self.max_age=max_age
        self.control=nn.Linear(d,3*heads);self.value=nn.Linear(d,d);self.out=nn.Linear(d+heads,d)
        self.clock_raw=nn.Parameter(torch.zeros(heads));self.correction_raw=nn.Parameter(torch.full((heads,),-2.))
        with torch.no_grad():self.control.bias.view(heads,3)[:,:2]=-2.

    def forward(self,x):
        batch,n,d=x.shape
        z=self.control(x).view(batch,n,self.heads,3).permute(0,2,3,1)
        left=torch.sigmoid(-z[:,:,0]);right=torch.sigmoid(-z[:,:,1]);load=torch.sigmoid(z[:,:,2])
        values=self.value(x).view(batch,n,self.heads,d//self.heads).transpose(1,2)
        values=torch.cat([values,torch.ones_like(values[...,:1])],dim=-1)
        mixed=causal_contacts(load,left,right,values)
        if self.variant!='static':
            survival=torch.sigmoid(self.clock_raw)[None,:,None]
            correction=clock_messages(load,left,right,values,survival,self.max_age,self.variant=='shared')
            mixed=mixed+F.softplus(self.correction_raw)[None,:,None,None]*correction
        mass=mixed[...,-1:];message=mixed[...,:-1]/(1+mass)
        message=message.transpose(1,2).reshape(batch,n,d)
        return self.out(torch.cat([message,torch.log1p(mass.squeeze(-1)).transpose(1,2)],dim=-1))


class CausalAttention(nn.Module):
    def __init__(self,d,heads):
        super().__init__();self.heads=heads;self.qkv=nn.Linear(d,3*d);self.out=nn.Linear(d,d)
    def forward(self,x):
        b,n,d=x.shape;q,k,v=self.qkv(x).view(b,n,3,self.heads,d//self.heads).permute(2,0,3,1,4)
        y=F.scaled_dot_product_attention(q,k,v,is_causal=True)
        return self.out(y.transpose(1,2).reshape(b,n,d))


class CausalModel(nn.Module):
    def __init__(self,variant,vocab,classes,d=32,heads=4,max_age=8):
        super().__init__();self.embed=nn.Embedding(vocab,d);self.conv=nn.Conv1d(d,d,5)
        self.norms=nn.ModuleList([nn.LayerNorm(d) for _ in range(4)])
        self.mixer=CausalAttention(d,heads) if variant=='attention' else CausalTransportMixer(d,heads,variant,max_age)
        self.mlp=nn.Sequential(nn.Linear(d,2*d),nn.GELU(),nn.Linear(2*d,d));self.head=nn.Linear(d,classes)
    def forward(self,tokens):
        x=self.embed(tokens);x=x+F.gelu(self.conv(F.pad(self.norms[0](x).transpose(1,2),(4,0))).transpose(1,2))
        x=x+self.mixer(self.norms[1](x));x=x+self.mlp(self.norms[2](x))
        return self.head(self.norms[3](x))
