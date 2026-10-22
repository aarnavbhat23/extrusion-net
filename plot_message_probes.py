"""Show conditioning separately from shared-controller generalization."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results/message_probe_study.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
for kind,color in (('alternating','#146a87'),('random','#ba692d')):
    for probes,style in ((1,'--'),(2,'-')):
        med=[]
        for size in (8,16,32,64):
            rows=[x['mean_variance_over_noise_variance'] for x in r['local_information']
                  if x['pattern']==kind and x['values_per_load']==probes and x['size']==size]
            med.append(np.median([x for x in rows if x is not None]))
        axes[0].plot([8,16,32,64],med,style,marker='o',color=color,label=f'{kind}, {2*probes} responses')
axes[0].set(xscale='log',yscale='log',xlabel='Number of loads',ylabel='Median local variance / noise variance',
            title='Free gates: identification can be unstable')
axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
for i,kind in enumerate(('uniform','random','alternating')):
    values=[x['test_gate_rmse'] for x in r['shared_controller'] if x['pattern']==kind]
    axes[1].scatter(np.full(10,i),values,alpha=.65,s=28)
    axes[1].plot([i-.25,i+.25],[np.mean(values)]*2,color='black',lw=2)
axes[1].set(xticks=[0,1,2],xticklabels=['Uniform','Random','Alternating'],
            ylabel='Unseen-context gate RMSE',title='Shared controller: 68 observed scalars')
axes[1].grid(axis='y',alpha=.2)
fig.suptitle('Exploratory synthetic probes — ten mechanisms/seeds per setting',fontsize=12)
fig.savefig(ROOT/'results/message_probe_study.png',dpi=170)
