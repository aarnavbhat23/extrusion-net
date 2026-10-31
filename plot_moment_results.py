"""Plot all-seed means and ranges for the stopped retrieval candidate."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results/moment_retrieval_study.json').read_text())['rows']
loads=[8,16,32,64]
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for kind,label in [('atom','Moment atom'),('mean','Averaging, matched state'),('delta','Delta rule'),('attention','Full attention')]:
    values=np.array([[r['tests'][str(n)]['accuracy']*100 for n in loads] for r in rows if r['kind']==kind])
    line=axes[0].plot(loads,values.mean(axis=0),marker='o',label=label)[0]
    axes[0].fill_between(loads,values.min(axis=0),values.max(axis=0),alpha=.15,color=line.get_color())
coverage=np.array([[r['tests'][str(n)]['certified_coverage']*100 for n in loads] for r in rows if r['kind']=='atom'])
axes[1].plot(loads,coverage.mean(axis=0),marker='o',color='tab:blue',label='Moment atom')
axes[1].fill_between(loads,coverage.min(axis=0),coverage.max(axis=0),alpha=.15,color='tab:blue')
for ax,gate,title,y in zip(axes,[95,90],['Unseen-key retrieval','Conditional certified coverage'],['Point accuracy (%)','Queries answered (%)']):
    ax.axhline(gate,color='0.35',linestyle='--',linewidth=1,label=f'Frozen gate: {gate}%')
    ax.set_xscale('log',base=2);ax.set_xticks(loads,labels=loads);ax.set_ylim(0,104)
    ax.set_xlabel('Stored associations (trained on 8)');ax.set_ylabel(y);ax.set_title(title)
    ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('The retrieval candidate fails the greater-load gate\nThree seeds; shading shows seed range, not a confidence interval',fontsize=12)
fig.savefig(ROOT/'results/moment_retrieval.png',dpi=160)
