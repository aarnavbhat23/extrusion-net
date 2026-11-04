"""Independent-engine errors and retained mechanism rejections."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results/external_engine_study.json').read_text())['rows']
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
matched=[r for r in rows if r['mode']=='matched']
for family,label in [('asynchronous','Asynchronous'),('synchronized','Synchronized'),('independent','Independent arms')]:
    means=[np.mean([t['rmse'][family] for r in matched if r['removal']==mu for t in r['targets']]) for mu in (.15,.5,1.)]
    axes[0].plot([.15,.5,1.],means,marker='o',label=label)
axes[0].set_xlabel('Removal / step rate');axes[0].set_ylabel('Mean reach-probability RMSE')
axes[0].set_title('Unmodified external simulator');axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
labels=['Matched','Capture-dependent\nresidence','Inside-facing\neviction']
for i,mode in enumerate(('matched','protected','eviction')):
    selected=[r for r in rows if r['mode']==mode];fraction=np.mean([r['rejected'] for r in selected])
    axes[1].bar(i,fraction,color='tab:blue' if i==0 else 'tab:orange')
    axes[1].text(i,fraction+.025,f"{sum(r['rejected'] for r in selected)}/{len(selected)}",ha='center')
axes[1].set_xticks(range(3),labels,fontsize=8);axes[1].set_ylim(0,1.15);axes[1].set_ylabel('Fraction rejected on fresh targets')
axes[1].set_title('Model-conditional interval check')
fig.suptitle('Restricted reach transfer: known source and two known barrier locations',fontsize=12)
fig.savefig(ROOT/'results/external_engine.png',dpi=160)
