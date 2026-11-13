"""Show complete intervention outcomes and retained inference failures."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
point=json.loads((ROOT/'results/occupancy_study.json').read_text())['rows']
bank=json.loads((ROOT/'results/bank_acquisition_study.json').read_text())['rows']
continuous=json.loads((ROOT/'results/continuous_occupancy_study.json').read_text())['rows']
refined=json.loads((ROOT/'results/continuous_occupancy_refinement.json').read_text())['rows']
witness=json.loads((ROOT/'results/occupancy_ambiguity_diagnostic.json').read_text())['rows']
fig,axes=plt.subplots(1,3,figsize=(14,4.6))
x=[r['target_rmse']['published_dlem'] for r in point];y=[r['target_rmse']['capture'] for r in point]
for on_bank,color,label in ((True,'#168376','Bank members (18)'),(False,'#ba743e','Off-grid (3)')):
    ids=[i for i,r in enumerate(point) if r['on_finite_bank']==on_bank]
    axes[0].scatter(np.array(x)[ids],np.array(y)[ids],c=color,s=28,label=label)
axes[0].legend(fontsize=7,loc='upper left');limit=max(max(x),max(y))*1.08
axes[0].plot([0,limit],[0,limit],color='gray',lw=1,ls='--');axes[0].set_xlim(0,limit);axes[0].set_ylim(0,limit)
axes[0].set_xlabel('Published recurrence target RMSE');axes[0].set_ylabel('Capture-model target RMSE')
axes[0].set_title('All 21 cases: median ratio 1.93×\nFrozen 2× gate fails',fontsize=10)
policies=['round_robin','global_disagreement','target_disagreement']
cost=[np.mean([r['result']['quoted_cost'] for r in bank if r['on_bank'] and r['policy']==p]) for p in policies]
axes[1].bar(range(3),cost,color=['#8d979c','#49658b','#168376']);axes[1].set_xticks(range(3),['Round robin','Global','Target'],rotation=15)
axes[1].set_ylabel('Mean hypothetical observation cost');axes[1].set_title('Target policy fails cost gate\nAll 54 on-grid cases attain precision',fontsize=10)
indices=np.arange(6)
axes[2].plot(indices,[r['continuous']['max_target_width']*100 for r in continuous],'o-',label='Original outer',color='#8d979c')
axes[2].plot(indices,[r['continuous']['max_target_width']*100 for r in refined],'o-',label='Refined outer',color='#49658b')
axes[2].plot(indices,[r['attained_absence_span']*100 for r in witness],'o-',label='Feasible absence span',color='#168376')
axes[2].axhline(2,color='#a44848',ls='--',lw=1,label='2-point target')
axes[2].set_xticks(indices,['1a','1b','2a','2b','3a','3b']);axes[2].set_xlabel('Fresh mechanism and repeat')
axes[2].set_ylabel('Probability width (percentage points)');axes[2].set_ylim(bottom=0)
axes[2].set_title('Continuous-domain validity, weak precision\nRefinement is the same-data reanalysis',fontsize=10)
axes[2].legend(fontsize=7)
for axis in axes:axis.grid(axis='y',alpha=.2)
fig.suptitle('Controlled random-delay snapshots: useful diagnostics, failed advancement gates',fontsize=12)
fig.text(.5,.012,'Single-motor simulator studies; no measured chromatin, general dLEM superiority or foundational learning claim.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.04,1,.94));fig.savefig(ROOT/'results/occupancy_research.png',dpi=180)
