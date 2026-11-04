"""Plot timing-family and speed-flexibility tradeoffs without hiding noise failures."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
first=json.loads((ROOT/'results/asynchronous_transport_study.json').read_text())['rows']
second=json.loads((ROOT/'results/asymmetric_clock_study.json').read_text())['rows']
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
noise=[.0001,.001,.01]
for kind,label in [('selected','Band-selected timing'),('synchronized','Forced synchronized'),('static_oracle','Static marginal oracle')]:
    values=[]
    for sigma in noise:
        rows=[r for r in first if r['noise_sd']==sigma and r['true_family']=='asynchronous']
        values.append(np.mean([r['selected_heldout_rmse'] if kind=='selected' else r['heldout_rmse'][kind] for r in rows]))
    axes[0].plot(noise,values,marker='o',label=label)
for kind,label in [('unequal','Unequal-speed fit'),('equal','Equal-speed fit, same observations')]:
    values=[np.mean([r['heldout_rmse'][kind] for r in second if r['noise_sd']==sigma and r['right_rate']!=1]) for sigma in noise[:2]]
    axes[1].plot(noise[:2],values,marker='o',label=label)
for ax in axes:
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Band noise SD');ax.set_ylabel('Mean held-out message RMSE')
    ax.grid(alpha=.2);ax.legend(fontsize=8)
axes[0].set_title('Timing selection: payoff fades at high noise')
axes[1].set_title('Extra speed freedom: moderate-noise gate fails')
fig.suptitle('Conditional transport learning: means over development mechanisms',fontsize=12)
fig.savefig(ROOT/'results/asynchronous_transport.png',dpi=160)
