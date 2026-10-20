"""Show the transparency-dependent probe budget and length transfer."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np
from colored_recovery import choose_colors

ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results/colored_recovery_study.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
ceilings=[.2,.5,.8,.95];budgets=[2*choose_colors(p,4096)[0] for p in ceilings]
axes[0].bar([str(p) for p in ceilings],budgets,color='#146a87')
for i,n in enumerate(budgets):axes[0].text(i,n+5,str(n),ha='center',fontsize=10)
axes[0].set(xlabel='Maximum pass probability',ylabel='Vector responses',title='A uniform guarantee has a transparency cost',ylim=(0,550))
for rho in ceilings:
    samples=[[x['gate_rmse'] for x in r['rows'] if x['size']==n and x['rho']==rho and x['noise_sd']==.001] for n in (32,128,512)]
    means=[np.mean(x) for x in samples];sem=[np.std(x,ddof=1)/np.sqrt(len(x)) for x in samples]
    axes[1].errorbar([32,128,512],means,yerr=sem,fmt='o-',capsize=2,label=f'ceiling {rho}')
axes[1].set(xlabel='Number of loads',ylabel='Mean gate RMSE ± SEM',xscale='log',title='Noisy recovery; five mechanisms per setting')
axes[1].set_xticks([32,128,512],labels=['32','128','512']);axes[1].xaxis.set_minor_locator(NullLocator())
axes[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0));axes[1].grid(alpha=.2);axes[1].legend(fontsize=8)
fig.suptitle('Controlled unsaturated operator, known loading and source-invariant gates',fontsize=12)
fig.savefig(ROOT/'results/colored_recovery_study.png',dpi=170)
