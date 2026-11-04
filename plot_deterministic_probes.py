"""Show the restricted global guarantee separately from observed recovery."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from deterministic_probes import recovery_constants

ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results/deterministic_probe_study.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
x=np.linspace(0,.09,200);axes[0].plot(x,[recovery_constants(a)['lower'] for a in x],label='Proved global lower constant',color='#146a87')
rho=[0,.02,.05,.08]
observed=[min(a['smallest_singular_value'] for a in r['local_and_secant_checks'] if a['rho']==p and a['design']=='coded') for p in rho]
axes[0].plot(rho,observed,'o--',color='#ba692d',label='Smallest observed local singular value')
axes[0].axhline(0,color='black',lw=.8);axes[0].set(xlabel='Maximum allowed pass probability',ylabel='Lower sensitivity',title='Proof and numerical evidence are distinct')
axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
for kind,color in (('coded','#146a87'),('random','#ba692d')):
    samples=[[a['gate_rmse'] for a in r['free_gate_fits'] if a['size']==n and a['noise_sd']==.001 and a['design']==kind] for n in (8,16,32)]
    means=[np.mean(x) for x in samples];sem=[np.std(x,ddof=1)/np.sqrt(len(x)) for x in samples]
    axes[1].errorbar([8,16,32],means,yerr=sem,fmt='o-',capsize=2,color=color,label=kind)
axes[1].set(xlabel='Number of loads',ylabel='Mean gate RMSE ± SEM',title='Four noisy responses; pass probabilities ≤ .05')
axes[1].legend();axes[1].grid(alpha=.2)
fig.suptitle('Global strong-barrier theorem and exploratory recovery at noise SD 0.001',fontsize=12)
fig.savefig(ROOT/'results/deterministic_probe_study.png',dpi=170)
