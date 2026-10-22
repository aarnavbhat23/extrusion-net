import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=json.loads(Path('results/pulse_order.json').read_text())
fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
for d,color in enumerate(['#777777','#21618c','#b03a2e']):
 rows=[x for x in r['memory'] if x['degree_balanced']==d]
 ax[0].plot(range(1,7),[x['contrast'] for x in rows],marker='o',color=color,label=f'{2**(d+1)} slots')
 ax[1].plot(range(1,7),[x['robust_known_corruption'] for x in rows],marker='o',color=color,label=f'{2**(d+1)} slots')
ax[0].set(xlabel='Synthetic mechanism (fixed order)',ylabel='Endpoint survival contrast',title='Balancing ageing can suppress the state signal')
ax[0].legend(fontsize=8)
ref=json.loads(Path('results/pulse_test_refinement.json').read_text())
exact=[x['robust_fisher_rejection'] for x in ref['memory'] if x['degree_balanced']==1]
ax[1].plot(range(1,7),exact,marker='x',ls='--',color='#21618c',label='4 slots, exact test (development)')
ax[1].axhline(.8,color='black',ls='--',label='Frozen power gate')
ax[1].set(xlabel='Synthetic mechanism (fixed order)',ylabel='Power with bounded recording errors',ylim=(-.03,1.05),title='Primary four-slot design: 0/6 reach 80%')
ax[1].legend(fontsize=8)
fig.savefig('results/pulse_order.png',dpi=180)
