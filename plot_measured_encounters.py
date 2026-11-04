import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

r=json.loads(Path('results/measured_encounters.json').read_text())['analyses']['primary_N']
fig,ax=plt.subplots(1,2,figsize=(10.5,4),layout='constrained')
models=list(r['mean_metrics']);labels=['Constant','Published force law','Force logistic','Force cubic','State logistic','Same-input boosting']
values=[r['mean_metrics'][m]['log_loss'] for m in models]
ax[0].barh(labels[::-1],values[::-1],color=['#21618c' if m!='state_logistic' else '#b03a2e' for m in models][::-1])
ax[0].set(xlabel='Equal-batch, molecule-weighted log loss (lower is better)',title='Measured encounter transfer: candidate fails')
folds=r['folds'];baseline=np.array([f['metrics']['published_exponential']['log_loss'] for f in folds]);candidate=np.array([f['metrics']['state_logistic']['log_loss'] for f in folds])
for i in range(len(folds)):ax[1].plot([0,1],[baseline[i],candidate[i]],marker='o',color='#999999',alpha=.8)
ax[1].plot([0,1],[baseline.mean(),candidate.mean()],marker='o',color='#b03a2e',lw=3,label='Equal-batch mean')
ax[1].set(xticks=[0,1],xticklabels=['Published force law','State logistic'],ylabel='Held-out batch log loss',title='Each line is a held-out experimental batch')
ax[1].legend(fontsize=8)
fig.savefig('results/measured_encounters.png',dpi=180)
