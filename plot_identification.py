"""Render the structural illustration and retained intervention study."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from study_identification import contact_array,gate_family,observable_bands

ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261021)
q=rng.uniform(.4,.9,15);left=rng.uniform(.4,.7,14);right=rng.uniform(.4,.7,14)
q,s,t=observable_bands(contact_array(q,left,right));l,r=gate_family(q,s,t)
report=json.loads((ROOT/'results/identification_study.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
for k in (0,5,10):axes[0].plot(r[:,0],l[:,k]*q[k+1],label=f'Cut gate {k}')
axes[0].set(xlabel='First right gate across equivalent explanations',ylabel='Predicted local contact after directional cut',title='Same baseline contacts, different interventions')
axes[0].legend(fontsize=8)
for name,label in [('boundary','Boundary'),('uniform_random_site','Random'),('midpoint_sensitivity','Local sensitivity'),('largest_response_span','Response range')]:
    rows=[x for x in report['summary'] if x['sigma']==.01 and x['strategy']==name]
    axes[1].plot([x['length'] for x in rows],[x['mean_heldout_intervention_rmse'] for x in rows],marker='o',label=label)
axes[1].set(xlabel='Sequence length',ylabel='Held-out local intervention RMSE',title='One noisy calibration; exact baseline map')
axes[1].legend(fontsize=8)
fig.savefig(ROOT/'results/identification.png',dpi=180)
