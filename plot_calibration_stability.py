"""Plot retained stability bounds and noisy-design outcomes."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent


def main():
    budget=json.loads((ROOT/'results/calibration_budget_study.json').read_text())
    noisy=json.loads((ROOT/'results/robust_identification_study.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for count,placement,label,color in [(1,'boundary_cluster','One boundary calibration','#b35437'),
            (8,'boundary_cluster','Eight clustered calibrations','#9b6c9d'),
            (8,'evenly_spaced','Eight spaced calibrations','#247c8d')]:
        rows=[r for r in budget['rows'] if r['pass_rate']==.6 and r['budget']==count and r['placement']==placement]
        axes[0].plot([r['gates'] for r in rows],[r['mean_variance_over_sigma_squared'] for r in rows],
                     'o-',label=label,color=color)
    axes[0].set(xlabel='Gate pairs',ylabel='Mean local variance bound / noise variance',
                title='Exact identification can remain unstable')
    axes[0].legend(frameon=False,fontsize=9)
    methods=[('boundary','Boundary','#b35437'),('largest_response_span','Response span','#9b6c9d'),
             ('uniform_random_site','Random site','#777777'),('local_noise_aware','Local design','#b29131'),
             ('grid_robust_local','Grid-robust design','#247c8d')]
    for method,label,color in methods:
        rows=[r for r in noisy['summary'] if r['sigma']==.01]
        axes[1].plot([r['size'] for r in rows],[r['strategies'][method]['mean_heldout_rmse'] for r in rows],
                     'o-',label=label,color=color)
    axes[1].set(xlabel='Contact-matrix size',ylabel='Mean held-out intervention RMSE',
                title='No reliable one-intervention winner',xticks=[8,16,32])
    axes[1].legend(frameon=False,fontsize=9)
    for ax in axes:ax.grid(alpha=.18)
    fig.text(.02,.02,'Left: local Gaussian oracle model, pass rate 0.6. Right: 10 fresh synthetic mechanisms per size, noise 0.01; exploratory.',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,1))
    fig.savefig(ROOT/'results/calibration_stability.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
