"""Render retained matched and loading-drift results."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
common = json.loads((root/'results/loading_target.json').read_text())
drift = json.loads((root/'results/loading_drift.json').read_text())
refined = json.loads((root/'results/timestamp_refinement.json').read_text())
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
x = np.arange(3)
methods = [('absolute_uncertain', 'Absolute, uncertain f', '#397c8c'), ('absolute_known_factor', 'Absolute, exact f supplied', '#799669')]
for i, (method, label, color) in enumerate(methods):
    chosen = sorted([r for r in common['settings'] if not r['stress'] and r['n'] == 24000 and r['method'] == method], key=lambda r:r['mechanism'])
    axes[0].bar(x+(i-1)*.25, [100*r['median_width'] for r in chosen], .25, label=label, color=color)
axes[0].bar(x+.25, [100*r['median_width'] for r in refined['settings']], .25, label='Timestamp, refined outer', color='#be9153')
axes[0].axhline(10, color='black', ls='--', lw=1)
axes[0].set_xticks(x, ['M1', 'M2', 'M3'])
axes[0].set_ylabel('Median interval width (percentage points)')
axes[0].set_ylim(0, 14)
axes[0].set_title('A  Matched common-loading target', loc='left', fontsize=10)
axes[0].legend(fontsize=7)
for axis, method, title in [(axes[1], 'broad_ratio', 'B  Loading ratio allowed in [.8,1]'), (axes[2], 'calibrated_ratio', 'C  Extra ratio calibration within .02')]:
    matrix = np.array([[next(r['median_width'] for r in drift['settings'] if r['method'] == method and r['mechanism'] == m and r['ratio'] == k)*100 for k in [1., .9, .8]] for m in range(3)])
    im = axis.imshow(matrix, cmap='YlOrBr', vmin=0, vmax=16, aspect='auto')
    axis.set_xticks(range(3), ['1.0', '0.9', '0.8']); axis.set_yticks(range(3), ['F1', 'F2', 'F3'])
    axis.set_xlabel('True treatment/baseline loading ratio')
    axis.set_title(title, loc='left', fontsize=10)
    for i in range(3):
        for j in range(3):
            axis.text(j, i, f'{matrix[i,j]:.1f}', ha='center', va='center', color='white' if matrix[i,j] > 11 else 'black')
fig.colorbar(im, ax=axes[1:], label='Median interval width (pp)', shrink=.8)
fig.suptitle('h=1, 24,000 observations; A uses original mechanisms, B/C fresh mechanisms; calibration/timing costs differ', fontsize=10)
fig.savefig(root/'results/loading_research.png', dpi=180)
