"""Display reference sample requirements and sensitivity-contract failures."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
rows = json.loads((ROOT/'results/counted_loading.json').read_text())['settings']
budgets = [500, 2000, 8000, 32000]
normal = [[r for r in rows if not r['stress'] and r['budget'] == n] for n in budgets]
plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
low = [100*min(r['median_width'] for r in group) for group in normal]
high = [100*max(r['median_width'] for r in group) for group in normal]
axes[0].fill_between(budgets, low, high, alpha=.2, color='#2364aa', label='Range across nine cases')
axes[0].plot(budgets, high, 'o-', color='#2364aa', label='Largest median width')
axes[0].plot(budgets, low, 'o-', color='#2364aa', alpha=.5)
axes[0].axhline(10, color='#555555', linestyle='--', linewidth=1, label='Frozen precision threshold')
axes[0].set(xscale='log', xlabel='Reference episodes per arm', ylabel='Median interval width (percentage points)',
            title='A  Counts narrow the target interval')
axes[0].legend(fontsize=8)
for method, color, label in [('common_sensitivity', '#b85b26', 'Assume equal sensitivity (wrong)'),
                             ('bounded_sensitivity', '#49906e', 'Allow sensitivity ratio in [0.9, 1]')]:
    coverage = [min(r['coverage'] for r in rows if r['stress'] and r['method'] == method and r['budget'] == n) for n in budgets]
    axes[1].plot(budgets, coverage, 'o-', color=color, label=label)
axes[1].axhline(.95, color='#555555', linestyle='--', linewidth=1)
axes[1].set(xscale='log', ylim=(0, 1.05), xlabel='Reference episodes per arm', ylabel='Minimum observed coverage across nine cases',
            title='B  Reporter misspecification defeats coverage')
axes[1].legend(fontsize=8, loc='lower left')
fig.suptitle('Simulated counts; 24,000 retention episodes per arm; 100 repetitions per case', fontsize=11)
fig.savefig(ROOT/'results/counted_loading.png', dpi=190)
fig.savefig(ROOT/'results/counted_loading.svg')
