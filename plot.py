"""Draw results/generalization.{png,svg} and results/motors.{png,svg}."""

import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = {  # fixed categorical order; marker shape is a second, color-free cue
    "extrusion": ("Extrusion (two-arm)", "#2a78d6", "o"),
    "one_arm": ("One-arm ablation", "#eb6834", "s"),
    "attn_rope": ("Attention + RoPE", "#1baf7a", "^"),
    "attn_alibi": ("Attention + ALiBi", "#eda100", "D"),
    "bilstm": ("BiLSTM", "#e87ba4", "v"),
}
BLUES = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
PANELS = [
    ((4000, 3), "4,000 training sequences, 3 layers"),
    ((300, 3), "300 training sequences, 3 layers"),
    ((4000, 1), "4,000 training sequences, 1 layer"),
]

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK_2, "axes.labelcolor": INK_2, "xtick.color": INK_2,
    "ytick.color": INK_2, "text.color": INK, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "svg.fonttype": "none",
})


def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(f"results/{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def generalization():
    runs = [json.loads(l) for l in Path("results/runs.jsonl").read_text().splitlines() if l.strip()]
    groups = defaultdict(list)
    for r in runs:
        groups[(r["mixer"], r["n_train"], r.get("layers", 3))].append(r)
    lengths = [64, 128, 256, 512]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
    for ax, (setting, title) in zip(axes, PANELS):
        for mixer, (label, color, marker) in SERIES.items():
            rs = groups.get((mixer, *setting))
            if not rs:
                continue
            acc = np.array([[r["test_acc"][str(n)] for n in lengths] for r in rs])
            ax.fill_between(lengths, acc.min(0), acc.max(0), color=color, alpha=0.18, linewidth=0)
            ax.plot(lengths, acc.mean(0), color=color, linewidth=2, marker=marker, markersize=6,
                    markeredgecolor=SURFACE, markeredgewidth=1.5, label=label)
        ax.set_xscale("log", base=2)
        ax.set_xticks(lengths, [str(n) for n in lengths])
        ax.minorticks_off()
        ax.set_ylim(0.5, 1.01)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.grid(axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xlabel("Test length (trained on 24-64)")
    axes[0].set_ylabel("Promoter accuracy")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.06))
    fig.text(0.5, -0.04, "Lines: mean of 3 seeds. Bands: min-max across seeds. The 1-layer panel compares only the two extrusion variants.",
             ha="center", color=INK_2, fontsize=9)
    save(fig, "generalization")


def motors():
    tables = json.loads(Path("results/motors.json").read_text())
    layer1 = tables[0]
    tokens = list(layer1)
    names = {".": ". filler", "E1": "E1 enhancer", "E2": "E2 enhancer", "E3": "E3 enhancer",
             "P": "P promoter", ">": "> forward CTCF", "<": "< reverse CTCF"}
    heads = len(layer1[tokens[0]])
    grid = np.array([[layer1[t][h][k] for h in range(heads) for k in range(3)] for t in tokens])
    fig, ax = plt.subplots(figsize=(10, 3.9))
    cmap = LinearSegmentedColormap.from_list("blues", BLUES)
    ax.imshow(grid, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = grid[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8.5, color="#ffffff" if v > 0.55 else INK)
    sub = ["stall\nleft arm", "stall\nright arm", "load"] * heads
    ax.set_xticks(range(len(sub)), sub, fontsize=8)
    ax.set_yticks(range(len(tokens)), [names[t] for t in tokens])
    ax.tick_params(length=0)
    for h in range(heads):
        ax.text(3 * h + 1, -0.9, f"Motor {h + 1}", ha="center", va="bottom", fontsize=10, color=INK)
        if h:
            ax.axvline(3 * h - 0.5, color=SURFACE, linewidth=4)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title("Layer 1 learned motor probabilities, averaged by token", loc="left", fontsize=10, pad=34)
    fig.text(0.125, -0.03, "Never told which tokens are barriers: motors 1, 2 and 4 learned that '>' stalls left arms and '<' stalls right arms.",
             color=INK_2, fontsize=9)
    save(fig, "motors")


if __name__ == "__main__":
    generalization()
    motors()
