"""Draw results/generalization, results/generality and results/motors (png and svg)."""

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
    "extrusion_linear": ("Extrusion, linear O(n)", "#008300", "X"),
    "one_arm": ("One-arm ablation", "#eb6834", "s"),
    "attn_rope": ("Attention + RoPE", "#1baf7a", "^"),
    "attn_alibi": ("Attention + ALiBi", "#eda100", "D"),
    "bilstm": ("BiLSTM", "#e87ba4", "v"),
}
BLUES = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
PANELS = [
    (("genome", 4000, 3), "4,000 training sequences, 3 layers"),
    (("genome", 300, 3), "300 training sequences, 3 layers"),
    (("genome", 4000, 1), "4,000 training sequences, 1 layer"),
]
GENERALITY = [
    (("dyck", 4000, 3), "Nested brackets (Dyck): innermost open bracket", (0.8, 1.005), None),
    (("recall", 4000, 3), "Associative recall: value of the queried key", (0.0, 1.02), 0.125),
]
LENGTHS = [64, 128, 256, 512]

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK_2, "axes.labelcolor": INK_2, "xtick.color": INK_2,
    "ytick.color": INK_2, "text.color": INK, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "svg.fonttype": "none",
})


def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(f"results/{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def load_groups():
    runs = [json.loads(l) for l in Path("results/runs.jsonl").read_text().splitlines() if l.strip()]
    groups = defaultdict(list)
    for r in runs:
        groups[(r["mixer"], r.get("task", "genome"), r["n_train"], r.get("layers", 3))].append(r)
    return groups


def draw_panel(ax, groups, setting, title):
    for mixer, (label, color, marker) in SERIES.items():
        rs = groups.get((mixer, *setting))
        if not rs:
            continue
        acc = np.array([[r["test_acc"][str(n)] for n in LENGTHS] for r in rs])
        ax.fill_between(LENGTHS, acc.min(0), acc.max(0), color=color, alpha=0.18, linewidth=0)
        dashed = (0, (5, 2.5)) if mixer == "extrusion_linear" else "-"  # it often sits on top of extrusion
        ax.plot(LENGTHS, acc.mean(0), color=color, linewidth=2, linestyle=dashed, marker=marker, markersize=6,
                markeredgecolor=SURFACE, markeredgewidth=1.5, label=label)
    ax.set_xscale("log", base=2)
    ax.set_xticks(LENGTHS, [str(n) for n in LENGTHS])
    ax.minorticks_off()
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("Test length (trained on 24-64)")


def shared_legend(fig, axes):
    seen = {}
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            seen.setdefault(l, h)
    order = [v[0] for v in SERIES.values() if v[0] in seen]
    fig.legend([seen[l] for l in order], order, loc="upper center", ncol=len(order), frameon=False,
               bbox_to_anchor=(0.5, 1.07), fontsize=9)


def generalization():
    groups = load_groups()
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    for ax, (setting, title) in zip(axes, PANELS):
        draw_panel(ax, groups, setting, title)
    axes[0].set_ylim(0.5, 1.01)
    axes[0].set_ylabel("Promoter accuracy")
    shared_legend(fig, axes)
    fig.text(0.5, -0.04, "Toy genome task. Lines: mean of 3 seeds. Bands: min-max across seeds. Panels show only the mixers run in that setting.",
             ha="center", color=INK_2, fontsize=9)
    save(fig, "generalization")


def generality():
    groups = load_groups()
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
    for ax, (setting, title, ylim, chance) in zip(axes, GENERALITY):
        draw_panel(ax, groups, setting, title)
        ax.set_ylim(*ylim)
        if chance:
            ax.axhline(chance, color=INK_2, linewidth=1, linestyle=(0, (4, 3)))
            ax.text(LENGTHS[0], chance + 0.02, "chance (1 in 8)", ha="left", color=INK_2, fontsize=8.5)
    axes[0].set_ylabel("Accuracy")
    shared_legend(fig, axes)
    fig.text(0.5, -0.04, "Tasks the model was not designed for. 4,000 training sequences, 3 layers, 3 seeds; bands show min-max across seeds.",
             ha="center", color=INK_2, fontsize=9)
    save(fig, "generality")


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
    generality()
    motors()
