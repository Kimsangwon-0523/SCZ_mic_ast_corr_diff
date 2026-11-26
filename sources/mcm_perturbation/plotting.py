"""Draw sample scores and comparison brackets for each astrocyte dataset."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

from .effects import CONDITIONS, summarize_score_effects

COLORS = {"Resting-MCM": "#526C82", "LPS-MCM": "#C05330", "Hemin-MCM": "#8166A0"}
INK = "#202D3A"
MUTED = "#596675"
PAIR = "#A3ADB8"
RULE = "#DFE5EB"
ASTROCYTE_PROGRAMS = ("ACHES", "SNAP-a")
FIGURE_NUMBERS = {"GSE282665": "24", "GSE228927": "25"}


def style_axis(ax, axis_label):
    ax.set_ylabel(axis_label, fontsize=10.5, color=INK, labelpad=9)
    ax.margins(y=0.19)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
    ax.yaxis.set_major_formatter(FuncFormatter(
        lambda value, _: f"{value:.3f}".rstrip("0").rstrip(".").replace("-", "−")
    ))
    ax.yaxis.grid(True, color=RULE, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", direction="out", length=3, width=0.9, pad=5, colors=INK)
    ax.tick_params(axis="y", labelcolor=MUTED, labelsize=10)
    ax.spines[["right", "top"]].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set(color=INK, linewidth=1.15)


def draw_panel(ax, scores, letter, title, conditions, axis_label):
    """Draw one point per sample, a mean for each condition, and documented pairs."""
    if set(scores["condition"]) != set(conditions):
        raise ValueError(f"Unexpected or missing conditions for {title}")
    positions = {}
    for i, condition in enumerate(conditions):
        rows = scores.loc[scores["condition"].eq(condition)].sort_values("sample_id")
        offsets = np.linspace(-0.09, 0.09, len(rows)) if len(rows) > 1 else [0.0]
        x = i + np.asarray(offsets)
        positions.update(zip(rows["sample_id"], x))
        color = COLORS[condition]
        ax.scatter(x, rows["score"], c=color, s=67, linewidth=0.9, edgecolors="white", zorder=4)
        ax.hlines(rows["score"].mean(), i - 0.21, i + 0.21, color=INK, linewidth=2.1, zorder=3)

    for _, rows in scores.dropna(subset=["pair_id"]).groupby("pair_id"):
        if len(rows) != 2 or set(rows["condition"]) != set(CONDITIONS[:2]):
            raise ValueError(f"Invalid recipient pair for {title}")
        ax.plot(rows["sample_id"].map(positions), rows["score"], color=PAIR, linewidth=1.45, zorder=2)

    labels = [f"{condition}\n(n={scores['condition'].eq(condition).sum()})"
              for condition in conditions]
    ax.set_xticks(range(len(conditions)), labels, fontsize=10, color=MUTED)
    for tick, condition in zip(ax.get_xticklabels(), conditions):
        tick.set_color(COLORS[condition])
    ax.set_xlim(-0.42, len(conditions) - 0.58)
    style_axis(ax, axis_label)
    ax.text(0, 1.075, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", color=MUTED)
    ax.annotate(title, xy=(0, 1.075), xycoords="axes fraction", xytext=(17, 0),
                textcoords="offset points", fontsize=12, fontweight="bold", color=INK)


def _bracket_label(q):
    if not np.isfinite(q):
        return "NA"
    stars = "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"
    return f"{stars}   q={q:.4g}"


def draw_comparison_brackets(ax, effects, conditions):
    """Place shorter comparisons below longer ones, above the observed scores."""
    bottom, top = ax.get_ylim()
    data_ticks = [value for value in ax.get_yticks() if bottom <= value <= top]
    height_fraction = 0.72 if len(effects) > 1 else 0.80
    ax.set_ylim(bottom, bottom + (top - bottom) / height_fraction)
    transform = ax.get_xaxis_transform()
    for i, row in enumerate(effects.itertuples(index=False)):
        x_reference = conditions.index(row.reference)
        x_treatment = conditions.index(row.treatment)
        y = (0.74 + 0.135 * i) if len(effects) > 1 else 0.84
        color = COLORS[row.treatment]
        ax.plot([x_reference, x_reference, x_treatment, x_treatment], [y, y + 0.026, y + 0.026, y],
                color=color, linewidth=1.4, transform=transform, clip_on=False, zorder=5)
        ax.text((x_reference + x_treatment) / 2, y + 0.040, _bracket_label(row.q_value),
                color=INK, fontsize=10, fontweight="semibold", ha="center", va="bottom",
                transform=transform, zorder=6)
    # Keep the bracket area visually clear without changing the data scale.
    ax.set_yticks(data_ticks)


def plot_dataset_significance(scores, coverage, effects, dataset):
    """Make the publication figure with comparisons directly above the score columns."""
    mouse = dataset == "GSE228927"
    conditions = CONDITIONS if mouse else CONDITIONS[:2]
    planned = scores.loc[scores["planned_sample"] & scores["dataset"].eq(dataset) & scores["recipient"].eq("astrocyte")]
    effects = effects.loc[effects["dataset"].eq(dataset) & effects["recipient"].eq("astrocyte")]
    fig, axes = plt.subplots(1, len(ASTROCYTE_PROGRAMS), figsize=(10.0, 5.8))
    fig.subplots_adjust(left=0.108, right=0.966, bottom=0.252, top=0.771, wspace=0.38)
    for i, (ax, program) in enumerate(zip(axes, ASTROCYTE_PROGRAMS)):
        subset = planned.loc[planned["program"].eq(program)]
        panel_effects = effects.loc[effects["program"].eq(program)]
        title = program
        if panel_effects["low_coverage"].any():
            title += " †"
        axis_label = "Module score" if program == "ACHES" else "Program score"
        draw_panel(ax, subset, "AB"[i], title, conditions, axis_label)
        draw_comparison_brackets(ax, panel_effects, conditions)

    species = "Mouse astrocytes" if mouse else "Human iPSC-derived astrocytes"
    scale = "log₂(FPKM + 1)" if mouse else "log₂(TPM + 1)"
    fig.text(0.065, 0.935, f"{dataset}  /  {species}", fontsize=21, fontweight="bold", color=INK)
    fig.text(0.065, 0.883, f"Microglia-conditioned medium (MCM) exposure  ·  {scale}",
             fontsize=11, color=MUTED)
    fig.add_artist(Line2D([0.065, 0.966], [0.855, 0.855], transform=fig.transFigure,
                          linewidth=0.8, color=RULE))
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=6, color=COLORS["Resting-MCM"],
               label="Individual sample"),
        Line2D([], [], color=INK, linewidth=2, label="Condition mean"),
    ]
    if not mouse:
        handles.append(Line2D([], [], color=PAIR, linewidth=1.45, label="Matched recipient line"))
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.060, 0.149), ncol=len(handles),
               frameon=False, labelcolor=MUTED, fontsize=10, columnspacing=2)
    fig.text(0.065, 0.109, "BH-FDR:  * q < 0.05    ** q < 0.01    *** q < 0.001    ns q ≥ 0.05",
             fontsize=9.5, color=INK)
    if mouse:
        snap = coverage.loc[coverage["dataset"].eq(dataset) & coverage["program"].eq("SNAP-a")].iloc[0]
        fig.text(0.065, 0.067, "Welch tests assume independent groups. Across-condition pairing is unavailable.",
                 fontsize=9, color=MUTED)
        fig.text(0.065, 0.031, f"† Exploratory SNAP-a: {snap['gene_coverage']:.0%} of source genes "
                 f"({snap['weight_mass_coverage']:.0%} of absolute loading weight).",
                 fontsize=9, color=MUTED)
    else:
        fig.text(0.065, 0.065, "Two-sided paired t tests use within-line score differences.",
                 fontsize=9, color=MUTED)
    return fig


def save_program_score_figures(scores, coverage, results_dir, figures_dir):
    """Write effect tables and numbered figures to their respective directories."""
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    effects = summarize_score_effects(scores)
    effects.to_csv(results_dir / "program_score_effects_ci.tsv", sep="\t", index=False, float_format="%.12g")
    with plt.rc_context({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none",
                         "figure.facecolor": "white", "axes.facecolor": "white"}):
        for dataset, number in FIGURE_NUMBERS.items():
            figure = plot_dataset_significance(scores, coverage, effects, dataset)
            try:
                for suffix in ("pdf", "png", "svg"):
                    path = figures_dir / f"{number}.{suffix}"
                    figure.savefig(path, dpi=300, facecolor="white")
                    if suffix == "svg":
                        text = path.read_text(encoding="utf-8")
                        path.write_text("\n".join(line.rstrip() for line in text.splitlines()) + "\n", encoding="utf-8")
            finally:
                plt.close(figure)
