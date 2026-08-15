#!/usr/bin/env python3
"""Generate the theorem-faithful vector figures used by the manuscript."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from decimal import Decimal
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "forr4_n4096_matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyArrowPatch, Rectangle


PROJECT = Path(__file__).resolve().parents[1]
FIGURES = PROJECT / "figures"
ROUND4 = (
    PROJECT
    / "code"
    / "source_snapshot"
    / "open_problem_forr4_passive_floor_consolidation_round_4"
)
LEDGER = ROUND4 / "artifacts" / "q64_complete_outward_ledger.json"

INK = "#1f2933"
BLUE = "#315d7d"
BLUE_DARK = "#173f5f"
OCHRE = "#d8a64d"
CREAM = "#f7f3e8"
GRAY = "#d9dee3"
RED = "#9a4f43"
GREEN = "#557a63"
PHASE_CMAP = ListedColormap([CREAM, BLUE_DARK])


def configure() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.0,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK,
            "axes.titlecolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.transparent": False,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(
        FIGURES / f"{stem}.pdf",
        bbox_inches="tight",
        pad_inches=0.03,
        metadata={"Creator": "figures/generate_figures.py", "CreationDate": None},
    )
    svg_path = FIGURES / f"{stem}.svg"
    fig.savefig(
        svg_path,
        bbox_inches="tight",
        pad_inches=0.03,
        metadata={"Creator": "figures/generate_figures.py", "Date": None},
    )
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    plt.close(fig)


def sylvester_signs(order: int) -> np.ndarray:
    matrix = np.ones((1, 1), dtype=np.int8)
    while matrix.shape[0] < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    if matrix.shape != (order, order):
        raise ValueError("order must be a power of two")
    return matrix


def signed_pair(signs: np.ndarray, permutation: np.ndarray, sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q = signs.shape[0]
    left = np.empty((q, q), dtype=np.int8)
    right = np.empty((q, q), dtype=np.int8)
    for b in range(q):
        left[:, b] = signs[:, permutation[b]] * sigma[b]
        right[permutation[b], :] = sigma[b] * signs[b, :]
    return left, right


def signed_permutation_plant(q: int = 64, seed: int = 4096) -> list[np.ndarray]:
    signs = sylvester_signs(q)
    rng = np.random.default_rng(seed)
    pairs = []
    for _ in range(3):
        permutation = rng.permutation(q)
        sigma = rng.choice(np.array([-1, 1], dtype=np.int8), size=q)
        pairs.append(signed_pair(signs, permutation, sigma))
    return [
        pairs[0][0],
        pairs[0][1] * pairs[1][0],
        pairs[1][1] * pairs[2][0],
        pairs[2][1],
    ]


def fwht_integer(vector: np.ndarray) -> np.ndarray:
    out = np.asarray(vector, dtype=np.int64).copy()
    width = 1
    while width < out.size:
        for start in range(0, out.size, 2 * width):
            left = out[start : start + width].copy()
            right = out[start + width : start + 2 * width].copy()
            out[start : start + width] = left + right
            out[start + width : start + 2 * width] = left - right
        width *= 2
    return out


def exact_plant_numerator(plant: list[np.ndarray]) -> tuple[int, int]:
    q = plant[0].shape[0]
    vector = plant[3].reshape(-1).astype(np.int64)
    for block in (plant[2], plant[1], plant[0]):
        vector = fwht_integer(vector)
        vector *= block.reshape(-1)
    numerator = int(vector.sum())
    denominator = q**5
    return numerator, denominator


def phase_values(block: np.ndarray) -> np.ndarray:
    return (block < 0).astype(np.int8)


def draw_phase_icon(ax: plt.Axes, x: float, y: float, size: float, block: np.ndarray, label: str) -> None:
    values = phase_values(block)
    q = values.shape[0]
    cell = size / q
    for row in range(q):
        for column in range(q):
            color = BLUE_DARK if values[row, column] else CREAM
            ax.add_patch(
                Rectangle(
                    (x + column * cell, y + (q - 1 - row) * cell),
                    cell,
                    cell,
                    facecolor=color,
                    edgecolor="none",
                )
            )
    ax.add_patch(Rectangle((x, y), size, size, fill=False, edgecolor=INK, linewidth=0.7))
    ax.text(x + size / 2, y - 0.15, label, ha="center", va="top", fontsize=8)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], **kwargs) -> None:
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=0.9,
        color=kwargs.pop("color", INK),
        **kwargs,
    )
    ax.add_patch(patch)


def active_protocol_figure(plant: list[np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(7.15, 2.45))
    ax.set_xlim(0, 13.0)
    ax.set_ylim(0, 4.25)
    ax.axis("off")

    ax.text(0.05, 2.13, r"$|+\rangle_{\rm path}|u\rangle$", ha="left", va="center")
    ax.plot([1.55, 1.55], [1.08, 3.17], color=INK, linewidth=0.9)
    arrow(ax, (1.15, 2.13), (1.52, 2.13))
    arrow(ax, (1.55, 2.13), (1.92, 3.17))
    arrow(ax, (1.55, 2.13), (1.92, 1.08))
    ax.text(1.48, 2.38, "path\nsplitter", ha="right", va="bottom", fontsize=6.5)

    icon_size = 0.67
    draw_phase_icon(ax, 2.00, 2.83, icon_size, plant[0][::10, ::10][:6, :6], r"$D_1$")
    draw_phase_icon(ax, 5.05, 2.83, icon_size, plant[1][::10, ::10][:6, :6], r"$D_2$")
    draw_phase_icon(ax, 2.00, 0.74, icon_size, plant[3][::10, ::10][:6, :6], r"$D_4$")
    draw_phase_icon(ax, 5.05, 0.74, icon_size, plant[2][::10, ::10][:6, :6], r"$D_3$")

    for y in (3.17, 1.08):
        arrow(ax, (2.70, y), (3.48, y))
        ax.add_patch(Rectangle((3.50, y - 0.27), 0.68, 0.54, facecolor=GRAY, edgecolor=INK, linewidth=0.7))
        ax.text(3.84, y, r"$H_{4096}$", ha="center", va="center", fontsize=7)
        arrow(ax, (4.20, y), (5.02, y))

    arrow(ax, (5.75, 3.17), (7.18, 3.17))
    arrow(ax, (5.75, 1.08), (6.28, 1.08))
    ax.add_patch(Rectangle((6.30, 0.81), 0.68, 0.54, facecolor=GRAY, edgecolor=INK, linewidth=0.7))
    ax.text(6.64, 1.08, r"$H_{4096}$", ha="center", va="center", fontsize=7)
    arrow(ax, (7.00, 1.08), (7.18, 1.08))

    ax.text(6.15, 3.48, r"$|L_x\rangle$", ha="center", va="bottom")
    ax.text(6.15, 0.77, r"$|R_x\rangle$", ha="center", va="top")
    ax.text(3.85, 4.05, "two charged phase-grid crossings", ha="center", va="top", color=RED, fontsize=7)
    ax.text(3.85, 0.03, "two charged phase-grid crossings", ha="center", va="bottom", color=RED, fontsize=7)

    ax.plot([7.55, 7.55], [1.08, 3.17], color=INK, linewidth=0.9)
    arrow(ax, (7.18, 3.17), (7.52, 2.25))
    arrow(ax, (7.18, 1.08), (7.52, 2.02))
    ax.add_patch(Rectangle((7.78, 1.72), 1.02, 0.82, facecolor="white", edgecolor=INK, linewidth=0.8))
    ax.text(8.29, 2.13, "path-$X$\nreceiver", ha="center", va="center", fontsize=7)
    arrow(ax, (7.58, 2.13), (7.76, 2.13))
    arrow(ax, (8.82, 2.13), (9.50, 2.13))
    ax.text(9.02, 2.44, r"$X\in\{\pm1\}$", ha="center", va="bottom")
    ax.text(8.29, 1.52, "mode summed", ha="center", va="top", fontsize=6.5)

    ax.add_patch(Rectangle((9.55, 1.55), 1.50, 1.16, facecolor="#eef3f6", edgecolor=BLUE, linewidth=0.8))
    ax.text(10.30, 2.13, "repeat three\nindependent\nflags", ha="center", va="center", fontsize=6.3)
    arrow(ax, (11.07, 2.13), (11.42, 2.13))
    ax.add_patch(Rectangle((11.45, 1.65), 1.38, 0.96, facecolor="#f4eee8", edgecolor=RED, linewidth=0.8))
    ax.text(12.14, 2.13, "majority\ndecision", ha="center", va="center", fontsize=7)

    ax.text(
        8.25,
        3.63,
        r"$\langle L_x|R_x\rangle=F_{4,H}(x)$",
        ha="center",
        va="center",
        fontsize=9,
    )
    ax.text(
        11.25,
        3.63,
        r"hard dose $=3\times2=6$",
        ha="center",
        va="center",
        fontsize=8,
    )
    save_figure(fig, "active_protocol")


def phase_plant_figure(plant: list[np.ndarray]) -> None:
    numerator, denominator = exact_plant_numerator(plant)
    if numerator != denominator:
        raise AssertionError(f"plant identity failed: {numerator}/{denominator}")

    fig, axes = plt.subplots(1, 4, figsize=(7.15, 2.45))
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.25, top=0.80, wspace=0.30)
    for index, (axis, block) in enumerate(zip(axes, plant), start=1):
        values = phase_values(block)
        coordinates = np.arange(values.shape[0] + 1)
        axis.pcolormesh(
            coordinates,
            coordinates,
            values,
            cmap=PHASE_CMAP,
            vmin=0,
            vmax=1,
            shading="flat",
            rasterized=False,
        )
        axis.set_aspect("equal")
        axis.set_xticks([0, 32, 64])
        axis.set_yticks([0, 32, 64])
        axis.tick_params(length=2, labelsize=6)
        axis.set_title(rf"$z^{{({index})}}$", pad=3, fontsize=9)
        if index == 1:
            axis.set_ylabel("row coordinate $a$")
        axis.set_xlabel("column coordinate $b$")
        if index < 4:
            axis.text(
                1.08,
                0.5,
                rf"$P_{index}$",
                transform=axis.transAxes,
                ha="center",
                va="center",
                fontsize=7,
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=GRAY),
            )
    fig.text(0.5, 0.04, r"phase:  $0$ = light,  $\pi$ = dark", ha="center", va="bottom", fontsize=7)
    fig.text(
        0.5,
        0.965,
        r"one deterministic $q=64$ signed-permutation plant, verified $F_{4,H}(z)=1$",
        ha="center",
        va="top",
        fontsize=8,
    )
    save_figure(fig, "signed_permutation_phase_grids")


def rounded(value: Decimal, digits: int) -> str:
    return f"{value:.{digits}f}"


def certificate_figure() -> None:
    payload = json.loads(LEDGER.read_text(encoding="utf-8"))
    result = payload["result"]
    status_counts = dict(result["status_counts"])
    supported = int(result["supported_balanced_entries"])
    excluded = int(result["excluded_unbalanced_high_sector_entries"])
    records = int(result["high_sector_profile_splits"])
    if sum(status_counts.values()) != supported or supported + excluded != records:
        raise AssertionError("registry totals are inconsistent")

    perron = Decimal(str(result["collatz_perron_upper"]))
    promise = Decimal(str(result["promise_loss_upper"]))
    total = Decimal(str(result["total_upper"]))
    passive_error = (Decimal(1) - total) / Decimal(2)
    active_error = Decimal(81) / Decimal(256)
    threshold = Decimal(1) / Decimal(3)

    fig = plt.figure(figsize=(7.15, 4.65))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.35], hspace=0.48, wspace=0.42)
    flow = fig.add_subplot(grid[0, :])
    bars = fig.add_subplot(grid[1, 0])
    errors = fig.add_subplot(grid[1, 1])

    flow.set_xlim(0, 12)
    flow.set_ylim(0, 2.4)
    flow.axis("off")
    boxes = [
        (0.10, 0.75, 1.70, 0.95, "hard pair\n$q=64$, $\\beta=19/25$"),
        (2.08, 0.75, 1.70, 0.95, "888 certified\ncoefficients"),
        (4.06, 0.75, 1.70, 0.95, "$210\\times210$\nstate matrix"),
        (6.04, 0.75, 1.70, 0.95, f"Perron upper\n$\\leq{rounded(perron, 10)}$"),
        (8.02, 0.75, 1.70, 0.95, f"promise upper\n$\\leq{rounded(promise, 10)}$"),
        (10.00, 0.75, 1.70, 0.95, f"transcript TV\n$\\leq{rounded(total, 10)}$"),
    ]
    for index, (x, y, width, height, text_value) in enumerate(boxes):
        face = "#eef3f6" if index not in (4, 5) else "#f5efe7"
        edge = BLUE if index not in (4, 5) else OCHRE
        flow.add_patch(Rectangle((x, y), width, height, facecolor=face, edgecolor=edge, linewidth=0.9))
        flow.text(x + width / 2, y + height / 2, text_value, ha="center", va="center", fontsize=6.3)
        if index < len(boxes) - 1:
            arrow(flow, (x + width + 0.06, y + height / 2), (boxes[index + 1][0] - 0.06, y + height / 2))
    flow.text(
        6.89,
        0.42,
        "number-sector-incoherent adaptivity: frontier multiplier 1",
        ha="center",
        va="top",
        fontsize=7,
        color=GREEN,
    )
    flow.plot([6.04, 7.74], [0.55, 0.55], color=GREEN, linewidth=1.1)
    flow.text(0.02, 2.18, "a", fontsize=10, fontweight="bold", va="top")
    flow.text(0.35, 2.18, "certificate assembly", fontsize=8, va="top")

    grouped = [
        ("Inherited", status_counts["proved_nonuniversal_inherited"]),
        ("Local Walsh", status_counts["proved_masked_local_walsh"]),
        ("Residual chain", status_counts["proved_final_residual_chain"]),
        ("Quintic slice", status_counts["proved_masked_quintic_slice"]),
        ("Four-cubic", status_counts["proved_masked_four_cubic_incidence"]),
    ]
    grouped_total = sum(value for _, value in grouped)
    grouped.append(("Other certified", supported - grouped_total))
    names = [name for name, _ in grouped][::-1]
    values = [value for _, value in grouped][::-1]
    bars.barh(names, values, color=[GRAY, OCHRE, BLUE, "#6f8ea5", "#9bb0c0", BLUE_DARK], height=0.64)
    bars.set_xlim(0, 480)
    bars.set_xlabel("certified registry entries")
    bars.spines[["top", "right"]].set_visible(False)
    bars.tick_params(axis="y", length=0, labelsize=7)
    bars.grid(axis="x", color="#e7eaed", linewidth=0.6)
    bars.set_axisbelow(True)
    for index, value in enumerate(values):
        bars.text(value + 7, index, str(value), va="center", ha="left", fontsize=7)
    bars.set_title("b   registry composition", loc="left", fontsize=8, pad=8)

    errors.set_xlim(0, 0.42)
    errors.set_ylim(-0.65, 1.65)
    errors.set_yticks([0, 1], ["active protocol", "passive obstruction"])
    errors.set_xlabel("decision error probability")
    errors.spines[["left", "top", "right"]].set_visible(False)
    errors.tick_params(axis="y", length=0, labelsize=7)
    errors.grid(axis="x", color="#e7eaed", linewidth=0.6)
    errors.set_axisbelow(True)
    errors.axvline(float(threshold), color=INK, linestyle="--", linewidth=1.0)
    errors.text(float(threshold), 1.55, r"required $1/3$", ha="center", va="bottom", fontsize=7)
    errors.scatter([float(active_error)], [0], s=45, color=BLUE_DARK, zorder=3)
    errors.scatter([float(passive_error)], [1], s=45, color=RED, zorder=3)
    errors.arrow(
        float(passive_error),
        1,
        0.045,
        0,
        width=0.006,
        head_width=0.10,
        head_length=0.012,
        color=RED,
        length_includes_head=True,
        zorder=2,
    )
    errors.text(float(active_error) - 0.005, -0.20, r"$81/256=0.31640625$", ha="right", va="top", fontsize=7)
    errors.text(
        float(passive_error),
        1.18,
        rf"$\geq{rounded(passive_error, 10)}$",
        ha="center",
        va="bottom",
        fontsize=7,
    )
    errors.set_title("c   certified decision errors", loc="left", fontsize=8, pad=8)

    fig.text(
        0.5,
        0.01,
        "Registry bar lengths count proof entries; they do not measure numerical contribution.",
        ha="center",
        va="bottom",
        fontsize=6.5,
        color=INK,
    )
    save_figure(fig, "passive_certificate")


def check_outputs() -> None:
    plant = signed_permutation_plant()
    numerator, denominator = exact_plant_numerator(plant)
    if numerator != denominator:
        raise AssertionError("the displayed q=64 plant is not exact")
    payload = json.loads(LEDGER.read_text(encoding="utf-8"))
    result = payload["result"]
    if int(result["supported_balanced_entries"]) != 888:
        raise AssertionError("ledger registry count changed")
    if (
        int(result["excluded_unbalanced_high_sector_incidence_records"]),
        int(result["excluded_unbalanced_high_sector_undirected_edges"]),
    ) != (272, 136):
        raise AssertionError("number-sector audit changed")
    for stem in ("active_protocol", "signed_permutation_phase_grids", "passive_certificate"):
        for suffix in ("pdf", "svg"):
            path = FIGURES / f"{stem}.{suffix}"
            if not path.is_file() or path.stat().st_size == 0:
                raise AssertionError(f"missing figure output: {path.name}")
    certificate_svg = (FIGURES / "passive_certificate.svg").read_text(encoding="utf-8")
    for required in (
        "0.2587440964",
        "0.2609692248",
        "0.3695153876",
        "number-sector-incoherent adaptivity",
    ):
        if required not in certificate_svg:
            raise AssertionError(f"stale passive-certificate figure: {required}")
    print(
        "PASS figure contract: exact q=64 plant, 2D phase grids, "
        "888-entry balanced ledger, 272/136 scope audit, and six vector outputs"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="check committed outputs without rewriting them")
    args = parser.parse_args()
    configure()
    if args.check:
        check_outputs()
        return
    plant = signed_permutation_plant()
    active_protocol_figure(plant)
    phase_plant_figure(plant)
    certificate_figure()
    check_outputs()


if __name__ == "__main__":
    main()
