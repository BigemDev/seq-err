"""Calculate empirical mismatch error rates and save CSVs/figures.

Denominators contain ALL eligible base observations, including correct bases.
VCF loci and non-ACGT reference positions are excluded from both numerators
and denominators. For depth we store number of bases directly: this is valid
with min_base_qual > 0, where positions_at_depth * depth would overcount.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Mapping


def quality_rates(
    error_hist: Mapping[int, int], total_qual: Mapping[int, int]
) -> list[tuple[int, int, int, float, float | None]]:
    rows = []
    for q, n_bases in sorted(total_qual.items()):
        if n_bases <= 0:
            continue
        errors = int(error_hist.get(q, 0))
        if errors > n_bases:
            raise ValueError(f"Quality {q}: errors ({errors}) > eligible bases ({n_bases})")
        rate = errors / n_bases
        empirical_q = -10 * math.log10(rate) if rate > 0 else None
        rows.append((q, errors, int(n_bases), rate, empirical_q))
    return rows


def depth_rates(
    error_depth: Mapping[int, int], total_depth_bases: Mapping[int, int]
) -> list[tuple[int, int, int, float]]:
    rows = []
    for depth, n_bases in sorted(total_depth_bases.items()):
        if depth <= 0 or n_bases <= 0:
            continue
        errors = int(error_depth.get(depth, 0))
        if errors > n_bases:
            raise ValueError(f"Depth {depth}: errors ({errors}) > eligible bases ({n_bases})")
        rows.append((depth, errors, int(n_bases), errors / n_bases))
    return rows


def _write_csv(path: Path, header: list[str], rows: list[tuple]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for row in rows:
            writer.writerow(["" if val is None else val for val in row])


def plot_error_rates(
    label: str,
    error_hist: Mapping[int, int],
    total_qual: Mapping[int, int],
    error_depth: Mapping[int, int],
    total_depth_bases: Mapping[int, int],
    out_dir: str | Path,
) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    q_rows = quality_rates(error_hist, total_qual)
    d_rows = depth_rates(error_depth, total_depth_bases)

    qual_csv = out_dir / "error_rate_by_quality.csv"
    depth_csv = out_dir / "error_rate_by_depth.csv"
    _write_csv(qual_csv, ["base_qual", "errors", "total_bases", "error_rate", "empirical_q"], q_rows)
    _write_csv(depth_csv, ["depth", "errors", "total_bases", "error_rate"], d_rows)
    written = [qual_csv, depth_csv]

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[plots] Matplotlib unavailable: wrote CSV tables only")
        return written

    title = label.upper()
    pts = [(q, rate) for q, _, _, rate, _ in q_rows if rate > 0]
    if pts:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-")
        ax.set_yscale("log")
        ax.set_xlabel("Reported base quality (Phred)")
        ax.set_ylabel("Error rate (errors / eligible bases)")
        ax.set_title(f"{title}: error rate by base quality")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        path = out_dir / "error_rate_by_quality.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)

    pts = [(q, eq) for q, _, _, _, eq in q_rows if eq is not None]
    if pts:
        fig, ax = plt.subplots(figsize=(6, 6))
        top = max(max(q for q, _ in pts), max(eq for _, eq in pts)) + 2
        ax.plot([0, top], [0, top], "k--", alpha=0.5, label="ideal")
        ax.plot([q for q, _ in pts], [eq for _, eq in pts], "o", label=label)
        ax.set_xlim(0, top)
        ax.set_ylim(0, top)
        ax.set_xlabel("Reported quality (Phred)")
        ax.set_ylabel("Empirical quality (-10 log10(error rate))")
        ax.set_title(f"{title}: reported vs empirical quality")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        path = out_dir / "reported_vs_empirical_quality.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)

    pts = [(depth, rate) for depth, _, _, rate in d_rows]
    if pts:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-")
        ax.set_xlabel("Filtered depth")
        ax.set_ylabel("Error rate (errors / eligible bases)")
        ax.set_title(f"{title}: error rate by depth")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        path = out_dir / "error_rate_by_depth.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)

    return written

def plot_quality_depth_heatmap_from_counts(
    error_counts: Counter,
    base_counts: Counter,
    technology: str,
    output_dir: str = "histograms",
) -> Path | None:
    """Empirical mismatch rate for each (depth bin, base-quality bin).

    IMPORTANT: base_counts must contain ALL eligible base observations at
    callable, non-VCF sites, including correct observations. Do not pass
    the number of errors at each depth as the denominator.
    """
    if not base_counts:
        return None

    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

    output_path = Path(output_dir) / technology
    output_path.mkdir(parents=True, exist_ok=True)

    depth_bins = sorted({d for d, _ in base_counts})
    quality_bins = sorted({q for _, q in base_counts})
    col_index = {d: i for i, d in enumerate(depth_bins)}
    row_index = {q: i for i, q in enumerate(quality_bins)}

    # NaN means "no eligible bases". 0.0 means "bases, but zero errors".
    matrix = np.full((len(quality_bins), len(depth_bins)), np.nan, dtype=float)
    rows = []
    for (depth_bin, quality_bin), n_bases in sorted(base_counts.items()):
        if n_bases <= 0:
            continue
        n_errors = int(error_counts.get((depth_bin, quality_bin), 0))
        if n_errors > n_bases:
            raise ValueError(
                f"({depth_bin}, {quality_bin}): {n_errors} errors > {n_bases} bases"
            )
        rate = n_errors / n_bases
        matrix[row_index[quality_bin], col_index[depth_bin]] = rate
        rows.append((depth_bin, quality_bin, n_errors, n_bases, rate))

    # The CSV carries denominators, so a cell based on 5 bases cannot be
    # confused with one based on 5 million bases.
    csv_path = output_path / "error_rate_by_quality_depth.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["depth_bin", "base_quality_bin", "errors", "total_bases", "error_rate"])
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(11, 6))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgray")

    # Fixed color scale is comparable across technologies; the symlog
    # transform distinguishes small error rates without discarding zero.
    im = ax.imshow(
        np.ma.masked_invalid(matrix), aspect="auto", origin="lower",
        cmap=cmap, norm=SymLogNorm(linthresh=1e-5, vmin=0, vmax=1, base=10),
        interpolation="nearest",
    )

    # Avoid unreadable labels on large / varied-depth datasets.
    x_step = max(1, (len(depth_bins) + 15) // 16)
    y_step = max(1, (len(quality_bins) + 11) // 12)
    x_ticks = list(range(0, len(depth_bins), x_step))
    y_ticks = list(range(0, len(quality_bins), y_step))
    ax.set_xticks(x_ticks, [str(depth_bins[i]) for i in x_ticks], rotation=45, ha="right")
    ax.set_yticks(y_ticks, [str(quality_bins[i]) for i in y_ticks])
    ax.set_xlabel("Filtered depth bin (lower bound, step 5)")
    ax.set_ylabel("Base quality bin (lower bound, step 5)")
    ax.set_title(f"{technology.upper()}: mismatch rate by depth and base quality")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Observed mismatch rate (errors / eligible bases)")
    ax.text(0, -0.19, "Gray: no eligible bases  |  Zero rate: no observed mismatches",
            transform=ax.transAxes, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path / "base_quality_vs_depth_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path / "base_quality_vs_depth_heatmap.png"

