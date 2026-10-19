"""
calibration.py

Compares declared base quality (Phred Q) with the empirically observed error
rate at that Q, per technology, to assess whether vendors over- or
under-state quality.

"""
from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

MIN_ERRORS_FOR_PLOT = 30


@dataclass
class Calibration:
    technology: str
    n_bases: Counter      # base_qual -> count of counted bases
    n_errors: Counter     # base_qual -> count of counted errors

    def rows(self) -> list[tuple[int, int, int, float, float | None]]:
        """Per-Q ``(base_qual, n_bases, n_errors, error_rate, empirical_q)``.

        ``empirical_q`` is ``None`` (empty) when ``n_errors == 0``.
        """
        rows: list[tuple[int, int, int, float, float | None]] = []
        for q in sorted(self.n_bases):
            nb = self.n_bases[q]
            if nb == 0:
                continue
            ne = self.n_errors.get(q, 0)
            rate = ne / nb
            eq = -10 * math.log10(rate) if ne > 0 else None
            rows.append((q, nb, ne, rate, eq))
        return rows

    def write_csv(self, path: str | Path) -> None:
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["technology", "base_qual", "n_bases", "n_errors",
                             "error_rate", "empirical_q"])
            for q, nb, ne, rate, eq in self.rows():
                writer.writerow([self.technology, q, nb, ne, rate,
                                 "" if eq is None else eq])


def summarize_calibration(calibrations: dict[str, Calibration]) -> str:
    """One text table per technology with the observed/expected error-rate ratio."""
    lines = [f"{'technology':<12}{'n_bases':>10}{'n_errors':>10}"
             f"{'observed_error_rate':>20}{'expected_error_rate':>20}"
             f"{'observed/expected':>16}"]
    for tech, cal in calibrations.items():
        total_bases = sum(cal.n_bases.values())
        total_errors = sum(cal.n_errors.values())
        observed = total_errors / total_bases if total_bases else 0.0
        expected = (sum(nb * 10 ** (-q / 10) for q, nb in cal.n_bases.items())
                    / total_bases if total_bases else 0.0)
        ratio = observed / expected if expected else 0.0
        lines.append(
            f"{tech:<12}{total_bases:>10}{total_errors:>10}"
            f"{observed:>20.6f}{expected:>20.6f}{ratio:>16.6f}"
        )
    return "\n".join(lines)


def plot_calibration(calibrations: dict[str, Calibration], out_png: str | Path) -> None:
    """Declared Q (x) vs empirical Q (y), one series per technology, y=x reference."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    xs_all: list[int] = []
    ys_all: list[float] = []
    for tech, cal in calibrations.items():
        xs, ys = [], []
        for q, nb, ne, rate, eq in cal.rows():
            if ne >= MIN_ERRORS_FOR_PLOT and eq is not None:
                xs.append(q)
                ys.append(eq)
        if xs:
            ax.plot(xs, ys, marker="o", label=tech)
            xs_all.extend(xs)
            ys_all.extend(ys)
    if xs_all:
        lo = min(min(xs_all), min(ys_all)) - 1
        hi = max(max(xs_all), max(ys_all)) + 1
    else:
        lo, hi = 0, 40
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="gray", label="y = x")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("declared base quality (Q)")
    ax.set_ylabel("empirical base quality (Q)")
    ax.set_title("Base-quality calibration: declared vs empirical Q")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out_png))
    plt.close(fig)
