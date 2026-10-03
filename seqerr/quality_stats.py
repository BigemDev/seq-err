"""
quality_stats.py

Computes and compares base-quality distributions of sequencing errors across technologies.

"""
from __future__ import annotations

import csv
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .bam_reader import Mismatch


@dataclass
class QualityDistribution:
    technology: str
    n: int
    mean_qual: float
    median_qual: float
    stdev_qual: float
    histogram: Counter  # base_qual -> count

    @classmethod
    def from_mismatches(cls, technology: str, mismatches: Iterable[Mismatch]) -> "QualityDistribution":
        quals = [m.base_qual for m in mismatches if m.technology == technology]
        if not quals:
            return cls(technology, 0, 0.0, 0.0, 0.0, Counter())
        return cls(
            technology=technology,
            n=len(quals),
            mean_qual=statistics.fmean(quals),
            median_qual=statistics.median(quals),
            stdev_qual=statistics.pstdev(quals) if len(quals) > 1 else 0.0,
            histogram=Counter(quals),
        )

    @classmethod
    def from_histogram(cls, technology: str, histogram: Counter) -> "QualityDistribution":
        hist = Counter({q: c for q, c in histogram.items() if c > 0})
        n = sum(hist.values())
        if n == 0:
            return cls(technology, 0, 0.0, 0.0, 0.0, Counter())
        mean = sum(q * c for q, c in hist.items()) / n
        ordered = sorted(hist)
        def kth(k):
            seen = 0
            for q in ordered:
                seen += hist[q]
                if k < seen:
                    return q
        median = kth(n // 2) if n % 2 else (kth(n // 2 - 1) + kth(n // 2)) / 2
        var = sum(c * (q - mean) ** 2 for q, c in hist.items()) / n
        return cls(technology, n, mean, float(median), var ** 0.5 if n > 1 else 0.0, hist)

    def write_histogram_csv(self, out_path: str | Path) -> None:
        with open(out_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["technology", "base_qual", "count"])
            for q in sorted(self.histogram):
                writer.writerow([self.technology, q, self.histogram[q]])


def compare_technologies(
    mismatches: Iterable[Mismatch],
    technologies: list[str],
) -> dict[str, QualityDistribution]:
    mismatches = list(mismatches)
    return {
        tech: QualityDistribution.from_mismatches(tech, mismatches)
        for tech in technologies
    }


def summarize_comparison(distributions: dict[str, QualityDistribution]) -> str:
    lines = [f"{'technology':<12}{'n_errors':>10}{'mean_Q':>10}{'median_Q':>10}{'stdev_Q':>10}"]
    for tech, d in distributions.items():
        lines.append(
            f"{tech:<12}{d.n:>10}{d.mean_qual:>10.2f}{d.median_qual:>10.2f}{d.stdev_qual:>10.2f}"
        )
    return "\n".join(lines)