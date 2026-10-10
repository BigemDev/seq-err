"""
Modules:
    run_all             - map paired reads and run the comparison workflow
    bam_reader         - iterate alignments, extract per-base mismatches from CIGAR/MD
    vcf_reader          - load known variant positions from a VCF
    mismatch_classify   - split mismatches into "variant" vs "sequencing error"
    bqsr_compare        - compare base qualities at the same positions before/after BQSR
    quality_stats       - build quality-score distributions and compare technologies
    pipeline            - individual pipeline steps and CLI entry point

"""

from .bam_reader import Mismatch, iter_mismatches
from .vcf_reader import load_variant_positions
from .mismatch_classify import classify_mismatches
from .bqsr_compare import compare_bqsr
from .quality_stats import QualityDistribution, compare_technologies
from .run_all import map_paired_bwa_mem, run_all

__all__ = [
    "map_paired_bwa_mem",
    "Mismatch",
    "iter_mismatches",
    "load_variant_positions",
    "classify_mismatches",
    "compare_bqsr",
    "QualityDistribution",
    "compare_technologies",
    "run_all",
]
