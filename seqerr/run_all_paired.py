"""Compatibility entry point for the paired-read pipeline."""

from .run_all import (
    TechnologyResult,
    _extract_and_classify,
    _intersect_vcfs,
    _prepare_reference,
    _write_summary,
    add_threading_args,
    build_parser,
    ensure_bwa_index,
    main,
    map_paired_bwa_mem,
    plan_intervals,
    read_fai,
    run_all,
    run_haplotypecaller_parallel,
)

run_all_paired = run_all


if __name__ == "__main__":
    main()
