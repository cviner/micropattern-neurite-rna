#!/usr/bin/env bash
# Reproduce every result and figure from data/, then compare them with the numbers in the paper.
set -euo pipefail
cd "$(dirname "$0")"

PY=${PYTHON:-python3}
R_CMD=${R_CMD:-R}
run_r() { "$R_CMD" --vanilla --quiet -f "$1"; }

$PY python/fetch_external.py

$PY python/micropattern.py
$PY python/qc.py
$PY python/upstream.py
$PY python/proximal_distal_nulls.py
$PY python/depth_controls.py
$PY python/normalization_robustness.py
$PY python/proximal_distal.py
$PY python/neurite_soma.py
$PY python/figure2.py

run_r R/deseq2.R
run_r R/cooks.R
run_r R/design_check.R
run_r R/go_enrichment.R
run_r R/gosemsim.R
run_r R/middleton_reanalysis.R
run_r R/supplementary_figures.R
run_r R/figure2.R

$PY python/check_claims.py
