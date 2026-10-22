# Micropatterned soma, proximal neurite, and distal neurite RNA profiles

Code, data, and photomask design for:

> Kim J, Viner C, Rosario JG, Jeong H, Wu J, Lee D, Kim J. A micropatterning platform for profiling soma, proximal neurite, and distal neurite RNA.

Mouse E18 hippocampal neurons were cultured on micropatterned coverslips and microdissected at 16 days in vitro into soma, proximal neurite, and distal neurite fractions, alongside whole neurons from micropatterned and standard cultures. Twelve matched sample sets from two litter pools gave 59 SMART-Seq v4 / Nextera XT libraries sequenced on a NovaSeq 6000.

## Contents

- `data/Patterned_Neuron_rmdup_exon_counts.tsv`: deduplicated exon counts, 23,545 genes (including 92 ERCC reference rows; no spike-in was added) by 59 libraries. This is the GEO processed file and the input to every analysis.
- `data/Patterned_Neuron_counted_exon_coordinates_mm10.tsv`: mm10 exon coordinates for every counted gene, from the annotation used for counting.
- `data/published/`: Supplementary Files 1–4 and Supplementary Table 1, used to check the code against the published tables.
- `data/upstream/`: per-library read-pair counts before and after deduplication, and non-deduplicated counts for five libraries.
- `photomask/`: the micropattern GDSII design (Supplementary File 5).
- `python/`, `R/`: one script per analysis.
- `environment/`: Python requirements, the R package installer, and session records.
- `results/`: everything the scripts write, committed so readers can compare without rerunning.
- `figures/`: Figure 2 and Supplementary Figures 1–3.
- `claims.tsv`: every number in the paper with the script output that produces it. `python/check_claims.py` compares them and writes `results/claims_check.tsv`.

## Requirements

Python 3.14 with `environment/requirements.txt`. R 4.5.3 with Bioconductor 3.21; `environment/install_R_packages.R` installs the versions used (DESeq2 1.48.2, apeglm 1.30.0, gprofiler2 0.2.4, org.Mm.eg.db 3.21.0, GOSemSim, clusterProfiler, ggplot2, patchwork, ggrepel, pheatmap).

## Running

    PYTHON=python3.14 R_CMD=R ./run_all.sh

`python/fetch_external.py` downloads the Middleton et al. (2019) supplementary tables from the journal and checks their SHA-256. `R/middleton_reanalysis.R` downloads that study's counts from GEO. The full run takes about an hour, most of it DESeq2 and GSEA.

## Analyses

| Paper section | Script |
|---|---|
| Micropattern versus standard: PERMANOVA, classifier, paired Hotelling-like test, 8-of-9 consistency, bootstrap, subsampling, stress markers, PCA covariates, batch | `python/micropattern.py` |
| Sample correlations, outlier screen, library sizes | `python/qc.py` |
| Deduplication | `python/upstream.py` |
| Proximal versus distal: PERMANOVA and leave-one-pair-out, 10-of-12 consistency, polarity, axon and glial markers, PERMDISP, correlations, culture-12 omission | `python/proximal_distal.py` |
| Nulls for the 54-gene candidate set | `python/proximal_distal_nulls.py` |
| Depth matching and pseudo-distal controls | `python/depth_controls.py` |
| Normalization dependence | `python/normalization_robustness.py` |
| Neurite versus soma overlaps with Middleton et al. (2019); consistently detected genes; stress-marker order | `python/neurite_soma.py` |
| DESeq2 (Supplementary Files 1, 2, 4) | `R/deseq2.R` |
| Cook's distance note | `R/cooks.R` |
| DESeq2 design check | `R/design_check.R` |
| g:Profiler ORA, GSEA, Middleton-style enrichGO | `R/go_enrichment.R` |
| GOSemSim | `R/gosemsim.R` |
| Supplementary Note 1 (Middleton et al. reanalysis) | `R/middleton_reanalysis.R` |
| Supplementary Figures 1–3 | `R/supplementary_figures.R` |
| Figure 2 | `python/figure2.py` (panel data), `R/figure2.R` (figure) |

Permutation tests use fixed seeds. When several tests share one generator, the earlier draws are consumed first.

## License

Copyright (C) 2026 Jihoon Kim, Coby Viner, Jean G. Rosario, Hyejoong Jeong, Jingyu Wu, Daeyeon Lee, and Junhyong Kim. GNU General Public License, version 3 (`LICENSE`).
