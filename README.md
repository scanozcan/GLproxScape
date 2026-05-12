# GLproxScape

<!-- badges: start -->
[![R-CMD-check](https://github.com/scanozcan/GLproxScape/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scanozcan/GLproxScape/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

Spatial deconvolution of **dCas9-APEX2 proximity proteomics** into
transcription-factor and chromatin-regulator binding predictions.

GLproxScape models each guide RNA's biotinylation footprint as a Gaussian
labelling cone (default σ = 300 bp), forward-smears the per-region
proteomics enrichment into a continuous spatial track *s(x)*, normalises
by guide coverage *C(x)* to recover an occupancy estimate β(x) = s(x) /
max(C(x), c_min · max C), then deconvolves into TF binding events on a
JASPAR position-weight-matrix basis. A separate zone-based path handles
chromatin readers / writers / erasers / remodellers that lack a
sequence-specific motif, and an optional ChIP-Atlas overlay validates
predictions against independent ChIP-seq peaks.

The package is the analysis backbone behind the GLproxScape preprint
(Ozcan *et al.*, in preparation) and ships bundled reanalyses of the
Myers 2018 hTERT/MYC dataset (5 sgRNAs each) and the Mackenzie 2026
FOXP2 dataset (3 sgRNAs).

## Install

```r
# install.packages("remotes")
remotes::install_github("scanozcan/GLproxScape", build_vignettes = TRUE)
```

`build_vignettes = TRUE` builds the FOXP2 walkthrough during install so
`vignette("foxp2-mackenzie", package = "GLproxScape")` works. Omit it
for a faster install if you don't need the vignette HTML.

While the repo is private (paper review period), you'll need a GitHub
Personal Access Token with `repo` scope in `~/.Renviron` as
`GITHUB_PAT=ghp_...` before `remotes::install_github()` can authenticate.

## Quick start

```r
library(GLproxScape)

# Use the bundled FOXP2 (Mackenzie 2026) example data
inputs_dir <- system.file("extdata/examples/foxp2_mackenzie",
                          package = "GLproxScape")
inputs <- load_caspex_inputs(inputs_dir)

res <- run_caspex(
  gene             = "FOXP2",
  transcript       = "ENST00000901759",   # Mackenzie's HEK293 active TSS1
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  upstream         = 200,
  downstream       = 2000,
  out_dir          = tempfile("foxp2_run_"),
  weight_mode      = "lfc_signed",
  pval_thresh      = 0.5,
  motif_thresh     = 0.75,
  chipatlas        = FALSE                # set TRUE for ChIP-Atlas overlay
)

# Inspect the top-N predicted binding events
head(res$binding_events[order(-res$binding_events$weight), ], 10)
```

The full deck (gRNA layout, per-region heatmap, motif-track PDF, and
per-TF deconvolution detail pages) is written to `out_dir`. The
returned list also exposes `res$spatial_df`, `res$motif_results`,
`res$promoter_info`, etc. for programmatic post-processing.

For the diagnostic plot pack (TF-pair co-occurrence triangle, sigma
sensitivity, jackknife stability, etc.):

```r
extras <- run_caspex_extras(res, out_dir = file.path(res$out_dir, "extras"))
```

For chromatin readers / writers / erasers that lack a sequence-specific
motif (e.g. BRD4, KMT2A, SMARCA4), call the zone-based path:

```r
extras_epi <- run_caspex_epigenetic(
  res,
  epigenetic_factors = readLines(system.file(
    "extdata/databases/EpiGenes_main.csv", package = "GLproxScape"))
)
```

## Bundled example datasets

Under `inst/extdata/examples/`, resolvable at runtime via `system.file()`:

| Folder                | Locus  | Guides | Reference                       |
|-----------------------|--------|--------|---------------------------------|
| `foxp2_mackenzie/`    | FOXP2  | 3      | Mackenzie *et al.*, 2026        |
| `tert_myers/`         | hTERT  | 5      | Myers *et al.*, 2018            |
| `myc_myers/`          | MYC    | 5      | Myers *et al.*, 2018            |

Each folder is a self-contained input bundle: a `grnas.tsv` manifest
(region → protospacer + per-region file) plus per-region `Region*.txt`
proteomics tables with `logFC` + `P.Value` columns.

A fourth bundled resource under `inst/extdata/databases/` contains the
TF and chromatin-factor universes (`TFLibrary.txt`, `EpiGenes_main.csv`,
`EpiGenes_complexes.csv`) sourced from public TF / EpiFactors databases.

## Vignette

End-to-end walkthrough of the FOXP2 reanalysis, ready to drop into a
methods-paper supplementary materials:

```r
vignette("foxp2-mackenzie", package = "GLproxScape")
```

It covers: loading the bundled inputs, picking the right transcript
anchor for an alt-promoter dataset, running the full pipeline, inspecting
the binding-events table, comparing PWM-score weighting modes
(`"none"` / `"linear"` / `"log"`), and the optional ChIP-Atlas overlay.

## Workflow summary

The pipeline runs in one call (`run_caspex`) that internally chains:

1. **Gene lookup** — resolves the canonical (or explicitly pinned)
   Ensembl transcript and fetches the promoter sequence over
   `[-upstream, +downstream]` bp of the TSS.
2. **gRNA matching** — exact-match each protospacer (and its reverse
   complement) against the promoter sequence to anchor the labelling
   kernel.
3. **Spatial model** — per-protein per-region significance gating
   (`pval_thresh`, `min_regions`), then `compute_spatial` aggregates
   into a TF-level summary scoring composite / specificity.
4. **Motif scan** — for each TF in the deck roster, fetches the JASPAR
   position weight matrix and scans the promoter at
   `motif_thresh × max_score`.
5. **Coverage-aware deconvolution** — `β(x) = s(x) / max(C(x), cov_floor ·
   max(C))` thresholded into zones, with one event emitted per JASPAR hit
   inside each zone. Optional PWM-score weighting reshapes per-event
   amplitudes (`motif_score_weight = "none" | "linear" | "log"`).
6. **Edge guards + merge** — events outside the gRNA support cone
   (`edge_guard_frac`, `max_grna_distance`) are dropped; closely-spaced
   events within `merge_dist` are merged into amplitude-weighted
   clusters.
7. **(Optional) ChIP-Atlas overlay** — public ChIP-seq peaks for each
   deck TF rendered beneath the predicted event bubbles for independent
   validation. Per-SRX BEDs cached under `tools::R_user_dir("caspex",
   "cache")`.
8. **Plotting + CSV write-out** — multi-page deconvolution PDF, motif
   track, gRNA-positions plot, per-TF binding-events CSV, spatial
   predictions CSV.

The supplemental methods in the accompanying paper give the precise
mathematical framing for every step.

## Reproducibility

Every analysis the paper reports — Mackenzie FOXP2, Myers hTERT and MYC
— is reproducible from the bundled `inst/extdata/examples/` folders.
The relevant `run_caspex()` parameter settings live in the FOXP2
vignette and in each dataset's runner script under the paper's
companion analysis folder.

The `transcript = "canonical"` default in `lookup_gene()` ensures
TSS-relative coordinates are stable across Ensembl releases. For
datasets that target alternative promoters (e.g. Mackenzie FOXP2 uses
HEK293's "active TSS1", not the canonical FOXP2-201 TSS), the runner
explicitly pins the right ENST ID (e.g. `ENST00000901759`); the package
exposes `check_transcripts()` as a diagnostic helper to identify
the correct anchor.

## Citation

If you use GLproxScape in your work, please cite:

> Ozcan C., *et al.* **GLproxScape: spatial deconvolution of
> dCas9-APEX2 proximity proteomics into TF and chromatin-regulator
> binding predictions.** (in preparation, 2026)

A machine-readable `inst/CITATION` will be added when the paper is
accepted; until then, `citation("GLproxScape")` returns a generic entry
auto-derived from `DESCRIPTION`.

## License

[MIT](LICENSE). See [LICENSE.md](LICENSE.md) for the full text.

## Issues and contributions

While the repo is private during paper review, contributions are
welcome from invited collaborators. File issues at
[github.com/scanozcan/GLproxScape/issues](https://github.com/scanozcan/GLproxScape/issues).
After the paper is published the repo will be made public and the
contributing workflow will switch to standard fork-and-PR.
