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

## How to use it on your own data

GLproxScape expects each promoter-tiling experiment to live in a
folder with one **manifest** file plus one **proteomics table** per
gRNA region. Once the folder is laid out, the entire analysis is two
function calls.

### Folder layout

```
my_gene_analysis/
├── inputs/
│   ├── grnas.tsv        ← manifest: which gRNA is which region, and where its data lives
│   ├── Region1.txt      ← per-region proteomics table (one per sgRNA)
│   ├── Region2.txt
│   ├── Region3.txt
│   └── ...
└── analysis.R           ← your runner script (whatever you want to call it)
```

The folder name is up to you. What matters is that `inputs/grnas.tsv`
points at the `Region*.txt` files that sit next to it.

### `grnas.tsv` — the manifest

A tab-separated file with three required columns: `region`, `sequence`,
`data_file`. One row per gRNA region. Lines starting with `#` and blank
lines are ignored, so the top of the file can carry provenance notes.

```
# YOUR_GENE promoter tiling — replace with your sgRNA sequences.
# Lines beginning with `#` are ignored; useful for provenance notes
# (assembly, chromosome, strand, etc.).
region	sequence	data_file
R1	GAAATCCAGGAGTCATATAA	Region1.txt
R2	GGCAGTCAATATCATACCAG	Region2.txt
R3	AGTAGACAGGTCAACCATTG	Region3.txt
R4	GAATGGAGGCAGTGCTACTA	Region4.txt
R5	TCAGCACTATATACATATGG	Region5.txt
R6	TTCTGAAGAGATAGCAACAA	Region6.txt
R7	TTGATGCTCAATGGAGGTGT	Region7.txt
```

Notes:
- `region` IDs are arbitrary labels (`R1`, `R2`, ... or `5kb`,
  `proximal`, `TSS`, anything you like). They become the column
  prefixes in the engine's outputs and the lane labels on plots.
- `sequence` is the protospacer (typically 20 bp). The PAM is optional;
  the engine matches against the promoter sequence with or without it.
  Every listed region must have a real gRNA sequence — the pipeline
  cannot place a region in space without one.
- `data_file` is the per-region proteomics table filename, resolved
  relative to the manifest's folder.

### `Region*.txt` — the per-region proteomics tables

Tab-separated tables with at minimum a protein-name column, a logFC
column, and a p-value column. An optional moderated-t column is also
recognised. Column-name matching is case-insensitive and tolerates
common aliases:

```
name	logFC	P.Value	t
POLL	1.648799	0.0164533	2.981229
ETFA	1.523198	0.0643181	2.125487
ZNF286B	0.993416	0.179416	1.463346
GPC6	0.965345	0.064466	2.124047
...
```

Recognised name aliases:
- protein column: `name`, `gene`, `protein`, `symbol`
- logFC column: `logFC`, `log2FC`, `log_fc`
- p-value column: `P.Value`, `pvalue`, `p_value`, `p`
- moderated t (optional): `t`, `t_stat`, `moderated_t`

If you don't have a `t` column, GLproxScape derives a signed z-score
from the p-value as the default weight (`weight_mode = "z"`). If you
do have moderated-t from limma, set `weight_mode = "mod_t"` in
`run_caspex()` to use it.

The TF vs. chromatin-factor classification is now controlled by the
`tf_universe` and `epi_universe` arguments to `run_caspex()` (typically
read from the bundled `TFLibrary.txt` and `EpiGenes_main.csv` —
see the next section), not by a column inside the per-region tables.

The protein column should hold HGNC symbols. Anything else (Ensembl
IDs, UniProt accessions) won't intersect cleanly with the JASPAR
motif database or the bundled TF / EpiFactors universes.

### Running the analysis

Once the folder is laid out, the runner is short. From an R session
opened in the parent folder:

```r
library(GLproxScape)

# 1) Load the inputs
inputs <- load_caspex_inputs("my_gene_analysis/inputs")

# 2) Run the full pipeline
res <- run_caspex(
  gene             = "YOUR_GENE",            # HGNC symbol
  transcript       = "canonical",          # or "ENST..." for an alt-promoter
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  upstream         = 3250,                  # bp window upstream of TSS
  downstream       = 100,                   # bp window downstream of TSS
  out_dir          = "my_gene_analysis/caspex_output",
  weight_mode      = "z",
  motif_thresh     = 0.75,                  # JASPAR PWM threshold (frac of max)
  chipatlas        = TRUE                   # set FALSE if you don't want the overlay
)

# 3) Optional: diagnostic plot pack (sigma sensitivity, jackknife, TF co-occurrence, ...)
run_caspex_extras(res, out_dir = file.path(res$out_dir, "extras"))

# 4) Optional: chromatin-factor zone-based deck (BRD4, KMT2A, SMARCA4, ...)
run_caspex_epigenetic(
  res,
  epigenetic_factors = readLines(system.file(
    "extdata/databases/EpiGenes_main.csv", package = "GLproxScape")),
  out_dir            = file.path(res$out_dir, "epigenetic")
)
```

After this, `my_gene_analysis/caspex_output/` contains the
binding-deconvolution PDF, the per-region heatmap, the gRNA-positions
plot, and the predictions CSVs.

### Picking the right transcript anchor with `check_transcripts()`

Most genes have several alternative promoters, and the canonical
Ensembl transcript isn't always the one your sgRNAs target. Running
the pipeline against the wrong transcript can silently produce zero
sgRNA matches (and therefore zero meaningful predictions). Use
`check_transcripts()` BEFORE picking the `transcript = "ENST..."`
argument:

```r
library(GLproxScape)
df <- check_transcripts(
  gene          = "YOUR_GENE",
  manifest_path = "my_gene_analysis/inputs"
)
head(df, 5)
```

The function prints a per-transcript table to the console — one row
per Ensembl transcript, with the TSS coordinate, biotype, canonical
flag, and which sgRNAs matched at which TSS-relative bp position. It
finishes with a one-line recommendation like:

```
=== Recommendation ===
Best transcript: ENST00000XXXXXXX  (N / N sgRNAs matched)
Use in run_caspex():  transcript = "ENST00000XXXXXXX"
```

Paste the recommended ENST into `run_caspex(transcript = ...)` and
you're done.

If the recommendation looks wrong, inspect the full data.frame
returned by `check_transcripts()` and pick the ENST whose TSS
coordinate places your sgRNAs where the original publication says
they sit. The canonical Ensembl transcript isn't always the one your
experiment targeted — FOXP2's HEK293 "active TSS1" (Mackenzie 2026),
for example, sits on a non-canonical alt-promoter hundreds of kb
upstream of the canonical FOXP2-201 TSS. `check_transcripts()`
surfaces this difference in its per-transcript table; the
"recommendation" line at the bottom just reports the transcript with
the most sgRNA matches, not necessarily the transcript you want to
cite. If no transcript matches any sgRNAs, double-check the HGNC
symbol, the species, and whether your sgRNA sequences are for the
genome build Ensembl is currently serving.

## `run_caspex()` parameter reference

Full list of every argument, its default, and the allowed values. Use
this as a quick lookup; the same content lives (with longer prose) in
`?run_caspex`.

#### Required inputs

- `gene` — HGNC symbol (e.g. `"FOXP2"`).
- `grnas` — named character vector of protospacer sequences (17-23 bp).
- `data_files` — named character vector of per-region proteomics table paths.

#### Annotation universes

- `tf_universe` = `NULL` — character vector of HGNC TF symbols (e.g.
  `readLines(system.file("extdata/databases/TFLibrary.txt",
  package = "GLproxScape"))`). NULL falls back to a `TFDatabase`-style
  column in the input files if present.
- `epi_universe` = `NULL` — character vector of chromatin-factor HGNC
  symbols (typically read from `EpiGenes_main.csv`).
- `tfs_only` = `TRUE` — restrict the spatial model to rows with
  `isTF = TRUE`. `FALSE` includes every protein.

#### Promoter window

- `species` = `"homo_sapiens"` — Ensembl species token.
- `transcript` = `"canonical"` — `"canonical"` | `"ENST..."` | `NA`
  (legacy gene-level union; not recommended).
- `upstream` = `2500` — bp upstream of TSS.
- `downstream` = `500` — bp downstream of TSS.

#### Spatial model + TF selection

- `pval_thresh` = `0.05` — per-region p-value gate.
- `min_regions` = `2` — minimum regions for a TF to count.
- `min_lfc` = `0` — optional logFC floor (set positive to ignore
  mildly negative values).
- `top_n` = `25` — number of TFs rendered on the spatial track plot.
- `motif_tfs` = `NULL` — explicit TF list for the motif scan; NULL
  uses the engine's `select_motif_tfs()` cut.
- `n_common`, `n_shared`, `n_specific` = `20` each — per-bucket caps
  on the motif-TF selection (top-N common across regions, top-N
  shared-focal, top-N per-region-specific).

#### Motif scan

- `motif_thresh` = `0.80` — JASPAR PWM threshold as fraction of max
  log-odds. Typical relaxation: `0.75`.
- `motif_scan_pool` = `"selected"` — `"selected"` (44-TF default cut)
  | `"spatial_all"` (also scan every spatial-model TF outside that cut;
  results live in `result$motif_results_extra`).
- `motif_score_weight` = `"none"` — `"none"` (binary threshold filter)
  | `"linear"` (amplitude × score_frac) | `"log"` (amplitude ×
  2^(score_frac - 1)).

#### Deconvolution kernel + filters

- `kernel_sigma` = `300` — Gaussian labelling-kernel width in bp.
- `min_weight_frac` = `0.15` — events below this fraction of the
  local peak amplitude are pruned.
- `min_peak_dist` = `150` — bp separation in the no-motif fallback
  peak detector.
- `merge_dist` = `100` — motif hits within this many bp are merged
  into one cluster.
- `coverage_correct` = `TRUE` — must be `TRUE`; `FALSE` (smoothed-NNLS)
  was retired in v0.1.0 and errors loudly.
- `cov_floor` = `0.05` — relative floor on the coverage denominator.
  Effective amplification cap = `1 / cov_floor` (~20×).
- `edge_guard_frac` = `0.25` — fraction-of-max-coverage floor for the
  in-support beta mask.
- `zone_peak_frac` = `0.50` — per-zone beta floor for motif retention
  (0 disables).
- `max_events_per_tf` = `30` — top-N cap per TF after merging
  (`Inf` disables).
- `merge_position` = `"argmax"` — `"argmax"` (snap to strongest motif)
  | `"centroid"` (amplitude-weighted mean).
- `max_grna_distance` = `NULL` — hard geometric cap on event-to-gRNA
  distance in bp; `NULL` resolves to `kernel_sigma` at runtime, `Inf`
  disables.
- `edge_grna_weight_cap` = `NULL` — drop events whose boundary-gRNA
  weight share exceeds this fraction in (0, 1]; `NULL` disables.

#### Region-weight mode

- `weight_mode` = `"z"` — `"z"` (signed z from p-value, default) |
  `"mod_t"` (limma moderated t) | `"lfc_pos"` | `"lfc_signed"` |
  `"lfc_x_negp"`.
- `signal_weight` = `NULL` — back-compat alias; if non-NULL overrides
  `weight_mode` for the signal track only.

#### Bootstrap diagnostic

- `position_stability` = `"none"` — `"none"` | `"wild_bootstrap"`
  (Rademacher Wild bootstrap on NNLS residuals, adds four columns to
  `binding_events`).
- `n_bootstrap` = `200L` — bootstrap draws.

#### ChIP-Atlas overlay

- `chipatlas` = `FALSE` — `TRUE` to fetch and render ChIP-seq peaks.
- `chipatlas_threshold` = `"05"` — `"05"` (Q<1e-5) | `"10"` | `"20"`.
- `chipatlas_max_experiments` = `100` — SRX cap per TF.
- `special_interest_gene` = `NULL` — character vector of TFs that
  bypass the SRX cap.
- `special_interest_cap` = `NULL` — optional integer cap for
  special-interest TFs; `NULL` = scan all SRX.
- `chipatlas_quiet` = `TRUE` — suppress per-SRX download messages.

#### Detail deck filters

- `detail_top_n` = `100` — number of TFs on the per-TF deconvolution
  detail PDF.
- `deconv_min_motif_hits` = `0` — minimum JASPAR hits in window for
  a TF to appear on the detail deck (0 disables).
- `deconv_min_max_lfc` = `0` — minimum max per-region logFC for a
  TF to appear on the detail deck (0 disables). Composes AND-style
  with `deconv_min_motif_hits`.

#### Output writing

- `out_dir` = `"caspex_output"` — output directory (created if absent).
- `save_plots` = `TRUE` — write the PDF deck.
- `plot_width` = `10`, `plot_height` = `8` — PDF dimensions in inches.

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

End-to-end walkthrough of the FOXP2 reanalysis:

```r
vignette("foxp2-mackenzie", package = "GLproxScape")
```

It covers loading the bundled inputs, picking the right transcript
anchor for an alt-promoter dataset, running the full pipeline,
inspecting the binding-events table, comparing PWM-score weighting
modes (`"none"` / `"linear"` / `"log"`), and the optional ChIP-Atlas
overlay.

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
