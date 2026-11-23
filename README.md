<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)"  srcset="man/figures/GLproxScape-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="man/figures/GLproxScape-light.svg">
    <img alt="GLproxScape" src="man/figures/GLproxScape-light.svg" width="360">
  </picture>
</p>

# GLproxScape

<!-- badges: start -->
[![R-CMD-check](https://github.com/scanozcan/GLproxScape/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scanozcan/GLproxScape/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Data: Zenodo](https://img.shields.io/badge/data-10.5281%2Fzenodo.20337799-blue)](https://doi.org/10.5281/zenodo.20337799)
[![Publication: TBD](https://img.shields.io/badge/publication-TBD-lightgrey)](#)
<!-- badges: end -->

Spatial reconstruction of chromatin occupancy landscapes from
**genomic locus proximity proteomics** (e.g. dCas9-APEX2 / CASPEX,
C-BERST, CasID).

GLproxScape treats per-region enrichments from tiled-sgRNA experiments
as indirect spatial measurements and reconstructs the latent
chromatin-binding landscape via a Gaussian labelling-kernel forward
model. Each guide's biotinylation footprint is modelled as a Gaussian
cone (default σ = 300 bp); per-region enrichment values are
forward-smeared into a continuous spatial track *s(x)*, normalised by
guide coverage *C(x)* to give the coverage-normalized track
s̃(x) = s(x) / max(C(x), c_min · max C), a kernel-regression estimate of
enrichment along the locus. Each protein is classified against curated
transcription-factor and chromatin-regulator lists and routed through
one of three analytical paths:

* **Motif-anchored event calling** for transcription factors with a
  motif entry. Motif sites inside enriched zones of the
  coverage-normalized track become binding events with bp-level
  positions and per-event weights. Motifs come from JASPAR (default),
  HOCOMOCO v12 or CIS-BP 3.10 (`motif_search_engine`).
* **Motif-free peak prediction** for transcription factors absent from
  the motif database. Peaks are called directly on the spatial signal
  *s(x)*, so motif-less TFs still receive position predictions (flagged
  `motif_based = FALSE` in the events table).
* **Zone-based detection** for chromatin regulators that lack
  sequence-specific motifs. Broad occupancy zones are detected directly
  on the labelling intensity, enabling recovery of overlapping members
  of multi-subunit complexes (e.g. MLL4/WBP7).

An optional ChIP-Atlas overlay cross-references each prediction
against independent public ChIP-seq peaks for orthogonal validation,
including tagged pulldowns (FLAG / HA / GFP / biotin-tagged TFs). With
`cell_type` set, experiments from the experimental cell type are
reported separately, the histone-mark deck uses the same cell type, and
human runs are annotated with the cell line's RNA expression from the
Human Protein Atlas.

> **New in version 0.2.0 (December 2026):** CIS-BP 3.10 motif engine
> (human, mouse and other species), tagged ChIP-Atlas pulldowns, one
> `cell_type` argument for all cell-type-aware steps, and cell-line
> expression annotation. See
> [README_update_December2026.md](README_update_December2026.md).

The package is the analysis backbone behind the GLproxScape paper
(Ozcan *et al.* 2026) and ships a bundled reanalysis of the
MacKenzie *et al.* (2026) FOXP2 dataset (3 sgRNAs) as an end-to-end
example.

<p align="center">
  <img alt="GLproxScape workflow overview" src="man/figures/glproxscape-figure1.png" width="900">
</p>

<p align="center">
  <em>(A) Experimental workflow: CASPEX
  (dCas9-APEX2) is targeted to a locus by tiled sgRNAs, biotinylated proteins
  are identified by mass spectrometry, and GLproxScape converts the per-region
  enrichments into binding events and chromatin-regulator zones. (B) Inputs,
  processing and outputs. (C) The spatial model: each sgRNA labels a Gaussian
  neighborhood, the weighted kernels sum to the signal s(x), and motif sites in
  enriched zones of the coverage-normalized track become binding events.</em>
</p>

## Install

```r
# install.packages("remotes")
remotes::install_github("scanozcan/GLproxScape", build_vignettes = TRUE)
```

`build_vignettes = TRUE` builds the FOXP2 walkthrough during install so
`vignette("foxp2-mackenzie", package = "GLproxScape")` is available;
drop it for a faster install if you don't need the vignette HTML.

**Requirements.** `run_glproxscape()` queries live web services on every
run: Ensembl (gene lookup + promoter sequence) and JASPAR (motif
PWMs), so an internet connection is required. Optional services are
fetched only when enabled and cached locally after the first call —
ChIP-Atlas peaks (`chipatlas = TRUE`), the HOCOMOCO bundle
(`motif_search_engine = "hocomoco"`), the CIS-BP species archive
(`motif_search_engine = "cisbp"`) and the Human Protein Atlas cell-line
table (`cell_type` on a human run), all cached under
`tools::R_user_dir("caspex", "cache")`.

## Quick start

```r
library(GLproxScape)

# Use the bundled FOXP2 example dataset (MacKenzie et al., 2026)
inputs_dir <- system.file("extdata/examples/foxp2_mackenzie",
                          package = "GLproxScape")
inputs <- load_glproxscape_inputs(inputs_dir)

res <- run_glproxscape(
  gene             = "FOXP2",
  transcript       = "ENST00000901759",   # MacKenzie's HEK293 active TSS1
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  upstream         = 200,
  downstream       = 2000,
  out_dir          = tempfile("foxp2_run_"),
  weight_mode      = "lfc_signed",
  pval_thresh      = 0.5,
  motif_thresh     = 0.75,
  chipatlas        = FALSE,               # set TRUE for the ChIP-Atlas overlay (slow on a cold cache)
  cell_type        = "HEK293",            # optional: cell-type-matched ChIP-Atlas, histone marks, HPA expression
  extras           = TRUE,                # diagnostic plot pack       -> <out_dir>/extras/
  epigenetic       = TRUE                 # chromatin-factor zone deck -> <out_dir>/epigenetic/
)

# Inspect the top-N predicted binding events
head(res$binding_events[order(-res$binding_events$weight), ], 10)
```

`chipatlas = FALSE` keeps this Quick start fast; flip to `TRUE` to
add the orthogonal-validation track from public ChIP-seq peaks (one
download per scanned TF on first run; cached locally thereafter).

A single `run_glproxscape()` call writes the full deliverable set to
`out_dir`: the main TF deck (gRNA layout, per-region heatmap,
motif-track PDF, per-TF binding-event detail pages, binding-events
CSV), the diagnostic plot pack under `<out_dir>/extras/`, and the
zone-based chromatin-factor deck under `<out_dir>/epigenetic/`. The
returned object exposes the corresponding R objects as
`res$spatial_df`, `res$motif_results`, `res$binding_events`,
`res$extras_result`, and `res$epigenetic_result`. Set
`extras = FALSE` or `epigenetic = FALSE` to skip either downstream
phase; pass phase-specific parameters via
`extras_args = list(...)` / `epigenetic_args = list(...)`.

The two phase functions are also exported and accept the `res` object
directly — useful for re-running just one phase while keeping the
main TF deck cached in memory:

```r
extras_out <- run_glproxscape_extras(res, out_dir = file.path(res$out_dir, "extras"))
epi_out    <- run_glproxscape_epigenetic(res, out_dir = file.path(res$out_dir, "epigenetic"))
```

## How to use it on your own data

GLproxScape expects each promoter-tiling experiment to live in a
folder with one **manifest** file plus one **proteomics table** per
gRNA region. Once the folder is laid out, the entire analysis is two
function calls — `load_glproxscape_inputs()` to read the manifest, then
`run_glproxscape()` for everything else.

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

Each table is the **differential proteomics result for one sgRNA
region** — the enrichment of that region's proximity-labelled
pulldown against a control (e.g. a no-guide / non-targeting
condition), computed by you upstream of GLproxScape (the bundled
FOXP2 example uses a per-sgRNA logFC vs. no-guide; see the paper for
the full design). GLproxScape consumes these per-region contrasts as
its spatial measurements — it does not compute differential
enrichment from raw abundances for you.

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
`run_glproxscape()` to use it.

The TF vs. chromatin-factor classification is handled inside
`run_glproxscape()` by auto-loading the bundled `TFLibrary.txt` and
`EpiGenes_main.csv`. To pin a custom universe, pass a character vector
to `tf_universe` or `epi_universe` and it overrides the bundled
default. No TF column is needed inside the per-region tables.

The protein column should hold HGNC symbols. Anything else (Ensembl
IDs, UniProt accessions) won't intersect cleanly with the JASPAR
motif database or the bundled TF / EpiFactors universes.

### Picking the right transcript anchor with `check_transcripts()`

Do this **before** running the pipeline — the transcript you pick
sets the TSS that every coordinate is measured from. Most genes have
several alternative promoters, and the canonical Ensembl transcript
isn't always the one your sgRNAs target. Running against the wrong
transcript can silently produce zero sgRNA matches (and therefore
zero meaningful predictions). Use `check_transcripts()` to choose
the `transcript = "ENST..."` argument:

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
Use in run_glproxscape():  transcript = "ENST00000XXXXXXX"
```

Paste the recommended ENST into `run_glproxscape(transcript = ...)` and
you're done.

The recommendation just picks the transcript with the most sgRNA
matches — usually right, but not always. The canonical Ensembl
transcript isn't necessarily the one your experiment targeted: a
gene can have several promoters, and the one active in your cell
type may sit on a non-canonical alt-promoter far from the canonical
TSS. (FOXP2's HEK293 "active TSS1" (MacKenzie *et al.*, 2026), for
example, sits hundreds of kb upstream of the canonical FOXP2-201
TSS.) You know which promoter your sgRNAs were designed against, so
when the recommendation looks off, inspect the full data.frame and
pick the ENST whose TSS coordinate places your sgRNAs where you
designed them to land. If no transcript matches any sgRNAs,
double-check the HGNC symbol, the species, and whether your sgRNA
sequences come from the genome build Ensembl is currently serving.

### Running the analysis

Once the folder is laid out, the runner is short. From an R session
opened in the parent folder:

```r
library(GLproxScape)

# 1) Load the inputs
inputs <- load_glproxscape_inputs("my_gene_analysis/inputs")

# 2) Run the full pipeline
res <- run_glproxscape(
  gene             = "YOUR_GENE",            # HGNC symbol
  transcript       = "canonical",          # or "ENST..." for an alt-promoter
  grnas            = inputs$grnas,
  data_files       = inputs$data_files,
  upstream         = 3250,                  # bp window upstream of TSS
  downstream       = 100,                   # bp window downstream of TSS
  out_dir          = "my_gene_analysis/glproxscape_output",
  weight_mode      = "z",
  motif_thresh     = 0.75,                  # JASPAR PWM threshold (frac of max)
  chipatlas        = TRUE,                  # set FALSE if you don't want the overlay
  extras           = TRUE,                  # diagnostic plot pack       -> <out_dir>/extras/
  epigenetic       = TRUE                   # chromatin-factor zone deck -> <out_dir>/epigenetic/
)
```

**Setting the promoter window.** `upstream` and `downstream` define
the analysis window as `[-upstream, +downstream]` bp around the TSS
of the chosen transcript. The window must span the genomic stretch
your sgRNAs tile across — any sgRNA whose match position falls
outside it is silently dropped, the same failure mode as picking the
wrong transcript. The `check_transcripts()` table above reports each
sgRNA's TSS-relative bp position for your chosen ENST; size the
window to comfortably bracket that range (the example's asymmetric
`3250 / 100` simply reflects sgRNAs tiling upstream of the TSS).
There is no single right value — match it to where your guides sit,
not to the defaults.

`extras` and `epigenetic` are TRUE by default; set either to `FALSE`
to skip the corresponding phase. Pass phase-specific parameters via
`extras_args = list(...)` / `epigenetic_args = list(...)` without
bloating the main signature.

After this, `my_gene_analysis/glproxscape_output/` contains the
binding-event deck (`binding_events.pdf`), the per-region heatmap, the gRNA-positions
plot, and the predictions CSVs.

**Choosing the region weights.** `weight_mode` sets how each region's
per-protein weight is computed from its `logFC` and p-value; only the
positive part drives spatial placement. The default `"z"` is the signed z
from the p-value (`sign(logFC) × qnorm(p / 2, lower.tail = FALSE)`), which
keeps the weight scale consistent across datasets with or without a `t`
column. `"mod_t"` uses the limma moderated t when present (wider scale, same
event ranking), `"lfc_pos"` and `"lfc_signed"` use the effect size only, and
`"lfc_x_negp"` is the legacy `max(logFC, 0) × −log10(p)`.

**Setting the cell type.** `cell_type` (e.g. `"HEK293"`) is the experimental
cell line or cell type, passed to every cell-type-aware step: ChIP-Atlas
experiments from this cell type are scanned first and reported separately
(`chipatlas_support_by_cell_type.csv`, `res$chipatlas_support`); the
histone-mark deck of `run_glproxscape_epigenetic()` uses it as its cell-type
filter; and for `species = "homo_sapiens"` the detected TFs are annotated
with the cell line's RNA expression from the Human Protein Atlas
(`tf_expression_annotation.csv`, `tf_expression_families.pdf`,
`res$expression`). For other species, use `annotate_tf_expression()` with
`load_cell_line_expression(file = ...)`. Without `cell_type`, these steps
are skipped and ChIP-Atlas support is not split by cell type.

### Reading the output

The primary result is the binding-events table — written to CSV and
returned as `res$binding_events`. One row per predicted binding event:

```r
head(res$binding_events[order(-res$binding_events$weight), ], 10)
```

Key columns:
- `tf` — the predicted factor (HGNC symbol).
- `position` — **TSS-relative bp** of the event: negative is
  upstream of the TSS, positive is downstream. This is the same
  coordinate frame as the `check_transcripts()` sgRNA positions, so
  events land in the span your guides tile.
- `weight` — the event weight (coverage-normalized enrichment at the
  event, optionally scaled by the motif score).
  Higher = stronger predicted binding; rank by this to read off the
  top calls.
- `motif_based` — `TRUE` if the event is anchored to a JASPAR/HOCOMOCO
  motif hit (sequence-specific TF path), `FALSE` if it came from the
  no-motif peak detector.
- `distance_to_nearest_grna` / `local_coverage` — how far the event
  sits from the nearest sgRNA and the labelling coverage there;
  events far from any guide or in low-coverage gaps are less reliable.

The zone-based predictions for epigenetic factors live separately on
`res$epigenetic_result` and under `<out_dir>/epigenetic/`. The binding-event deck visualises the
same events as bubbles along the promoter; the per-region heatmap and
gRNA-positions plot are diagnostics for the underlying signal.

## Full example with all options

All arguments of `run_glproxscape()` with their defaults, using the example
from "Running the analysis" and grouped by step. Values marked "default" can
be left out; the others turn on optional features (CIS-BP motifs,
cell-type-aware ChIP-Atlas, histone marks and expression, the positional
bootstrap). Full descriptions: `?run_glproxscape`.

```r
res <- run_glproxscape(
  # -- Inputs ---------------------------------------------------------------
  gene                  = "YOUR_GENE",         # HGNC symbol
  grnas                 = inputs$grnas,
  data_files            = inputs$data_files,
  tf_universe           = NULL,                # default: bundled TF list (Lambert et al. 2018)
  epi_universe          = NULL,                # default: bundled EpiFactors list
  species               = "homo_sapiens",      # default; e.g. "mus_musculus"
  transcript            = "canonical",         # or "ENST..." for an alt-promoter
  cell_type             = "HEK293",            # experimental cell line: ChIP-Atlas by cell type,
                                               #   histone-mark filter, HPA expression (human)
  out_dir               = "my_gene_analysis/glproxscape_output",

  # -- Promoter window ------------------------------------------------------
  upstream              = 3250,                # bp upstream of the TSS
  downstream            = 100,                 # bp downstream of the TSS

  # -- Protein filters and TF selection ---------------------------------------
  pval_thresh           = 0.05,                # default
  min_regions           = 2,                   # default
  min_lfc               = 0,                   # default
  tfs_only              = TRUE,                # default
  top_n                 = 25,                  # default
  motif_tfs             = NULL,                # default: chosen from n_common / n_shared / n_specific
  n_common              = 20,                  # default
  n_shared              = 20,                  # default
  n_specific            = 20,                  # default

  # -- Region weights ----------------------------------------------------------
  weight_mode           = "z",                 # default; "mod_t", "lfc_pos", "lfc_signed", "lfc_x_negp"
  signal_weight         = NULL,                # default: same as weight_mode

  # -- Motif scan ----------------------------------------------------------------
  motif_search_engine   = "cisbp",             # default "jaspar"; or "hocomoco"
  motif_thresh          = 0.75,                # default 0.80 (fraction of max PWM score)
  motif_score_weight    = "none",              # default; "linear", "log"
  motif_scan_pool       = "selected",          # default; "spatial_all"
  hocomoco_version      = "v12",               # default (used with "hocomoco")
  hocomoco_species      = NULL,                # default: from species
  cisbp_status          = "direct",            # default; "inferred" for species with few motifs
  cisbp_prefer          = "in_vitro",          # default; "none"
  cisbp_pick            = "medoid",            # default; "ic"
  cisbp_species         = NULL,                # default: from species

  # -- Labeling kernel and event calling -----------------------------------------
  kernel_sigma          = 300,                 # default (bp)
  cov_floor             = 0.05,                # default
  edge_guard_frac       = 0.25,                # default
  min_weight_frac       = 0.15,                # default
  zone_peak_frac        = 0.5,                 # default
  min_peak_dist         = 150,                 # default (motif-free peaks)
  merge_dist            = 100,                 # default
  merge_position        = "argmax",            # default; "centroid"
  max_events_per_tf     = 30,                  # default
  max_grna_distance     = NULL,                # default: kernel_sigma
  edge_grna_weight_cap  = NULL,                # default: no cap

  # -- Positional stability (optional) ----------------------------------------------
  position_stability    = "wild_bootstrap",    # default "none"
  n_bootstrap           = 200L,                # default

  # -- ChIP-Atlas overlay ----------------------------------------------------------------
  chipatlas             = TRUE,                # default FALSE
  chipatlas_threshold   = "05",                # default; "10", "20"
  chipatlas_max_experiments = 100,             # default (newest per TF; Inf = all)
  chipatlas_tagged      = TRUE,                # default; FALSE = TF-antibody experiments only
  chipatlas_genome      = NULL,                # default: from species
  special_interest_gene = c("CTCF"),           # default NULL; TFs scanned without the cap
  special_interest_cap  = NULL,                # default: no cap for those TFs
  chipatlas_quiet       = TRUE,                # default

  # -- Detail deck --------------------------------------------------------------------------
  detail_top_n          = 100,                 # default
  deconv_min_motif_hits = 0,                   # default
  deconv_min_max_lfc    = 0,                   # default

  # -- Output ---------------------------------------------------------------------------------
  save_plots            = TRUE,                # default
  plot_width            = 10,                  # default
  plot_height           = 8,                   # default

  # -- Downstream phases ------------------------------------------------------------------------
  extras                = TRUE,                # default; diagnostic plot pack -> <out_dir>/extras/
  epigenetic            = TRUE,                # default; chromatin-factor zones -> <out_dir>/epigenetic/
  extras_args           = list(),              # default
  epigenetic_args       = list()               # default; e.g. list(histone_cell_type = "HEK293T")
)
```

## Bundled example data

The package ships one end-to-end example dataset under
`inst/extdata/examples/foxp2_mackenzie/` — a reanalysis of
MacKenzie *et al.* (2026) at the FOXP2 promoter using 3 sgRNAs.
Resolve the path at runtime with
`system.file("extdata/examples/foxp2_mackenzie", package = "GLproxScape")`.
The folder is a self-contained input bundle: a `grnas.tsv` manifest
(region → protospacer + per-region file) plus per-region `Region*.txt`
proteomics tables with `logFC` and `P.Value` columns.

The TF and chromatin-factor universes used by `run_glproxscape()` are
bundled separately under `inst/extdata/databases/`
(`TFLibrary.txt`, `EpiGenes_main.csv`, `EpiGenes_complexes.csv`),
sourced from the Lambert *et al.* TF list and the EpiFactors database
respectively. These are loaded automatically; override by passing
`tf_universe` / `epi_universe` explicitly.

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

The pipeline runs in one call (`run_glproxscape`) that internally chains:

1. **Gene lookup** — resolves the canonical (or explicitly pinned)
   Ensembl transcript and fetches the promoter sequence over
   `[-upstream, +downstream]` bp of the TSS.
2. **gRNA matching** — exact-match each protospacer (and its reverse
   complement) against the promoter sequence to anchor the labelling
   kernel.
3. **Spatial model** — per-protein per-region significance gating
   (`pval_thresh`, `min_regions`), then `compute_spatial` aggregates
   into a TF-level summary scoring composite / specificity.
4. **Motif scan** — for each TF in the deck roster, fetches the
   position weight matrix from JASPAR (default) or HOCOMOCO v12
   (`motif_search_engine = "hocomoco"`) and scans the promoter at
   `motif_thresh × max_score`. The same fractional threshold is used
   for both engines; HOCOMOCO matrices are typically higher information
   content so the absolute stringency is slightly stricter.
5. **Coverage-normalized event calling** — `s̃(x) = s(x) / max(C(x), cov_floor ·
   max(C))` thresholded into zones, with one event emitted per motif hit
   inside each zone whose local value reaches the zone threshold. Optional PWM-score weighting reshapes per-event
   amplitudes (`motif_score_weight = "none" | "linear" | "log"`).
6. **Edge guards + merge** — events outside the gRNA support cone
   (`edge_guard_frac`, `max_grna_distance`) are dropped; closely-spaced
   events within `merge_dist` are merged (weights summed; position = the
   strongest motif by default).
7. **(Optional) ChIP-Atlas overlay** — public ChIP-seq peaks for each
   deck TF rendered beneath the predicted event bubbles for independent
   validation. Per-SRX BEDs cached under `tools::R_user_dir("caspex",
   "cache")`.
8. **Plotting + CSV write-out** — multi-page binding-event PDF, motif
   track, gRNA-positions plot, per-TF binding-events CSV, spatial
   predictions CSV.
9. **Diagnostic extras** — `run_glproxscape_extras()` produces the
   per-TF one-pagers, σ-sensitivity sweep, gRNA jackknife, TF-pair
   co-occurrence triangle, and other diagnostic plots in
   `<out_dir>/extras/`. Skipped when `extras = FALSE`.
10. **Chromatin-factor zone deck** — `run_glproxscape_epigenetic()` runs
    the zone-based path for epigenetic factors that lack a
    sequence-specific motif, plus the histone-marks
    locus map and (when both EpiGenes CSVs are bundled or supplied) a
    per-complex deck showing co-localising members. Output in
    `<out_dir>/epigenetic/`. Skipped when `epigenetic = FALSE`.

The Methods section of the accompanying paper gives the precise
mathematical framing for every step.

## Reproducibility

The MacKenzie FOXP2 analysis reported in the paper is reproducible
end-to-end from the bundled `inst/extdata/examples/foxp2_mackenzie/`
folder; the corresponding `run_glproxscape()` parameter settings live in
the FOXP2 vignette. The other published paper analyses (Myers *et al.*
TERT and MYC, Pizzolato *et al.* FOXQ1, Gao *et al.* mouse Ripk3)
ship as a Zenodo companion deposit alongside the manuscript
([doi:10.5281/zenodo.20337799](https://doi.org/10.5281/zenodo.20337799));
only the FOXP2 example is bundled inside the package itself.

The `transcript = "canonical"` default in `lookup_gene()` keeps
TSS-relative coordinates stable across Ensembl releases. For
datasets that target alternative promoters — e.g. MacKenzie FOXP2's
HEK293 "active TSS1", which sits on a non-canonical alt-promoter
rather than the canonical FOXP2-201 TSS — pin the right ENST ID
explicitly (`transcript = "ENST00000901759"`); the diagnostic helper
`check_transcripts()` surfaces the per-transcript sgRNA-match table
when you're picking it.

### Cross-species runs

Non-human datasets work end-to-end by passing the appropriate
`species` argument (e.g. `species = "mus_musculus"` for the Gao
*et al.* mouse Ripk3 reanalysis). `run_glproxscape()` then auto-routes
the Ensembl gene lookup, promoter sequence fetch, and ChIP-Atlas
overlay to the matching genome. Ensembl currently serves GRCm39
for mouse while ChIP-Atlas's mouse corpus remains in mm10 (GRCm38);
the engine detects this assembly-frame mismatch and lifts gene
coordinates via the Ensembl archive at `e102.rest.ensembl.org`
before filtering ChIP-Atlas peaks, so the overlay aligns
correctly. Supported species → assembly mappings include
`homo_sapiens` → `hg38`, `mus_musculus` → `mm10`,
`rattus_norvegicus` → `rn7`, plus fly, worm, zebrafish, chicken,
yeast, and a few others; override with `chipatlas_genome = "..."`
when a specific older build is needed.

## Citation

If you use GLproxScape in your work, please cite:

> Ozcan SC, Sergi B, Yildirim B, Cagiral U, Gonen M, Acilan C.
> **GLproxScape reconstructs spatial chromatin occupancy landscapes
> from tiled genomic locus proteomics.** 2026.

This citation will be updated with the full publication details once the
paper is published.

## License

[MIT](LICENSE). See [LICENSE.md](LICENSE.md) for the full text.

## Issues and contributions

Bug reports, feature requests, and pull requests are welcome at
[github.com/scanozcan/GLproxScape](https://github.com/scanozcan/GLproxScape).
