# FOXP2 promoter — Mackenzie 2026 reanalysis inputs

This folder contains the per-region proteomics tables and sgRNA
manifest for the GLproxScape reanalysis of the Mackenzie *et al.* 2026
FOXP2 dataset. It is the canonical example dataset shipped with the
package and is used end-to-end in the FOXP2 vignette.

## Files

- `grnas.tsv` — sgRNA manifest (region ID → protospacer → per-region
  data file). Three sgRNAs targeting the FOXP2 promoter; genomic
  coordinates listed in the manifest header.
- `Region1.txt` (g1), `Region2.txt` (g2), `Region3.txt` (g3) — per-sgRNA
  protein-level enrichment tables with columns `name`, `logFC`,
  `P.Value`, `t`.

## Original publication

> Mackenzie *et al.* (2026). *Mapping transcription factor binding
> at promoters by dCas9-APEX2 proximity labelling.*
> Molecular & Cellular Proteomics (in press).

The paper is open access; the supplementary file `mmc3` contains the
authors' Proteome Discoverer 2.5 grouped-abundance output for an
18-plex TMTpro experiment (3 sgRNAs × 3 replicates plus NoG+ and NoG-
controls). The Region*.txt tables in this folder are NOT verbatim
copies of any paper supplementary; they are this project's re-processed
per-sgRNA contrasts derived from that grouped-abundance matrix using a
limma pipeline (script reference below).

## How the Region*.txt tables were produced

The processing script is `mackenzie_reanalysis/input_prep/0-build-mackenzie-inputs.R`
in the GLproxScape paper's companion Zenodo deposit. In one paragraph:

The 18 TMT channels are demultiplexed using the channel-to-sample
mapping from the paper's supplementary methods (3 reps × {g1, g2, g3,
NoG+, NoG-} plus pool and empty channels). Channel intensities are
pool-anchored: each channel is scaled so its column-sum matches the
pool channel (TMT-126 / 133N), correcting for residual loading
differences not absorbed by Proteome Discoverer's internal
normalization. The matrix is collapsed from accessions to gene symbols
by summation, log2-transformed with a floor on zeros, optionally
filtered for streptavidin/bead contaminants (proteins enriched in
NoG+ vs NoG- at adj.P < 0.10 and logFC > 0.5), and fed to limma
`lmFit + eBayes(contrasts.fit(...))` with three per-sgRNA contrasts:
g1 − NoG+, g2 − NoG+, g3 − NoG+. The output of each contrast is
written to the corresponding `RegionN.txt` (logFC, moderated-t,
P.Value, gene symbol). All three Region files share the same protein
universe, so `load_caspex_inputs()` can stack them into one long
data.frame.

## Reproducibility note

This is the same dataset used in Figure 2 of the GLproxScape preprint
(Ozcan *et al.*, in preparation). To reproduce the figure end-to-end
from the bundled inputs:

```r
library(GLproxScape)
inputs_dir <- system.file("extdata/examples/foxp2_mackenzie",
                          package = "GLproxScape")
inputs <- load_caspex_inputs(inputs_dir)
res <- run_caspex(
  gene       = "FOXP2",
  transcript = "ENST00000901759",   # Mackenzie's HEK293 "active TSS1"
  grnas      = inputs$grnas,
  data_files = inputs$data_files,
  upstream   = 200,
  downstream = 2000,
  weight_mode  = "lfc_signed",
  pval_thresh  = 0.5,
  motif_thresh = 0.75,
  out_dir      = tempfile("foxp2_run_")
)
```

`run_caspex()` auto-loads the bundled TF and epigenetic-factor universes
(`inst/extdata/databases/TFLibrary.txt`,
`inst/extdata/databases/EpiGenes_main.csv`) when `tf_universe` and
`epi_universe` are left at their `NULL` defaults — no manual list-loading
needed. To pin a custom universe, pass a character vector to either
argument and it overrides the bundled default.

The paper's Zenodo deposit will additionally contain (i) the
`mmc3_proteome_discoverer_output.xlsx` raw input to the pre-processing
script, (ii) the pre-processing script itself, and (iii) the runner
scripts that drive the full paper-figure regeneration via the
package's exported functions.

## License

Re-distributed under the package's MIT license. Original raw data
remains under the journal's open-access terms; please cite Mackenzie
*et al.* 2026 when using these tables.
