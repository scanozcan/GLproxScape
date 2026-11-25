# GLproxScape update — December 2026 (version 0.2.0)

This update adds three features to GLproxScape and one argument that ties
the cell-type-aware steps together. Existing analyses are unchanged unless
the new options are used: the default motif source is still JASPAR, and all
cell-type steps are skipped when `cell_type` is not given.

## 1. CIS-BP motif engine

`run_glproxscape(..., motif_search_engine = "cisbp")` uses CIS-BP 3.10
(Weirauch et al. 2014; Lambert et al. 2019), which includes the CodeBook
motifs for previously uncharacterized human TFs (Jolma et al. 2026), as a
third motif source next to JASPAR and HOCOMOCO v12.

- **Download and cache.** The CIS-BP species archive (TF information and
  PWMs) is downloaded once and cached under
  `tools::R_user_dir("glproxscape", "cache")`. If the automatic download is
  blocked, download the archive from
  <https://cisbp.ccbr.utoronto.ca/bulk.php> and import it with
  `download_cisbp_bundle(species, from_file = "archive.zip")`.
- **One motif per TF.** CIS-BP lists several motifs per TF, including
  motifs imported from JASPAR and HOCOMOCO. With `cisbp_prefer = "in_vitro"`
  (default), the motif is chosen from the best available tier:
  1. high-throughput in vitro motifs (HT-SELEX, GHT-SELEX, PBM, SMiLE-seq),
  2. other direct (experimental) motifs,
  3. motifs imported from JASPAR or HOCOMOCO.

  Within the tier, the motif most similar to the TF's other motifs is used
  (`cisbp_pick = "medoid"`; `"ic"` picks the highest information content).
- **Same threshold meaning.** CIS-BP probability matrices are converted to
  log-odds exactly as HOCOMOCO matrices, so `motif_thresh` means the same
  for both.
- **Species.** The CIS-BP species follows `species` (e.g.
  `"mus_musculus"` → `"Mus_musculus"`), so CIS-BP also covers species that
  HOCOMOCO does not. For species with few experimental motifs,
  `cisbp_status = "inferred"` uses motifs inferred from DNA-binding-domain
  similarity when no direct motif exists.
- **New functions.** `fetch_cisbp_pwm()`, `list_cisbp_tfs()`,
  `download_cisbp_bundle()`, `clear_cisbp_memory_cache()`.

```r
res <- run_glproxscape(gene = "FOXP2", grnas = grnas, data_files = data_files,
                  motif_search_engine = "cisbp")
fetch_cisbp_pwm("MAX")$source   # which CIS-BP source the motif came from
```

## 2. Tagged ChIP-Atlas pulldowns

ChIP-Atlas files many ChIP-seq experiments of tagged TFs (FLAG, HA, V5, GFP,
biotin) under "Epitope tags", "GFP" or "Biotin" instead of the TF. These
experiments are now included in the ChIP-Atlas overlay by default, together
with the TF-antibody experiments (about 875 extra human experiments covering
384 TFs). Each peak table has a `pulldown` column (`"native"` or `"tagged"`).
Set `chipatlas_tagged = FALSE` to use TF-antibody experiments only. The
pulldown tables for hg38 and mm10 ship with the package.

## 3. Cell-line expression (human)

For human runs with `cell_type`, the TFs with predicted events are annotated
with the cell line's RNA expression from the Human Protein Atlas (nTPM,
about 1,200 cell lines; the table is downloaded once and cached).

- TFs with motif-anchored events are grouped by DNA-binding-domain family
  (Lambert et al. 2018; table shipped with the package). C2H2 zinc fingers
  and TFs without an assigned family are not grouped, because their members
  do not share motifs.
- Within a family, family members whose events overlap (Jaccard index of
  events ≥ 0.3, events within 25 bp counted as the same site) share sites,
  and the most highly expressed of these is flagged as the preferred
  candidate.
- The annotation does not filter or re-weight predictions.
- Outputs: `tf_expression_annotation.csv`, `tf_family_shared_sites.csv`,
  `tf_expression_families.pdf` and `res$expression`.
- Other species or own data: `annotate_tf_expression(res,
  load_cell_line_expression(file = "my_tpm.tsv"))`.
- **New functions.** `load_cell_line_expression()`, `load_tf_families()`,
  `annotate_tf_expression()`, `plot_expression_families()`.

## One `cell_type` argument

`run_glproxscape(..., cell_type = "HEK293")` passes the experimental cell type to
every step that uses it:

| Step | What `cell_type` does |
|---|---|
| ChIP-Atlas overlay | experiments from this cell type (HEK293 sublines such as 293, 293T and 293FT count as one cell line) are scanned first (kept within the per-TF cap) and flagged (`cell_type_matched`); event support from matched vs other cell types is written to `chipatlas_support_by_cell_type.csv` (`chipatlas_support_by_cell_type()`) |
| Histone-mark deck (`run_glproxscape_epigenetic()`) | used as the cell-type filter (replaces the former `"HEK293T"` default of `histone_cell_type`, which still overrides it when given) |
| Expression annotation | Human Protein Atlas cell line (human runs only) |

## Renamed functions

The API now uses GLproxScape names. The old names still work and call the new
functions, with a deprecation message.

| Old name | New name |
|---|---|
| `run_caspex()` | `run_glproxscape()` |
| `run_caspex_epigenetic()` | `run_glproxscape_epigenetic()` |
| `run_caspex_extras()` | `run_glproxscape_extras()` |
| `load_caspex_inputs()` | `load_glproxscape_inputs()` |
| `plot_binding_deconvolution()` | `plot_binding_events()` |

Also renamed: the per-TF binding-event deck is now written as
`binding_events.pdf` (was `binding_deconvolution.pdf`); plot titles read
"<TF> - binding events"; the y axis of the binding-event plots reads
"GLproxScape signal (...)"; the default output folder is
`glproxscape_output` (was `caspex_output`); the source files are `R/glproxscape_*.R`. The cache
folder is now `tools::R_user_dir("glproxscape", "cache")` (was `"caspex"`); an
existing `caspex` cache is moved there on first use and a link is left at the old
path, so downloaded ChIP-Atlas, HOCOMOCO, CIS-BP and Human Protein Atlas files are
reused.

## Fixes

- ChIP-Atlas: the cached `experimentList.tab` was rejected and downloaded
  again because recent releases start with non-SRA records; the file check
  now looks beyond the first line. The download may take up to one hour on
  slow connections instead of timing out after 60 s.
- Ensembl requests now wait up to 60 s and are retried up to three times.
- Cell-type matching (ChIP-Atlas TF scan and histone-mark deck) ignores
  case, spaces and punctuation ("MIA PaCa-2" = "MiaPaca-2" = "miapaca2") and
  accepts subline suffixes after a separator ("HeLa" matches "HeLa-S3"), but
  no longer matches names that only share a prefix ("PANC1" no longer
  matches "PANC 10.05", "HeLa" no longer matches "HEL"). HEK293 sublines
  (293, 293T, 293FT) count as one cell line.
- The epigenetic deck now uses the same ChIP-Atlas genome as the TF scan
  (previously always hg38).
- Mouse ChIP-Atlas (mm10): the Ensembl archive REST server used for the
  GRCm39 -> GRCm38 lift-over (e102.rest.ensembl.org) has been retired, which
  left mouse runs without ChIP-Atlas peaks. The TSS is now mapped with the
  Ensembl /map endpoint of the current REST server.

## References

- Weirauch MT et al. Determination and inference of eukaryotic transcription
  factor sequence specificity. *Cell* 2014;158:1431–1443.
- Lambert SA et al. Similarity regression predicts evolution of transcription
  factor sequence specificity. *Nat Genet* 2019;51:981–989.
- Jolma A et al. An expanded codebook of human transcription factor
  DNA-binding specificity. *Nature* 2026;657:275–283.
- Lambert SA et al. The human transcription factors. *Cell* 2018;172:650–665.
- Uhlén M et al. Tissue-based map of the human proteome. *Science*
  2015;347:1260419.
- Oki S et al. ChIP-Atlas: a data-mining suite powered by full integration
  of public ChIP-seq data. *EMBO Rep* 2018;19:e46255.
