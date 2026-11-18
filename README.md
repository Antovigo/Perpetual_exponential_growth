# Perpetual exponential growth

## List of variants from the experiments

The `variants/` folder has every variant detected in the sequenced populations:

| File | Contents |
|---|---|
| `variants/variant_frequencies.tsv.gz` | Frequency of each variant in each population at each sequenced time point (one row per variant, population and time point) |
| `variants/variant_catalogue.tsv` | One row per variant and population: location, effect, maximum frequency, and when it appeared and fixed |

Main columns:

- `pop`: population sampled; `origin`: population the variant came from (differs from `pop` for variants carried over between populations)
- `gen`: time point, in units of 2,500 generations (0 to 24, i.e. 0 to 60,000 generations)
- `freq` (frequencies table only): variant frequency, from 0 to 1
- `chrom`: `U00096` (*E. coli* MG1655 chromosome) or `JX871397` (prophage Φ80); `pos`, `ref`, `alt`: position and alleles
- `type`: `SNP`, `INS`, `DEL`, `SUB`, `REP` (repeat-length change) or `REV` (reversion); `cat`: effect class (missense, synonymous, nonsense, frameshift, intergenic…)
- `context`: affected gene, or `upstream/downstream` genes for intergenic variants; `effect`: amino-acid change or equivalent; `product`: gene product
- `max_f`: maximum frequency reached; `fixation`: time point at which the variant fixed (`Inf` if it never fixed); time point at which it appeared: `apparition` in the frequencies table, `gen` in the catalogue
- `vid`: variant identifier (gene and number); `tid`: identifier of the variant's trajectory in one population

Both tables are produced by `analysis/genomics/04_gene_attribution/C2. variant_catalogue.ipynb`
(as `d_labeled.tsv` and `v.tsv`), from the breseq variant calls in `datasets/genomics/vcf/`.

## Figure code and data

Notebooks and datasets needed to regenerate the figures of the manuscript.

### Layout

```
R/                 preamble.R (plotting helpers, psave) and theme.R, sourced by every R notebook
analysis/<topic>/  analysis notebooks; plots they make go to Plots/ (with GG/*.rds ggplot objects)
datasets/<topic>/  input data (committed) and intermediate tables (written by the notebooks, gitignored)
variants/          ready-to-use variant tables (see above)
figures/           make_sticking_figures.ipynb, make_collage.ipynb, static/ images, output/ final figures
```

`figures/output/` already holds reference copies of the figures (PNG) as the code produces them, as do
`analysis/genomics/07_mutation_statistics/Plots/prom_mut_plot.png`,
`analysis/genomics/08_convergence/Plots/multiplicity_vs_null.png` and
`analysis/growth_curves/2025-03-23_AF_gamma/Plots/growth_curve_reproductibility.png`;
running the notebooks overwrites them.

All paths are relative to each notebook's folder, so run a notebook from its own directory
(the Jupyter default).

### Regenerating the figures

Run the notebooks in this order; each one only needs the outputs of those above it.
All of them were run this way, from the committed data, to produce the notebook outputs and figures in the repository.

**Genomics**

1. `analysis/genomics/01_variant_calling/A1. vcf_to_mutations.ipynb`
1. `analysis/genomics/01_variant_calling/A2. coarse_filtering.ipynb`
1. `analysis/genomics/04_gene_attribution/C0. parse_gff3.ipynb`
1. `analysis/genomics/02_structural_variation/B2. unassigned_plots.ipynb`
1. `analysis/genomics/03_variant_filtering/B1. fix_continuity.ipynb`
1. `analysis/genomics/03_variant_filtering/B2. disentangle_variants.ipynb`
1. `analysis/genomics/03_variant_filtering/B3. background_sequence_coverage.ipynb`
1. `analysis/genomics/03_variant_filtering/B4. mutation_filtering.ipynb`
1. `analysis/genomics/04_gene_attribution/C1. gene_attribution.ipynb`
1. `analysis/genomics/04_gene_attribution/C2. variant_catalogue.ipynb`
1. `analysis/genomics/05_genome_mapping/E1. affected_genes_map.ipynb`
1. `analysis/genomics/06_mutators/F0. variant_to_mutator_bg.ipynb`
1. `analysis/genomics/06_mutators/F1. mutators.ipynb`
1. `analysis/genomics/06_mutators/F3. counterfactual_effects.ipynb`
1. `analysis/genomics/07_mutation_statistics/G1. mutation_statistics.ipynb`
1. `analysis/genomics/08_convergence/N0. null_model.ipynb`
1. `analysis/genomics/09_phage_fhuA/P1. fhuA.ipynb`
1. `analysis/genomics/09_phage_fhuA/P2. phage_dynamics.ipynb`
1. `analysis/genomics/10_mutation_traces/experiment_overview.ipynb`
1. `analysis/genomics/10_mutation_traces/mutation_traces.ipynb`

**Chromosome replication**

1. `analysis/chromosome_replication/coverage_analysis/O2. process_coverage_data.ipynb`
1. `analysis/chromosome_replication/growth_curves/curve_fitting_MF_epsilon.ipynb`
1. `analysis/chromosome_replication/coverage_analysis/O3. plots_MF_epsilon.ipynb`

**Growth curves**

1. `analysis/growth_curves/2024-06-11_beta/Esch_beta_fitting.ipynb`
1. `analysis/growth_curves/2024-06-11_beta/Esch_beta_plots_pure_only.ipynb`
1. `analysis/growth_curves/2025-02-27_AF_alpha/curve_fitting_esch_AF_alpha.ipynb`
1. `analysis/growth_curves/2025-03-23_AF_gamma/curve_fitting_esch_AF_gamma.ipynb`
1. `analysis/growth_curves/2025-03-23_AF_gamma/growth_curve_plots_esch_AF_gamma.ipynb`
1. `analysis/growth_curves/2025-04-05_BC_alpha/BC_alpha_curve_fitting.ipynb`
1. `analysis/growth_curves/2025-04-05_BC_alpha/BC_alpha_curve_plots.ipynb`
1. `analysis/growth_curves/2025-04-15_DE_alpha/Esch_DE_fitting.ipynb`
1. `analysis/growth_curves/2025-04-15_DE_alpha/Esch_DE_plots.ipynb`

**Morphology**

1. `analysis/morphology/plots_phase.ipynb`

**Adhesion and motility**

1. `analysis/adhesion/fim_switch/L1. fim_switching_plots.ipynb`
1. `analysis/adhesion/motility_summary/motility_summary.ipynb`
1. `analysis/adhesion/stickiness/cv_plots_2.ipynb`

**Phage**

1. `analysis/phage/C20-21_resistance/plots.ipynb`
1. `analysis/phage/phage_collapse/phage_collapse_plot.ipynb`
1. `analysis/phage/phage_production/plots_esch_AF_phage_production.ipynb`

**Stationary survival**

1. `analysis/stationary_survival/ss_beta_plots.ipynb`

**Figures**

1. `figures/make_sticking_figures.ipynb`
1. `figures/make_collage.ipynb`

`figures/make_collage.ipynb` assembles the panels into the final figures in `figures/output/`:
machine overview (panels D–G), genomics, latent traits, chromosome replication, phage,
and the supplementary figures (mutators, structural variants, phage, separate growth dynamics).
The supplementary genomics figure (mutation multiplicity vs null) is
`analysis/genomics/08_convergence/Plots/multiplicity_vs_null.png`, made by `N0. null_model.ipynb` (Python). The growth-curve reproducibility figure is written by
`analysis/growth_curves/2025-03-23_AF_gamma/growth_curve_plots_esch_AF_gamma.ipynb`,
and the promoter-mutation figure by `analysis/genomics/07_mutation_statistics/G1. mutation_statistics.ipynb`.

### Upstream steps (included for provenance, not run here)

These steps need breseq outputs or raw sequencing reads, which are too large for this repository.
Their results are committed as the starting points listed under Datasets.

| Step | Code | Produces |
|---|---|---|
| Align reads and call variants with breseq (per sample, on a SLURM cluster) | `analysis/genomics/00_breseq/make_folders.sh`, `esch_run_breseq.sh` | breseq output folders |
| Collect VCF, GD, HTML and coverage files from the breseq outputs | `analysis/genomics/00_breseq/gather_breseq_output.py` | `datasets/genomics/vcf/` (VCF files) |
| Parse breseq's unassigned junction / missing-coverage evidence (HTML reports) | `analysis/genomics/02_structural_variation/B1. unassigned_evidence.ipynb` (parsing in a raw cell) | `datasets/genomics/02_structural_variation/results/unassigned_*.tsv` |
| Coverage of fixed variants at the endpoints (masked variants) | last cells of `analysis/genomics/03_variant_filtering/B3. background_sequence_coverage.ipynb` (skipped when the coverage files are absent) | `datasets/genomics/03_variant_filtering/intermediate/endpoint_coverages.tsv` |
| Count phage attB junctions in the raw reads | `analysis/genomics/09_phage_fhuA/P0. attB_occupancies.ipynb` | `datasets/genomics/09_phage_fhuA/intermediate/junctions_counts.tsv` |
| Extract fimS reads from the raw reads | `analysis/adhesion/fim_switch/extract_fims.py`, `batch_extract_fims_reads.sh` | `datasets/adhesion/fim_switch/FimS_Reads/` |
| breseq on the marker-frequency samples | `analysis/chromosome_replication/coverage_analysis/run_breseq_marker_frequency.sh` | per-sample coverage tables in `datasets/chromosome_replication/coverage_analysis/intermediate/` |

Note: the VCF files of samples EC04 and EC05 are swapped (sample mix-up); `03_variant_filtering/B2. disentangle_variants.ipynb`
corrects the order of time points 4 and 5 in population C.

### Datasets

Committed files are the starting points of the analysis. A few start downstream of the raw data
because the raw files are too large for a Git repository:

| Committed file | Stands in for |
|---|---|
| `datasets/genomics/03_variant_filtering/intermediate/endpoint_coverages.tsv` | breseq coverage tables of the LTEE endpoints (1.2 GB) |
| `datasets/genomics/02_structural_variation/results/unassigned_*.tsv` | breseq "unassigned evidence" HTML reports (221 MB; parsed by `02_structural_variation/B1. unassigned_evidence.ipynb`) |
| `datasets/chromosome_replication/coverage_analysis/intermediate/` | breseq coverage tables of the marker-frequency samples (1.1 GB) |
| `datasets/morphology/Cells/` | 32 GB of ND2/TIFF phase-contrast images, segmented by the morphology pipeline |
| `datasets/genomics/09_phage_fhuA/intermediate/junctions_counts.tsv` | the raw Illumina reads (counted by `09_phage_fhuA/P0. attB_occupancies.ipynb`) |

Other committed inputs:

- `datasets/genomics/vcf/`: breseq variant calls (VCF) for every sequenced population sample
- `datasets/genomics/reference_genome/`: MG1655 (U00096) and Φ80 (JX871397) sequences and annotations
- `datasets/genomics/external_datasets/`: RegulonDB transcription units and promoters, essential genes
- `datasets/genomics/06_mutators/Data_from_Good_2017/`: mutation trajectories from Good et al. 2017
- `datasets/genomics/10_mutation_traces/lineages.tsv`: hand-written lineage table
- `datasets/growth_curves/*/`: plate-reader exports, one folder per experiment
- `datasets/adhesion/`, `datasets/phage/`, `datasets/stationary_survival/`: plate-reader and plaque-count tables
- `figures/static/`: images made outside these notebooks (phage genome map, super-infection
  montage, replication sketch, soft-agar motility collage)

### Software

`environment.yml` describes a conda-forge environment with the Python and R packages the notebooks use
(R notebooks run on the IRkernel); see the comments at its top for setup. The saved panel files
(`Plots/GG/*.rds`) are ggplot2 4.x objects.

The figures use Noto Sans through R's generic `sans` font. On macOS, R draws `sans` as Helvetica unless told
otherwise; to match the published figures, install Noto Sans and add this to the environment's `Rprofile.site`:

```r
options(bitmapType = "cairo")  # PNG through cairo/fontconfig
local({ f <- file.path(Sys.getenv("HOME"), "Library/Fonts", paste0("NotoSans-", c("Regular", "Bold", "Italic", "BoldItalic"), ".ttf"))
        systemfonts::register_font("sans", f[1], f[2], f[3], f[4]) })  # ggsave's PNG device (ragg)
```
