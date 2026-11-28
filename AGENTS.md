# Guide for coding agents

Code and data behind the figures of "Perpetual exponential growth". Read `README.md` for the layout, the
run order and the provenance of the datasets; this file is about changing a figure.

## How a figure is made

1. An analysis notebook (`analysis/<topic>/...`) draws each panel and saves it with `psave('Plots/<name>')`,
   which writes `Plots/<name>.png`, `Plots/PDF/<name>.pdf` and the ggplot object `Plots/GG/<name>.rds`.
2. `figures/make_collage.ipynb` reads the `.rds` panels, arranges them with patchwork (panel tags, sizes)
   and saves the figure to `figures/output/<figure>.png` (plus `output/PDF/`).
3. Figure 1 only: `figures/make_adhesion.py` draws panels A–C and stacks them above panels D–G
   (`output/adhesion_DG.svg`, from the collage).

To change how a panel looks, edit the analysis notebook that draws it, re-run that notebook, then re-run
the setup cells (first two) and the figure's section of `make_collage.ipynb`. To change the layout of a figure
(panel sizes, order, tags), edit only `make_collage.ipynb`. Each notebook only depends on files written by
the notebooks above it in the README's run order.

## Where each figure comes from

| Figure (`figures/output/`) | Panels drawn in |
|---|---|
| `adhesion.png` (Fig. 1) | A–C: `figures/make_adhesion.py` (simulations, schematics; C from `datasets/adhesion/biofilm_OD_traces.json`). D adhesion, E fimS, G doubling time: `figures/make_sticking_figures.ipynb` (→ `figures/panels/*.rds`). F fimH: `analysis/genomics/10_mutation_traces/mutation_traces.ipynb` |
| `genomics.png` (Fig. 2) | `10_mutation_traces/experiment_overview.ipynb`, `10_mutation_traces/mutation_traces.ipynb` (takeover, domestication scars), `06_mutators/F1. mutators.ipynb` (allele-frequency plots), `07_mutation_statistics/G1. mutation_statistics.ipynb` (enrichments) — all under `analysis/genomics/` |
| `latent_traits.png` (Fig. 3) | `analysis/stationary_survival/ss_beta_plots.ipynb`, `analysis/adhesion/motility_summary/motility_summary.ipynb`, `analysis/growth_curves/2024-06-11_beta/Esch_beta_plots_pure_only.ipynb`, `analysis/morphology/plots_phase.ipynb` |
| `phage.png` (Fig. 4) | `analysis/genomics/09_phage_fhuA/P1. fhuA.ipynb` (A), `09_phage_fhuA/P2. phage_dynamics.ipynb`, `analysis/phage/phage_production/plots_esch_AF_phage_production.ipynb`, `analysis/phage/phage_collapse/phage_collapse_plot.ipynb`; image `figures/static/montage_horiz.png`; F genome map: `analysis/phage/phage_identification/phage_genome_explorer.ipynb` (pygenomeviz, exported to `figures/static/genome_plot.png`) |
| `chromosome_replication.png` (Fig. 5) | `analysis/chromosome_replication/coverage_analysis/O3. plots_MF_epsilon.ipynb`; image `figures/static/multifork_replication_sketch.png` |
| `supp_mutators.png` | `G1. mutation_statistics.ipynb` (A), `F1. mutators.ipynb` (B) |
| `supp_structural.png` | `O3. plots_MF_epsilon.ipynb` (inversion in F), `analysis/genomics/02_structural_variation/B2. unassigned_plots.ipynb` (deletions) |
| `supp_phage.png` | `analysis/phage/C20-21_resistance/plots.ipynb`, `mutation_traces.ipynb` (fhuA vs mutators), `P2. phage_dynamics.ipynb` |
| `separate_dynamics.png` | the `*_plots` notebooks in `analysis/growth_curves/2025-03-23_AF_gamma/`, `2025-04-05_BC_alpha/`, `2025-04-15_DE_alpha/` |

Figures saved directly by an analysis notebook (no collage): `analysis/genomics/08_convergence/Plots/multiplicity_vs_null.png`
(`N0. null_model.ipynb`, Python), `analysis/genomics/07_mutation_statistics/Plots/prom_mut_plot.png` (`G1`),
`analysis/growth_curves/2025-03-23_AF_gamma/Plots/growth_curve_reproductibility.png`.

To find the producer of any panel, search the notebooks for the `.rds` name that `make_collage.ipynb` reads.

## Running things

- Environment: `environment.yml` (see its header). R notebooks use the IRkernel, Python notebooks a Python 3.11
  kernel; Noto Sans must be installed (README, "Software").
- Run a notebook headlessly from the repo root:
  `jupyter nbconvert --to notebook --execute --inplace "path/to/notebook.ipynb"` (it runs in the notebook's
  folder; all paths are relative to it).
- Figure 1: `cd figures && python make_adhesion.py` after `make_collage.ipynb` (needs `rsvg-convert`).
- Commit notebooks with their outputs (they document the results). Generated tables and plots are
  gitignored, except the reference figures listed in `.gitignore`; after changing a figure, regenerate its
  PNG so the committed copy matches the code.
- Look at the resulting PNG to check a change: panels share a fixed figure size, so text or legends can
  overflow or overlap silently.

## Conventions

- R notebooks start with `source('../../R/preamble.R')` (depth varies): `size(h, w)` sets the plot size in
  inches, `psave()` saves; the ggplot theme and palettes are in `R/theme.R`. In `make_collage.ipynb`,
  `letters(top, right, bottom, left)` (setup cell) styles and positions the panel tags.
  `ggplot2` 4.x objects are saved as `.rds`, so panels are restyled in the collage with `+ theme(...)`.
- Generations are often stored as time points: `gen` × 2,500 = generations.
- Sequencing sample names are `<experiment><population><generation>`, e.g. `EB21` = experiment E,
  population B, generation 21. The first letter is the experiment, *not* a population: `experiment == 'E'`
  is the main 24-point time course and keeps all eight populations, while the `A*` and `I*` samples are
  separate one-off sets. The population is the second letter. Do not mistake a filter on `experiment`
  for a filter on `pop`.
- `figures/make_adhesion.py` writes SVG by hand, in points (figure width 504 pt = 7 in), with
  styling constants at the top. Text is placed using widths measured for Noto Sans (constants named `W_*`,
  `LETTER_W`, `WL`): if you change a label's text or size, update its width too.
- Committed files under `datasets/` are inputs (raw or the earliest distributable step); never edit them
  to change a figure.
