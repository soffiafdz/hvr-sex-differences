# HVR Sex Differences in UK Biobank

Code and manuscript source for:

> Fernandez-Lozano S, Collins LD. **Hippocampal-to-Ventricle Ratio as a Head-Size-Independent Biomarker: Sex Differences and Cognitive Associations in 27,680 UK Biobank Participants.** *Human Brain Mapping* (submitted).

## Overview

This project examines whether the hippocampal-to-ventricle ratio (HVR), a self-normalizing brain metric, provides more consistent sex difference estimates than conventional head-size adjustment methods. Using UK Biobank structural MRI data (N = 27,680), we:

1. Compare hippocampal sex differences across four head-size adjustment methods (unadjusted, proportions, stereotaxic, residualized) and in samples matched on age and ICV
2. Fit GAMLSS normative centile curves stratified by sex for hippocampus, lateral ventricles, and HVR
3. Evaluate brain-cognition associations via multi-group structural equation modeling
4. Provide sex-specific normative reference tables for clinical use

## Data Availability

This study uses data from the [UK Biobank](https://www.ukbiobank.ac.uk/) (Application Number 35605). Individual-level data cannot be shared publicly under the UK Biobank Material Transfer Agreement. Researchers can apply for access at <https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access>.

The `data/`, `models/`, and `outputs/` directories are not included in this repository. To reproduce the analyses, you will need:

| Input | Path (in config) | Description |
|-------|-------------------|-------------|
| UK Biobank parquet | `data/raw/ukb676574` | Main phenotype dataset |
| Field list | `data/raw/fields.yaml` | Selected UKB field IDs |
| Field metadata | `data/raw/fields.tsv` | UKB field definitions |
| Preprocessed volumes | `data/raw/pp_vols/*.csv` | Per-subject ICC volumes and scale factors |
| AssemblyNet | `data/fst/ukb_assemblynet.fst` | AssemblyNet segmentation volumes |
| DARQ QC | `data/lists/darq_out_*.csv` | Deep-learning QC scores |
| Exclusion list | `data/lists/exclude_id.txt` | Subject IDs to exclude |
| Segmentations | `data/derivatives/hclvag_segmentations.csv` | HC/LV segmentation volumes |
| lavaan models | `models/lavaan/*.txt` | SEM/CFA model specifications |

All data paths are configured in `config/pipeline_config.yaml` under the `data:` section.

## Pipeline

### Dependency graph

```
01 Parse Volumes ─────┐
                      ├──→ 05 Adjust Head-size ──┬──→ 08 Normative Tables
02 Quality Control ───┘          │                ├──→ 09 Sex Differences
                                 │                └──→ 10 HVR Comparison ──┐
03 Parse Covariates ──→ 04 Clean SES/Cognitive    │                        │
                                 │                │                        │
                                 └──→ 06 Cognitive Factors ───────────────┤
                                           │                               │
                                           └──→ 11 SEM Analysis ─────────┤
                                                                          │
                        07 Demographics ──────────────────────────────────┤
                                                                          │
                                         12 Manuscript Objects ◄──────────┘
                                                   │
                                           Quarto Render
```

### Script summary

| Step | Script | Description | Key outputs |
|------|--------|-------------|-------------|
| 00 | `00_run_pipeline.R` | Orchestrator: loads utils, parses CLI args, runs steps in order | — |
| 01 | `01_parse_volumes.R` | Aggregate per-subject brain volumes from CSV files | `ukb-lng_icc-scale.csv` |
| 02 | `02_quality_control.R` | Combine DARQ and AssemblyNet QC scores; apply thresholds | `qc_darq-assemblynet.rds` |
| 03 | `03_parse_covariates.R` | Query UK Biobank parquet for demographics, ICD-10, education | `ukb_covars.fst` |
| 04 | `04_clean_ses_cognitive.R` | Clean SES indices, parse and score cognitive battery | `cog-tests.rds` |
| 05 | `05_adjust_headsize.R` | Apply 4 head-size adjustments; optimal sex-age-ICV matching | `hc-hvr_adj.rds` |
| 06 | `06_cognitive_factors.R` | Bifactor CFA; measurement invariance; test-retest reliability | `lat-cog_vals.rds` |
| 07 | `07_demographics.R` | Demographic comparison tables (primary, matched, sensitivity) | `demog_data.rds` |
| 08 | `08_normative_tables.R` | GAMLSS normative centile modeling (BCCG/NO families, site covariate) | `gamlss.rds`, centile tables |
| 09 | `09_sex_differences.R` | Cohen's d effect sizes across methods, samples, age bins | `sex_differences.rds` |
| 10 | `10_hvr_comparison.R` | HVR vs HC: age-sex interactions, incremental R-squared | `hvr_hc_comparison.rds` |
| 11 | `11_sem_analysis.R` | Multi-group SEM (HC, HVR, HC-residualized) with bootstrap | `sem_analysis.rds` |
| 12 | `12_manuscript_objects.R` | Assemble all results into `manuscript_env.rds` for Quarto | `manuscript_env.rds` |

### Generated outputs

The pipeline writes to three directories (all gitignored):

- **`data/derivatives/`** — Intermediate RDS/FST/CSV files (cleaned data, factor scores, effect sizes)
- **`models/`** — Fitted model objects (GAMLSS, lavaan CFA/SEM) and parameter extractions
- **`outputs/`** — Publication figures (PNG), tables (HTML/TEX), and rendered documents (PDF)

## Directory Structure

```
hvr_sex_differences/
├── R/
│   ├── 00_run_pipeline.R          # Master orchestration script
│   ├── 01–12_*.R                  # Pipeline steps (see table above)
│   └── utils/                     # Shared utility modules
│       ├── config.R               #   YAML config loader
│       ├── data_io.R              #   Safe read/write wrappers
│       ├── export.R               #   Publication figure/table export
│       ├── formatting.R           #   Numeric formatting (p, d, CI)
│       ├── gamlss.R               #   GAMLSS prediction and validation
│       ├── logging.R              #   Structured logging system
│       ├── plotting_core.R        #   Palettes, theme, save_plot
│       ├── plotting_figures.R     #   Publication figures (fig_*)
│       ├── plotting_pipeline.R    #   Exploratory plots (plot_*)
│       ├── statistics.R           #   Williams t-test, partial cor
│       ├── tables_core.R          #   GT styling, LaTeX post-processing
│       ├── tables_data.R          #   Data extraction for tables
│       ├── tables_gt.R            #   GT table builders
│       ├── tables_normative.R     #   GAMLSS normative table builders
│       └── validation.R           #   Column/range/factor validation
├── config/
│   └── pipeline_config.yaml       # Paths, parameters, settings
├── reports-src/
│   ├── manuscript.qmd             # Main manuscript (Quarto)
│   ├── supplementary.qmd          # Supplementary tables and figures
│   ├── _common.R                  # Shared setup (gt override, refs)
│   ├── _quarto.yml                # Quarto project config
│   ├── title.tex                  # LaTeX title partial
│   ├── references.bib             # Bibliography
│   ├── human-brain-mapping.csl    # HBM citation style
│   └── apa.csl                    # APA citation style
├── renv.lock                      # R package versions (renv::restore)
├── environment.yml                # Conda environment spec
└── LICENSE                        # MIT
```

## Prerequisites

- **R** >= 4.3
- **Quarto** >= 1.4 (for manuscript rendering)
- **Conda/Mamba** (optional, for environment management)

### R packages

All R package versions are pinned in `renv.lock`. To restore:

```r
install.packages("renv")
renv::restore()
```

Key dependencies: `data.table`, `lavaan`, `gamlss`, `gt`, `ggplot2`, `duckdb`, `fst`, `MatchIt`, `patchwork`, `ggtext`.

### Conda environment (optional)

```bash
conda env create -f environment.yml
conda activate hvr_sex_differences
```

## Running the Pipeline

### Full pipeline

```bash
Rscript R/00_run_pipeline.R
```

### Selective execution

```bash
Rscript R/00_run_pipeline.R --scripts=5,9    # Run only steps 5 and 9
Rscript R/00_run_pipeline.R --from=7         # Run from step 7 onward
Rscript R/00_run_pipeline.R --to=6           # Run steps 1-6
Rscript R/00_run_pipeline.R --force          # Force regenerate all outputs
```

Individual scripts can be run directly (e.g., `Rscript R/09_sex_differences.R`), though they expect upstream outputs to exist.

### Rendering the manuscript

The pipeline orchestrator runs step 12 automatically. To render manually:

```bash
Rscript R/12_manuscript_objects.R           # Pre-compute inline values
cd reports-src && quarto render             # Render PDF + supplementary
```

Output goes to `outputs/reports/`.

## Configuration

All parameters are centralized in `config/pipeline_config.yaml`:

- **`general`** — Random seed (1618), log level, core count
- **`data`** — Input/output file paths (all relative to project root)
- **`parameters`** — QC thresholds, ICD-10 exclusion patterns, cognitive test definitions, CFA/SEM estimators, GAMLSS families, adjustment methods, ROI definitions
- **`scripts`** — Per-script settings (force regeneration flags, matching tolerances, bootstrap config)
- **`manuscript`** — Figure dimensions and table specifications

## License

MIT License. See [LICENSE](LICENSE).
