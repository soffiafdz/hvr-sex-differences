# HVR Sex Differences in UK Biobank

Code and manuscript source for:

> Fernandez-Lozano S, Collins LD. **Hippocampal-to-Ventricle Ratio as a Head-Size-Independent Biomarker: Sex Differences and Cognitive Associations in 27,680 UK Biobank Participants.** *Human Brain Mapping* (submitted).

## Overview

This project examines whether the hippocampal-to-ventricle ratio (HVR), a self-normalizing brain metric, provides more consistent sex difference estimates than conventional head-size adjustment methods. Using UK Biobank structural MRI data (N = 27,680), we:

1. Compare hippocampal sex differences across four adjustment methods
2. Fit GAMLSS normative centile curves stratified by sex
3. Evaluate brain-cognition associations via structural equation modeling
4. Provide sex-specific normative reference tables for clinical use

## Data Availability

This study uses data from the [UK Biobank](https://www.ukbiobank.ac.uk/) (Application Number 35605). Individual-level data cannot be shared publicly under the UK Biobank Material Transfer Agreement. Researchers can apply for access at <https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access>.

The `data/` directory is not included in this repository. To reproduce the analyses, you will need:

- UK Biobank parquet dataset with the fields listed in `data/raw/fields.yaml`
- Preprocessed brain volumes from the UK Biobank imaging pipeline (`data/raw/pp_vols/`)
- AssemblyNet segmentation volumes (`data/fst/ukb_assemblynet.fst`)
- DARQ quality control scores (`data/lists/darq_out_*.csv`)

## Directory Structure

```
hvr_sex_differences/
├── R/
│   ├── 00_run_pipeline.R          # Master orchestration script
│   ├── 01_parse_volumes.R         # Aggregate preprocessed volumes
│   ├── 02_quality_control.R       # DARQ + AssemblyNet QC filtering
│   ├── 03_parse_covariates.R      # Extract UK Biobank covariates
│   ├── 04_clean_ses_cognitive.R   # Clean SES and cognitive data
│   ├── 05_adjust_headsize.R       # Head-size adjustment + matching
│   ├── 06_cognitive_factors.R     # CFA for cognitive factors
│   ├── 07_demographics.R          # Demographic tables
│   ├── 08_normative_tables.R      # GAMLSS normative modeling
│   ├── 09_sex_differences.R       # Sex difference analysis
│   ├── 10_hvr_comparison.R        # HVR vs HC comparison
│   ├── 11_sem_analysis.R          # Structural equation modeling
│   ├── 12_manuscript_objects.R    # Pre-compute manuscript environment
│   └── utils/                     # Shared utility modules
│       ├── config.R               #   Configuration management
│       ├── data_io.R              #   Safe I/O wrappers
│       ├── export.R               #   Publication figure/table export
│       ├── formatting.R           #   Numeric formatting helpers
│       ├── gamlss.R               #   GAMLSS utilities
│       ├── logging.R              #   Structured logging
│       ├── plotting_core.R        #   Palettes, themes, save_plot
│       ├── plotting_figures.R     #   Publication figure functions
│       ├── plotting_pipeline.R    #   Exploratory pipeline plots
│       ├── statistics.R           #   Statistical helpers
│       ├── tables_core.R          #   Table styling and LaTeX fixes
│       ├── tables_data.R          #   Data extraction for tables
│       ├── tables_gt.R            #   GT table formatting
│       ├── tables_normative.R     #   GAMLSS normative tables
│       └── validation.R           #   Data validation
├── config/
│   └── pipeline_config.yaml       # All paths, parameters, and settings
├── reports-src/
│   ├── manuscript.qmd             # Main manuscript (Quarto)
│   ├── supplementary.qmd          # Supplementary material
│   ├── _common.R                  # Shared Quarto setup
│   ├── _quarto.yml                # Quarto project config
│   ├── title.tex                  # LaTeX title template
│   ├── references.bib             # Bibliography
│   ├── human-brain-mapping.csl    # HBM citation style
│   └── apa.csl                    # APA citation style
├── renv.lock                      # R package versions
├── environment.yml                # Conda environment
└── LICENSE                        # MIT License
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

Key dependencies: `data.table`, `lavaan`, `gamlss`, `gt`, `ggplot2`, `duckdb`, `fst`, `MatchIt`.

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

### Individual scripts

```bash
Rscript R/00_run_pipeline.R --scripts=5,9    # Run only steps 5 and 9
Rscript R/00_run_pipeline.R --from=7         # Run from step 7 onward
Rscript R/00_run_pipeline.R --force          # Force regenerate all outputs
```

Scripts can also be run directly (e.g., `Rscript R/09_sex_differences.R`), though they expect upstream outputs to exist.

### Rendering the manuscript

First run step 12 to pre-compute manuscript objects, then render:

```bash
Rscript R/12_manuscript_objects.R
cd reports-src && quarto render
```

Rendered output goes to `outputs/reports/`.

## License

MIT License. See [LICENSE](LICENSE).
