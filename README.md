# NULISA proteomics analysis: serum and CSF biomarkers across genetic groups and treatments

Analysis code for NULISAseq (Alamar Biosciences) protein measurements in serum and cerebrospinal fluid
(CSF). It covers differential protein abundance between patient groups and healthy controls,
longitudinal changes under two treatments, multivariate analyses (PCA, k-means) and correlations of
protein levels with a clinical score.

| Step | File | Language | What it does |
|---|---|---|---|
| 1 | `analysis/01_differential_abundance_nulisaseqr.R` | R | NULISAseqR linear models (`lmNULISAseq`) for group comparisons and DrugB time points vs baseline; volcano plots (FDR and unadjusted) |
| 2a | `analysis/02a_druga_serum_csf_analysis_v1.ipynb` | Python | DrugA cohort, serum + CSF: baseline vs healthy controls, DrugA vs baseline, PCA / k-means, Venn diagrams, deltas (version 1) |
| 2b | `analysis/02b_druga_serum_csf_analysis_v2.ipynb` | Python | Same analysis as 2a (version 2) |
| 2c | `analysis/02c_druga_serum_csf_score_correlations.ipynb` | Python | DrugA cohort: group comparisons plus protein–clinical-score correlations |
| 3 | `analysis/03_csf_group_comparison_longitudinal.ipynb` | Python | CSF: genetic groups and longitudinal DrugA / DrugB samples, PCA, observed-vs-expected correlations |
| 4 | `analysis/04_genetic_groups_analysis.ipynb` | Python | Genetic groups and the DrugB time course: PCA / k-means, correlations with the clinical score |
| 5a | `analysis/05a_correlation_plots.ipynb` | Python | Correlation plots of protein levels vs clinical score (all data points), serum and CSF |
| 5b | `analysis/05b_correlation_plots_high_nfh_ab42.ipynb` | Python | As 5a, including high-NfH / Aβ42 subsets |

## Repository structure

```
├── analysis/                  # 01-05 (see table above)
├── config/config.example.yaml # template: copy to config.yaml (git-ignored) and edit
├── utils/                     # shared config loaders (R, Python)
├── envs/                      # conda environment + R package installer
└── data/                      # created locally, git-ignored
```

## Group labels

All code uses generic labels. Sample names and column names in your input files must use the same labels.

| Label | Meaning |
|---|---|
| `Control`, `HC` | Healthy controls |
| `Case` | All patients combined |
| `NoMut` | Patients without a known causal mutation |
| `MutBpos` / `MutB` | Carriers of mutation B (`MutB-…` marks samples from the DrugB study, e.g. `MutB-Baseline`) |
| `MutCpos` / `MutC` | Carriers of mutation C (e.g. `MutCpos-Baseline`) |
| `DrugA`, `DrugB` | Treatments; time points as e.g. `DrugB-3month`, `DrugA-12months` |
| `BL` | Baseline |
| `Score` | Clinical severity score |

The code selects samples by substring, e.g. `grepl("MutB", …)` matches both `MutBpos…` and `MutB-…`.
Keep this naming scheme so that each sample is assigned to the right group.

## Setup

```bash
git clone https://github.com/<user>/NULISA.git
cd NULISA

cp config/config.example.yaml config/config.yaml   # then edit: data folders + excluded samples
conda env create -f envs/nulisa_python.yml
Rscript envs/install_r_packages.R                  # includes NULISAseqR from GitHub
```

`config/config.yaml` is git-ignored. It holds the data folders and the names of samples excluded from
selected comparisons, so that real sample names never appear in the code.

## Input files

Each analysis reads its inputs from its own folder, set in the config.

| Config key | Used by | Files expected in that folder |
|---|---|---|
| `de_workdir` | 01 | `Raw Data.csv`: NULISAseq long-format export (`SampleName`, `Target`, `NPQ`) |
| `druga_workdir` | 02a–c | `Data_serum_filtered.csv`, `Data_CSF_filtered.csv`, `Data_*_filtered_BL_Ctrl.csv`, `Data_*_filtered_BL_DrugA.csv`, `Data_*_sigs_BL_HC.csv`, `Serum sigs BL DrugA.csv`, `CSF sigs BL DrugA.csv`; 02a/b also `delta serum.csv`, `Delta_Serum_data2.csv`, `Delta_Serum_sigs.csv`, `Marker_Score_ratio.csv`, `Prog_CSF_Protein.csv`, `Prog_Score_Serum.csv`; 02c also `Serum_Score_Corr.csv`, `CSF_Score_Corr.csv` |
| `longitudinal_workdir` | 03 | `Data_final_analysis.csv`, `Expected_Corr.xlsx`, `Corr_12months.xlsx`, `NPTX 3months Corr.xlsx` |
| `genetic_groups_workdir` | 04 | `Data_final_analysis.csv`, `NoMut_Correlations.xlsx`, `MutBpos_Correlations.xlsx`, `Case_Corr.xlsx`, `GeoMean_Corr.xlsx`, `Geomean_corr_Case_patients.xlsx` |
| `correlation_workdir` | 05a–b | `correlation.xlsx`, `Serum_Score_Corr.csv`, `CSF_Score_Corr.csv`; 05b also `correlation_highNfH_Aß42.xlsx`, `correlation281124.xlsx`, `correlation031224 male.xlsx` |

Outputs (tables and figures) are written to the same folder.

## Running the analysis

Run from the repository root:
```bash
Rscript analysis/01_differential_abundance_nulisaseqr.R
```
Open the notebooks with Jupyter (environment `nulisa_python`) and run them top to bottom. The first cell of
each notebook loads the config and switches to its data folder. Notebooks 02a–c write to the same folder,
so later runs overwrite earlier outputs with the same names.

**Statistics.** The R script reports FDR-adjusted, unadjusted and Bonferroni p-values. Volcano plots are
drawn for both FDR-adjusted and unadjusted p < 0.05.

## Data protection

This repository contains code only. Patient data, sample names, clinical scores and all results stay
local: `data/`, `config/config.yaml` and all tabular and image files are git-ignored. Notebooks are stored
without outputs; keep it that way before committing (`pip install nbstripout && nbstripout --install`).

## Reproducibility

Record package versions before publishing: `renv::snapshot()` (R) and
`conda env export --no-builds -n nulisa_python > envs/nulisa_python.lock.yml`.

## Citation

If you use this code, please cite: *<paper reference / DOI>*.

## License

MIT, see [LICENSE](LICENSE).
