# KARDEMUMMA

[![PyPI version](https://img.shields.io/pypi/v/kardemumma.svg)](https://pypi.org/project/kardemumma/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22882804.svg)](https://doi.org/10.5281/zenodo.22882804)

**KARDEMUMMA** stands for **K**ey **A**nalysis of **R**eproducible **D**ata for **E**fficient **M**onitoring in **U**nified **M**ass **S**pectrometry **M**ethods and **A**ssays.

This repository contains the Python package for processing and quality-checking targeted mass spectrometry outputs (for example Skyline/OpenSWATH-style exports). The tool is built based on targeted proteomics assay at KTH Royal intitute of technology and Science for Life Laboratory (SciLifeLab), Sweden. The aim of this tool is to provide a simplified analysis pipeline of plasma proteomics as well as bridging research and clinical applications. 

> The repository name has now been updated from `skyline_qc` to `kardemumma` throughout the project.

## Documentation

API reference (HTML and PDF) is attached to each [GitHub Release](https://github.com/thanadol-git/kardemumma/releases). Download and open in a browser:

- [📖 API Documentation (HTML)](https://github.com/thanadol-git/kardemumma/releases/latest/download/kardemumma_api.html)
- [📄 API Documentation (PDF)](https://github.com/thanadol-git/kardemumma/releases/latest/download/kardemumma_api.pdf)

## Installation

You can install the dependencies and set up the environment using [Conda](https://docs.conda.io/en/latest/):

1. Clone the repository:
   ```bash
   git clone https://github.com/thanadol-git/kardemumma.git
   cd kardemumma
   ```

2. Create the Conda environment using the `config.yml` file:
   ```bash
   conda env create -f environment.yml -p ./env
   
   ```

   Alternatively, if you want to use a unique environment name:
   ```bash
   conda env create -f environment.yml -n kardemumma
   ```

3. Activate the environment:
   ```bash
   conda activate ./env
   ```
   or, if you used an environment name:
   ```bash
   conda activate kardemumma
   ```

4. Install the package in editable mode:
   ```bash
   pip install -e .
   ```

## Available pipelines
- [ ] `SDRF generation for plasma proteomics`
- [ ] `Targeted PRM with ProteomEdge AB`
- [ ] `Targeted SRM`


## Requirement

### To run the script
- Result files exported from Skyline (`*.csv`) and a [targeted-SDRF](https://github.com/bigbio/sdrf-templates/issues/44) file (`.sdrf.tsv`) are required for analysis. Column names are validated in `src/kardemumma/importer.py`.

  **Skyline CSV — required columns** (`ImportSkylineFile`, `_SKYLINE_EXPECTED_COLS`):

  `Precursor`, `Replicate`, `File Name`, `Peptide Retention Time`, `Predicted Retention Time`, `Precursor Charge`, `Peptide Sequence`, `Peptide`, `Normalized Area`, `RatioLightToHeavy`, `Ratio Dot Product`, `Library Dot Product`, `Protein Name`

  | Precursor | Replicate | File Name | Peptide Retention Time | Predicted Retention Time | Precursor Charge | Peptide Sequence | Peptide | Normalized Area | RatioLightToHeavy | Ratio Dot Product | Library Dot Product | Protein Name |
  |---|---|---|---:|---:|---:|---|---|---:|---:|---:|---:|---|
  | 547.7823++ (light) | Sample_01 | Sample_01.raw | 12.45 | 12.30 | 2 | LVNELTEFAK | LVNELTEFAK | 1240000 | 0.98 | 0.95 | 0.97 | sp\|P01023\|A2M |
  | 547.7823++ (heavy) | Sample_01 | Sample_01.raw | 12.46 | 12.30 | 2 | LVNELTEFAK | LVNELTEF[K].AK | 1310000 | | | 0.96 | sp\|P01023\|A2M |
  | 490.5796+++ (light) | Sample_02 | Sample_02.raw | 13.35 | 13.83 | 3 | ELDKYGVSDYHK | ELDKYGVSDYHK | 850000 | 1.05 | 0.94 | 0.92 | iRT_Tag |

  **SDRF — required columns for merge** (`MergeFiles.merge_files`):

  - `comment[data file]` (joined to Skyline `File Name`)
  - all `characteristics[...]` columns present in the file
  - optional: any `factor value[...]` columns (merged when present)

  | source name | characteristics[plate] | characteristics[Sample] | comment[data file] | factor value[Sample] |
  |---|---|---|---|---|
  | PRM_20251218_MarthaReselection_Plate_9_A4 | Plate_9 | CHAPS | PRM_20251218_MarthaReselection_Plate_9_A4.raw | CHAPS |
  | PRM_20251218_MarthaReselection_Plate_9_B4 | Plate_9 | CHAPS | PRM_20251218_MarthaReselection_Plate_9_B4.raw | CHAPS |
  | PRM_20251218_MarthaReselection_Pool | Plate_9 | Pool | PRM_20251218_MarthaReselection_Pool.raw | Pool |

### CLI workflow

With a Skyline CSV and SDRF in place, run the three-step pipeline:

```bash
# Step 1 — import, QC summary, and cutoff guidance
kardemumma-preview \
  --skyline path/to/skyline_export.csv \
  --sdrf path/to/experiment.sdrf.tsv \
  --output results/

# Step 2 — apply filters, pool normalization, and batch correction
kardemumma-cutoff \
  --skyline path/to/skyline_export.csv \
  --sdrf path/to/experiment.sdrf.tsv \
  --output results/

# Step 3 — absolute quantification and report plots (reads results/ from step 2)
kardemumma-report \
  --output results/ \
  --sdrf path/to/experiment.sdrf.tsv
```

To generate a ratio-analysis Jupyter notebook instead, use `kardemumma-ipynb`:

```bash
kardemumma-ipynb \
  --skyline path/to/skyline_export.csv \
  --sdrf path/to/experiment.sdrf.tsv \
  --output notebooks/ \
  --group-a CHAPS \
  --group-b Pool
```

Run any command with `--help` for optional parameters (e.g. `--dotp`, `--light-cutoff`, `--pool-value`).

### Notes
- All required dependencies will be installed via Conda and pip as specified in `environment.yml`.
- Python 3.10 is recommended.
- For pip installs, make sure you have internet access.

## To Dos

### Version 0.1.3
- [ ] Unit testing 
- [ ] mzQC
- [ ] pmultiqc
- [ ] Logo banners

### Phase 1 — Python Package & PyPI Release (Version 0.x.x)

**1. Code & API clean-up**
- [x] -1. Move plot functions from DA4K notebook to `prm.py` and expose via `kdm.*` (`map_peptide_sequence`, `plot_peptide_concentration_by_group`, `plot_median_peptide_concentration_by_group`, `plot_all_median_peptide_concentration_by_group`, `plot_all_peptide_concentration_by_group`)
- [ ] 0. Remove 3 under-QC samples from analysis
- [ ] 1. Check with Yasset on how to set up targeted SDRF
- [ ] 2. Integrate `prm-slider` to work with transition levels
- [ ] 3. Work on SRM support (`sdrf.py` + new `srm.py` module)
- [ ] 4. Combine output layer with OpenMS formats
- [ ] 5. Audit all public functions — consistent naming, type hints, docstrings
- [ ] 6. Ensure `__init__.py` exports a clean, stable public API
- [ ] Create landing logo and banners.
- [ ] Remove 3 proejct ipynb

**2. Package metadata & build**
- [x] 7. Update `pyproject.toml`: add missing `pyteomics` dependency, bump version, add classifiers (`python_requires`, `install_requires`) ✔️ (done)
- [x] 8. Add `CHANGELOG.md` with initial release notes (included in release & GitHub Action)
- [x] 9. Add `LICENSE` file (MIT)
- [x] 10. Verify `pip install -e .` builds cleanly in a fresh environment
- [x] 11. Build distribution: `python -m build` → inspect `dist/`

**3. Testing & CI**
- [x] 12. Add unit tests with `pytest` for core modules (`prm.py`, `sdrf.py`)
- [x] 13. Add a GitHub Actions workflow (`.github/workflows/release.yml`) for releases (tests run on push/tag)
- [x] 14. Add a release workflow that publishes docs on version tag push (`.github/workflows/release_docs.yml`); PyPI publish pending

**4. PyPI release**
- [x] 15. Register package name on [PyPI](https://pypi.org) (check availability of `kardemumma`)
- [x] 16. Create API token on PyPI and store as `PYPI_API_TOKEN` GitHub secret
- [x] 17. Publish first release: `python -m twine upload dist/*` (or via GitHub Actions)
- [x] 18. Verify: `pip install kardemumma` works from PyPI

---

### Phase 2 — Nextflow Pipeline (Version 2.x.x)

**5. Pipeline design**
- [ ] 19. Define end-to-end workflow: raw input → SDRF validation → OpenSWATH/Skyline export → PRM QC → ratio/DA output
- [ ] 20. Sketch module boundaries as Nextflow `process` blocks (one process per major step)
- [ ] 21. Decide on container strategy: Docker images (or Singularity) per process, each with `kardemumma` installed from PyPI

**6. Implementation**
- [ ] 22. Scaffold repository structure: `nextflow/`, `modules/`, `conf/`, `assets/`
- [ ] 23. Write a `main.nf` entry workflow with configurable params (`--input`, `--outdir`, `--mode prm|srm`)
- [ ] 24. Implement individual processes wrapping `kardemumma` CLI calls or Python scripts
- [ ] 25. Add `nextflow.config` with profiles: `standard` (local), `cluster` (SLURM/HPC at SciLifeLab), `cloud`
- [ ] 26. Pin `kardemumma` version in each container/environment to match tested PyPI release

**7. Testing & docs**
- [ ] 27. Add small test dataset (synthetic or anonymised) to `tests/` for end-to-end pipeline testing
- [ ] 28. Add `nf-test` or a simple CI job that runs the pipeline on the test dataset
- [ ] 29. Write pipeline usage docs in `docs/pipeline.md` (input format, params, outputs)
- [ ] 30. Consider submission to [nf-core](https://nf-co.re) once pipeline is stable


## Issues
- iRT peptides: why do they contain Biognosys sequences?
- Oxidation
- Stats for PEP

## Key developers
- Thanadol Sutantiwanichkul
- Justin Sing
- Yuqi Zheng
- Khue Hua Tran Minh
- Maria-Jesus Iglesias Mareque
- Fredrik Edfors