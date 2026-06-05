# KARDEMUMMA

**KARDEMUMMA** stands for **K**ey **A**nalysis of **R**eproducible **D**ata for **E**fficient **M**onitoring in **U**nified **M**ass **S**pectrometry **M**ethods and **A**ssays.

This repository contains the Python package for processing and quality-checking targeted mass spectrometry outputs (for example Skyline/OpenSWATH-style exports). The tool is built based on targeted proteomics assay at KTH Royal intitute of technology and Science for Life Laboratory (SciLifeLab), Sweden. The aim of this tool is to provide a simplified analysis pipeline of plasma proteomics as well as bridging research and clinical applications. 

> The repository name has now been updated from `skyline_qc` to `kardemumma` throughout the project.

## Installation

You can install the dependencies and set up the environment using [Conda](https://docs.conda.io/en/latest/):

1. Clone the repository:
   ```bash
   git clone https://github.com/thanadol-git/kardemumma.git
   cd kardemumma
   ```

2. Create the environment using the provided `environment.yml` file:
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
   ```

## Available pipelines
- [ ] `SDRF generation for plasma proteomics`
- [ ] `Targeted PRM with ProteomEdge AB`
- [ ] `Targeted SRM`


## Requirement

### Notes
- All required dependencies will be installed via Conda and pip as specified in `environment.yml`.
- Python 3.10 is recommended.
- For pip installs, make sure you have internet access.

### Storage location for Edfors lab, 
The project should be located at the `hot storage` of the lab. One can find the project name inside. Within the folder, there shall be a raw folder where we keep the raw files from MS injections. The SDRF file should be located along side. Please take a look below. Constatnly, the raw file should be tested for it completeness. (TBD)

```
01_hot/
├── Project_ABC/
│   ├── raw/
│   │   ├── sample1.raw
│   │   ├── sample2.raw
│   │   ├── sample3.raw
│   │   └── ... (multiple raw files)
│   ├── Project_ABC.sdrf.tsv
│   └── Others/
└── Project_XYZ/
└── Project_XYA/
└── Project_XYB/
```




## To Dos

### Phase 1 — Python Package & PyPI Release

**1. Code & API clean-up**
- [ ] 0. Remove 3 under-QC samples from analysis
- [ ] 1. Check with Yasset on how to set up targeted SDRF
- [ ] 2. Integrate `prm-slider` to work with transition levels
- [ ] 3. Work on SRM support (`sdrf.py` + new `srm.py` module)
- [ ] 4. Combine output layer with OpenMS formats
- [ ] 5. Audit all public functions — consistent naming, type hints, docstrings
- [ ] 6. Ensure `__init__.py` exports a clean, stable public API

**2. Package metadata & build**
- [ ] 7. Update `pyproject.toml` (or `setup.cfg`): version, description, classifiers, `python_requires`, `install_requires`
- [ ] 8. Add `CHANGELOG.md` with initial release notes
- [ ] 9. Add `LICENSE` file if missing
- [ ] 10. Verify `pip install -e .` builds cleanly in a fresh environment
- [ ] 11. Build distribution: `python -m build` → inspect `dist/`

**3. Testing & CI**
- [ ] 12. Add unit tests with `pytest` for core modules (`prm.py`, `openswath.py`, `sdrf.py`, etc.)
- [ ] 13. Add a GitHub Actions workflow (`.github/workflows/ci.yml`) that runs tests on push/PR
- [ ] 14. Add a release workflow that publishes to PyPI on version tag push

**4. PyPI release**
- [ ] 15. Register package name on [PyPI](https://pypi.org) (check availability of `kardemumma`)
- [ ] 16. Create API token on PyPI and store as `PYPI_API_TOKEN` GitHub secret
- [ ] 17. Publish first release: `python -m twine upload dist/*` (or via GitHub Actions)
- [ ] 18. Verify: `pip install kardemumma` works from PyPI

---

### Phase 2 — Nextflow Pipeline

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
- Khue Hua Tran Minh
- Maria-Jesus Iglesias Mareque
- Fredrik Edfors