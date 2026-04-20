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
   conda env create -f environment.yml -n skyline-qc
   ```

3. Activate the environment:
   ```bash
   conda activate ./env
   ```
   or, if you used an environment name:
   ```bash
   conda activate skyline-qc
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

### Notes
- All required dependencies will be installed via Conda and pip as specified in `environment.yml`.
- Python 3.10 is recommended.
- For pip installs, make sure you have internet access.

### Storage location for Edfors lab, 
The project should be located at the `hot storage` of the lab. One can find the project name inside. Within the folder, there shall be a raw folder where we keep the raw files from MS injections. The SDRF file should be located along side. Please take a look below. Constatnly, the raw file should be tested for it completeness.
```
kardemumma/
└── Project ABC/
    ├── raw/
    │   ├── sample1.raw
    │   ├── sample2.raw
    │   ├── sample3.raw
    │   └── ... (multiple raw files)
    ├── Project_ABC.sdrf.tsv
    └── Others/
```




## To Dos
0. Remove 3 under-QC samples from analysis 
1. Check with Yasset on how to set up targeted SDRF.
2. Publish python package asap
3. Develop snakemake pipeline 
4. Combine with OpenMS.
5. Work a bit with SRM.


## Issues
- iRT peptides: why do they contain Biognosys sequences?
- Oxidation
- Stats for PEP

## Key developers
- Thanadol Sutantiwanichkul
- Justin Sing
- Khue Hua Tran Minh
- Maria-Jesus