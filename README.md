# skyline_qc

This will be a simple Python package to clean up Skyline results. 

## Installation

You can install the dependencies and set up the environment using [Conda](https://docs.conda.io/en/latest/):

1. Clone the repository:
   ```bash
   git clone https://github.com/thanadol-git/skyline_qc.git
   cd skyline_qc
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

### Notes
- All required dependencies will be installed via Conda and pip as specified in `environment.yml`.
- Python 3.10 is recommended.
- For pip installs, make sure you have internet access.

## To Dos
1. PRM tool from MARTHA project
2. Work a bit with SRM 
3. Combine with OpenMS 
