These figures are plotted with Ubuntu 18.04 Linux distribution.
Thus, these insctructions are for Linux operating system.


# Prerequisites


- Python run within Anaconda environment (downloaded from https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh)
- Conda environment file in this repository: [environment.yml](environment.yml)
- Change `prefix` line in [environment.yml](environment.yml) file to be in line with your conda environment (here it's `data`).
- Conda environment can be installed after installation of Anaconda with command:
    - `conda env create -f environment.yml`
- Activate conda environment
    - `conda activate data`

- Set environment variable `${SIMULATIONDATAROOTFOLDER}` with folder that you have stored all the individual simulation folder:
    - `export SIMULATIONDATAROOTFOLDER=/folder/that/contains/all/the/simulations`

- Set environment variable `${SIMULATIONFIGUREFOLDER}` with the folder that you want to store your figures:
    - `export SIMULATIONFIGUREFOLDER=/folder/where/you/want/to/store/the/figures`

- Go to directory of your liking where you can download needed source code repositories.

- Download plotting library:
    - `git clone https://github.com/JaakkoAhola/LES-03plotting`
    - `cd LES-03plotting`
    -  Checkout latest release tag. (For acp-2019-1182 manuscript, checkout release tag: v1.1), e.g.:
        - `git checkout v1.1`


# Prepare metadata for plotting


- Prepare  .csv data file about simulations:
    - `git clone https://github.com/JaakkoAhola/LES-ice-02postpros`
    - `cd LES-ice-02postpros`
    - Checkout latest release tag (For acp-2019-1182 manuscript, checkout release tag: v1.1), e.g.:
        - `git checkout v1.1`
    - `python prepareICEforPlotting.py`
    - `cd ..`

# Prepare plotting library (=this repository)


- Download this repository with command:
    - `git clone https://github.com/JaakkoAhola/LES-ice-03plotting`
- `cd LES-ice-03plotting`
- checkout the latest release of the repository (For acp-2019-1182 manuscript, checkout release tag: v1.1.2), e.g.:
    - `git checkout v1.1.2`

# Making figures


- Plot figures 1-7 with command:
    - `python manuscriptFigures.py`

- Plot figure 8 with command:
    - `python updraftAnalysis.py`

- Now the files should be in `${SIMULATIONFIGUREFOLDER}`.
