# These figures are plotted with Ubuntu 18.04 Linux distribution.
# Thus, these insctructions are for Linux operating system.
# Lines on this file that don't start with # are command line command that can/should be used

#####################
### PREREQUISITES ###
#####################

# Python run within Anaconda environment (downloaded from https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh)
# Conda environment file in this repository: environment.yml
# Change prefix line in environment file to be in line with your conda environment.
# Conda environment can be installed after installation of Anaconda with command:
conda env create -f environment.yml
# Activate conda environment.
conda activate data

# Set environment variable ${SIMULATIONDATAROOTFOLDER} with folder that you have stored all the individual simulation folder.
export SIMULATIONDATAROOTFOLDER=/folder/that/contains/all/the/simulations

# Set environment variable ${SIMULATIONFIGUREFOLDER} with the folder that you want to store your figures.
export SIMULATIONFIGUREFOLDER=/folder/where/you/want/to/store/the/figures

# Go to directory of your liking where you can download needed repositories.

# Download plotting library.
git clone https://github.com/JaakkoAhola/LES-03plotting
# Use latest release tag. (For acp-2019-1182 manuscript, use release tag: v1.0).

############################
### END OF PREREQUISITES ###
############################

# Prepare  .csv data file about simulations.
git clone https://github.com/JaakkoAhola/LES-ice-02postpros
# Use latest release tag (For acp-2019-1182 manuscript, use release tag: v1.0).
cd LES-ice-02postpros
python prepareICEforPlotting.py
cd ..


# Download this repository with command.
git clone https://github.com/JaakkoAhola/LES-ice-03plotting
# Use the latest release of repository # Use latest release tag (For acp-2019-1182 manuscript, use release tag: v1.0).

cd LES-ice-03plotting
# Plot figures with command.
python manuscriptFigures.py

# Now the files should be in $SIMULATIONFIGUREFOLDER.
