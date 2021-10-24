# Exponential Random Graph Models
This repo contains most of the code used for Chapter 4 of my PhD thesis. In a nutshell, it turns out that due to the infinite dimensional nature of the considered ERGMs they are precisely described by the mean field theory which is supported by the results of Monte Carlo simulations implemented in this repo. The code is structured as follows.

## Code
Every considered homogeneous ERGM apart from the MF model has a dedicated directory, namely ./q_star_model, ./triad_model and ./AKS_model, each of which contains */figures directory containing the figures from the thesis. Every image from the */figures directory has its corresponding (with the same name) python script in the "scripts" directory, which was used to generate the image. Most of these scripts rely on the data stored in */data directory in the form of .csv files with descriptive names. The data files needed for most of the scripts are already in the correct locations within */data directory, so running the scripts with Python 3 should (assuming the installation of all required packages) reproduce the images as they appear in the */figures directory. To run the actual simulations of the model one needs to compile the C++ code from the */main_sim, which relies on the files from the ./GraphOS directory. This can be done (assuming g++ is installed) by running the */main_sim/compile.sh script which places the binaries into */bin directory from which they should be launched in order for the results to be saved in their correct locations. Before that it is recommended to comment out the compilation of the unneeded code from */main_sim/compile.sh to reduce the time it takes to run the compilation script (should not take too long anyways). Finally, the ./general_images directory contains some figures that don't rely on any data and follows the same structure, i.e. .png files from */figures correspond to Python scripts in */scripts which produce them.
