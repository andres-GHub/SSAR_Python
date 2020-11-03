# Self-similar aftershock rates (SSAR) model code
Python code used to carry out the statistical analysis of magnitude correlations for the self-similar aftershock rates (SSAR) model presented in https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.022314, the catalog ("TimeMagLT.txt") used in the analysis, along with the data for the plots.

SSAR_mag_corr.py is the main program which carries out the statistical analysis described in the paper:

https://npg.copernicus.org/articles/27/1/2020/npg-27-1-2020.html
Magnitude correlations in a self-similar aftershock rates model of seismicity
Authors: Andres F. Zambrano Moreno and Jörn Davidsen

for a seismic catalog ("TimeMagLT.txt", which needs to be in the same folder as SSAR_mag_corr.py). The main output is the magnitude correlations file "dP_XXX.txt" which contains three colums:

m0 | dP | error (symmetric)

and shown as m0 vs dP plots. All other files, i.e. FUNC_XXXX.py, are functions used in the main program.

###########

This folder also contains the SSAR model catalog "TimeMagLT.txt" used in the analysis.

