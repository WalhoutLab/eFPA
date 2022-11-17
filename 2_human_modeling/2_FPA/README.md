This guidance shows how to reproduce the enhanced flux potential analysis (eFPA) results on the yeast SIMMER dataset (Figures 1-3 and corresponding supplementary figures). 

This module contains the following folders:
- [FPA](FPA): the main folder
- [bins](bins): some auxillary functions for network analysis
- [input/GEMs](input/GEMs): metabolic network models used in this study

### Reproducing the eFPA results for yeast SIMMER dataset

To reproduce the results, please go to the [working directory](FPA/modeling/).

The entire modeling and interpretation process was splitted into 15 steps. Each step is aimed at a specific analysis and is encoded in a '.m' file. If you want to reproduce an analysis, please run the corresponding '.m' in MATLAB. 

For example:
```
matlab < step1_analysis_flux_expression_correlation.m
```
It is recommended that one runs these scripts in the MATLAB graphical interface, instead of command lines. Some of the analyses are interactive so it is best to open the scripts in MATLAB and follow the comments in the scripts.

### Computation time
Reproducing the entire modeling can be computationally challenging. The randomization test of eFPA requires a high-performance computing cluster since it calculates millions of LP optimizations. Running eFPA against the 232 flux-measured reactions with a single distance boundary parameter is usually fast, e.g., 10 minutes on a regular laptop. If a titration of distance boundary is desired, we recommend running the scripts in a lab server. If a user is interested in reproducing the randomization with access to high-performance computing resources, please follow the instructions for running the special eFPA functions for computer clusters (our script was designed for clusters managed by LSF scheduler only). 
