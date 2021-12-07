# Systematic analysis of the relation between relative enzyme expression, flux, and the metabolic network
# NOTE: This repo is still under development! Not all functionalities are usable at present! You may encounter missing files as we are in the process of uploading and final reviewing the codes.
We developped the improved FPA (Flux Potential Analysis) algorithm, to systematically study the quatitative relation between enzyme expression and metabolic flux in terms of their relative levels and in the context of the metabolic network. improved FPA is applicable to any expression dataset as a tool to predict the relative metabolic flux under the principle of the local-pathway proxy role discovered in our study.

The repo includes four parts: (1) ...

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The programs were developed and tested in MATLAB R2019a/R2017a. [COnstraint-Based Reconstruction and Analysis (The COBRA Toolbox)](https://opencobra.github.io/cobratoolbox/stable/) is required to perform the analysis. Check [here](https://opencobra.github.io/cobratoolbox/stable/installation.html) for the installation guidance of COBRA Toolbox. The programs were tested for COBRA Toolbox - 2020 version, but should be compatible with an earlier version. 

The Linear Program (LP) and Mixed-Integer Linear Problem (MILP) solver used in the study was [gurobi](http://gurobi.com) 8.10. The built-in solver interface of COBRA Toolbox was used, so that we expect our program to also work with other supported solver in COBRA Toolbox. Please [see here](https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-installation) for furter information about solver availability. 

### Installing

This package doesn't require any installation or compiling. Please see the following section for reproducing the tissue expression integration result in this study, or perform the analysis on a desired metabolic network model and gene expression dataset. 

## Running the tests

The repo includes four independent modules, [iMAT++ algorithm](1_iMAT++), [original iMAT algorithm](2_iMAT), [Flux Potential Analysis (FPA)](3_FPA), and [metabolic distance calculator](MetabolicDistance). Please see the instruction within each module for running a test.

The followings are descriptions on each module (folder) listed.

[1_iMAT++](1_iMAT++): The iMAT++ module. This folder contains all the input, scripts and example codes for running the iMAT++ for <i>C. elegans</i> tissue-level integration and application on other models/datasets. 

[2_iMAT](2_iMAT): The original iMAT module. This folder contains all the input, scripts and example codes for running the original iMAT for <i>C. elegans</i> tissue-level integration. 

[3_FPA](3_FPA): The Flux Potential Analysis (FPA) module. This folder contains all the input, scripts and example codes for running the FPA for <i>C. elegans</i> tissue-level integration and application on other models/datasets. 

[MetabolicDistance](MetabolicDistance): The metabolic distance calculator module. This folder contains all the input, scripts and example codes for calculating the metabolic distance for any metabolic network model. The output distance matrix is required to run FPA. The output for the <i>C. elegans</i> models used in the reference study are available in the pertaining folders. 

[bins](bins): The shared functions for running the above mentioned analysis. These functions include modified version of some COBRA Toolbox functions and common new functions such as a molecular weight calculator.

[input](input): The shared input for running the above mentioned analysis. These inputs include <i>C. elegans</i> metabolic model and other input information.


## Contributing

Please contact us for reporting bugs or contributing purposes. Email contacting is encouraged. Please send to [Xuhang Li](mailto:xuhang.li@umassmed.edu) or [Safak Yilmaz](mailto:lutfu.yilmaz@umassmed.edu)


## Authors

* **Safak Yilmaz** - *Development of iMAT++, FPA and Metabolic Distance Calculator* - [lsafak](https://github.com/lsafak)
* **Xuhang Li** - *Matlab implementation of iMAT++/iMAT/FPA* - [XuhangLi](https://github.com/XuhangLi)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* TBD
