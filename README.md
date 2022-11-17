# Enzyme levels broadly influence flux throughout the metabolic network

## NOTICE 11/17/2022
This repo is under constant development. Please let us know if you encounter any problem (email or raising an issue are both welcome). In addition, we are developing a user-friendly version of eFPA. If you are interested, don't forget to Star our repo! We will work harder with the positive feedback!

## Intro
We developped the enhanced Flux Potential Analysis (eFPA) algorithm, to systematically study the quantitative relation between enzyme expression and metabolic flux in terms of their relative levels and in the context of the metabolic network. eFPA is applicable to any expression dataset as a tool to predict the relative metabolic flux under the principle of enzyme reach discovered in our study.

The repo includes four parts: (1) data and codes to reproduce the yeast analysis (2) data and codes to reproduce the human tissue analysis (3) metabolic distance calculatior for making the distance matrix from a given metabolic network model and (4) a tutorial to run eFPA on your own expression data and metabolic network model.

Please note, if you only want to use the eFPA to analyze your own data, you may directly go to [eFPA tutorial](4_eFPA_tutorial) that contains a tutorial of eFPA analysis. To reproduce our analyses presented in the paper, please refer to the sections 1 and 2. 

Please be advised that if you are going to run FPA on a model other than the models used in our study, you need to calculate the distance matrix for your model, which can be done by following instructions in [metabolic distance calculator](3_distance_calculation). This computation may be intensive on large models (i.e., over 2000 rxns x 2000 mets). You may need to use high-performance computation resources or reach out to us for assistance. 

We provide pre-calculated distance matrix for commonly used models (human, C. elegans and yeast, more to come!). please visit [our website](http://wormflux.umassmed.edu/download.php) for downloading.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The programs were developed and tested in MATLAB R2019a/R2020b. [COnstraint-Based Reconstruction and Analysis (The COBRA Toolbox)](https://opencobra.github.io/cobratoolbox/stable/) is required to perform the analysis. Check [here](https://opencobra.github.io/cobratoolbox/stable/installation.html) for the installation guidance of COBRA Toolbox. The programs were tested for COBRA Toolbox - 2020 version, but should be compatible with an earlier or later version. 

The Linear Program (LP) and Mixed-Integer Linear Problem (MILP) solver used in the study was [gurobi](http://gurobi.com) 9.0.0. The built-in solver interface of COBRA Toolbox was used, so that we expect our program to also work with other supported solver in COBRA Toolbox. Please [see here](https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-installation) for furter information about solver availability. 

### Large files
Several input files (e.g., the distance matrix of human model) are large in size (>100 MB), so they are not included in this repository. If you are interested in reproduce the results in our paper, please download the zipped full repository from [our website](http://wormflux.umassmed.edu/download.php) for the missing large files. The large files are only required to reproduce our yeast and human analysis, rather than applying eFPA to your custom dataset.

### Installing

This package doesn't require any installation or compiling. Please see the following section for reproducing the analyses this study, or use the eFPA as a tool for analyzing your own data.


## Running the tests

The repo includes four independent parts, [yeast modeling](1_yeast_modeling), [human tissue modeling](2_human_modeling), [metabolic distance calculator](3_distance_calculation), and [eFPA tutorial](4_eFPA_tutorial). Please see the instruction within each module for running a test.


The followings are descriptions on each module (folder) listed.

[yeast modeling](1_yeast_modeling): This folder contains all the input, scripts, example pipelines and results to reproduce our analyses related to the modeling in yeast.

[human tissue modeling](2_human_modeling): This folder contains all the input, scripts, example pipelines and results to reproduce our analyses related to the modeling in human.

[metabolic distance calculator](3_distance_calculation): This folder contains all the input, scripts and example codes for calculating the metabolic distance for any metabolic network model. The output distance matrix is required to run eFPA.

[eFPA tutorial](4_eFPA_tutorial): This folder contains a tutorial to perform eFPA on a given model and expression dataset.


## Contributing

Please contact us for reporting bugs or contributing purposes. Email contacting is encouraged. Please send to [Xuhang Li](mailto:xuhang.li@umassmed.edu) or [Safak Yilmaz](mailto:lutfu.yilmaz@umassmed.edu)


## Authors

* **Xuhang Li** [XuhangLi](https://github.com/XuhangLi)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* We thank Walhout lab members for their support in this project.
