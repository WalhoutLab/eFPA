# Predicting flux from enzyme expression and metabolic network architecture
# NOTE: This repo is still under development! Not all functionalities are usable at present! You may encounter missing files as we are in the process of uploading and cleaning up the codes.
We developped the enhanced Flux Potential Analysis (eFPA) algorithm, to systematically study the quatitative relation between enzyme expression and metabolic flux in terms of their relative levels and in the context of the metabolic network. eFPA is applicable to any expression dataset as a tool to predict the relative metabolic flux under the principle of the local-pathway enzyme reach discovered in our study.

The repo includes four parts: (1) data and codes to reproduce the yeast flux analysis (2) data and codes to reproduce the human tissue analysis (3) metabolic distance calculatior for making the distance matrix of a metabolic network model and (4) a tutorial to run eFPA on your own expression data and metabolic network model.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

The programs were developed and tested in MATLAB R2019a/R2020b. [COnstraint-Based Reconstruction and Analysis (The COBRA Toolbox)](https://opencobra.github.io/cobratoolbox/stable/) is required to perform the analysis. Check [here](https://opencobra.github.io/cobratoolbox/stable/installation.html) for the installation guidance of COBRA Toolbox. The programs were tested for COBRA Toolbox - 2020 version, but should be compatible with an earlier version. 

The Linear Program (LP) and Mixed-Integer Linear Problem (MILP) solver used in the study was [gurobi](http://gurobi.com) 8.10. The built-in solver interface of COBRA Toolbox was used, so that we expect our program to also work with other supported solver in COBRA Toolbox. Please [see here](https://opencobra.github.io/cobratoolbox/stable/installation.html#solver-installation) for furter information about solver availability. 

### Installing

This package doesn't require any installation or compiling. Please see the following section for reproducing the analyses this study, or use the improved FPA as a tool for analyzing your own data.


## Running the tests

The repo includes five independent parts, [yeast modeling](1_yeast_modeling), [C. elegans tissue metabolism modeling](2_C_elegans_modeling), [human tissue modeling](3_human_modeling), [metabolic distance calculator](4_distance_calculation), and [improved FPA tutorial](5_improvedFPA_tutorial). Please see the instruction within each module for running a test.

Please note, if you only want to use the improved FPA tool for analyzing your own data, you may directly go to [improved FPA tutorial](5_improvedFPA_tutorial) that contains a full tutorial of FPA analysis. To reproduce our analyses presented in the paper, please refer to the sections #1-3. Please be advised that if you are going to run FPA on a model other than the models used in our study, you need to calculate the distance matrix for your model, which can be achieved by instructions in [metabolic distance calculator](4_distance_calculation). This computation may be intensive on large models (i.e., > 2000 rxns x 2000 mets). You may need to use high-performance computation resources or reach out to us for assistance. 

The followings are descriptions on each module (folder) listed.

[yeast modeling](1_yeast_modeling): This folder contains all the input, scripts and example pipelines and results to reproduce our analyses related to the modeling in yeast.

[C. elegans tissue metabolism modeling](2_C_elegans_modeling): This folder contains all the input, scripts and example pipelines and results to reproduce our analyses related to the modeling in C. elegans.

[human tissue modeling](3_human_modeling): This folder contains all the input, scripts and example pipelines and results to reproduce our analyses related to the modeling in human.

[metabolic distance calculator](4_distance_calculation): This folder contains all the input, scripts and example codes for calculating the metabolic distance for any metabolic network model. The output distance matrix is required to run FPA.

[improved FPA tutorial](5_improvedFPA_tutorial): This folder contains a full tutorial to perform FPA on a given model and expression dataset.


## Contributing

Please contact us for reporting bugs or contributing purposes. Email contacting is encouraged. Please send to [Xuhang Li](mailto:xuhang.li@umassmed.edu) or [Safak Yilmaz](mailto:lutfu.yilmaz@umassmed.edu)


## Authors

* **Xuhang Li** [XuhangLi](https://github.com/XuhangLi)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* TBD
