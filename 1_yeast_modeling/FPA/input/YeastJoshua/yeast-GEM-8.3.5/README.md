# yeast-GEM: The consensus genome-scale metabolic model of _Saccharomyces cerevisiae_

[![DOI](https://zenodo.org/badge/52777598.svg)](https://zenodo.org/badge/latestdoi/52777598) [![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem.svg)](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem) [![Join the chat at https://gitter.im/SysBioChalmers/yeast-GEM](https://badges.gitter.im/SysBioChalmers/yeast-GEM.svg)](https://gitter.im/SysBioChalmers/yeast-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

* Brief Model Description:

This repository contains the current consensus genome-scale metabolic model of _Saccharomyces cerevisiae_. It is the continuation of the legacy project [yeastnet](https://sourceforge.net/projects/yeast/). For the latest release please [click here](https://github.com/SysBioChalmers/yeast-GEM/releases).

* Model KeyWords:

**GEM Category:** species; **Utilisation:** experimental data reconstruction, multi-omics integrative analysis, _in silico_ strain design, model template; **Field:** metabolic-network reconstruction; **Type of Model:** reconstruction, curated; **Model Source:** YeastMetabolicNetwork; **Omic Source:** genomics, metabolomics; **Taxonomy:** _Saccharomyces cerevisiae_; **Metabolic System:** general metabolism; **Bioreactor**; **Strain:** S288C; **Condition:** aerobic, glucose-limited, defined media;

* Last update: 2020-03-31

* Main Model Descriptors:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
|:-------:|:--------------:|:---------:|:----------:|:-----:|
|_Saccharomyces cerevisiae_|[Yeast 7.6](https://sourceforge.net/projects/yeast/)|3991|2691|1147|

This repository is administered by Benjamín J. Sánchez ([@BenjaSanchez](https://github.com/benjasanchez)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Installation

### Required Software - User:

* Matlab user:
  * A functional Matlab installation (MATLAB 7.3 or higher).
  * The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
* Python user:
  * Python 2.7, 3.4, 3.5 or 3.6
  * [cobrapy](https://github.com/opencobra/cobrapy)

### Required Software - Contributor:

* Both of the previous Matlab requirements.
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).
* A [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Dependencies - Recommended Software:
* For Matlab, the [libSBML MATLAB API](https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface/) (version 5.17.0 is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.

### Installation Instructions
* For users: Clone it from [`master`](https://github.com/SysBioChalmers/yeast-GEM) in the Github repo, or just download [the latest release](https://github.com/SysBioChalmers/yeast-GEM/releases).
* For contributors: Fork it to your Github account, and create a new branch from [`devel`](https://github.com/SysBioChalmers/yeast-GEM/tree/devel).

## Usage

Make sure to load/save the model with the corresponding wrapper functions!
* In Matlab:
  * Loading: `complementaryScripts/loadYeastModel.m`
  * Saving: `complementaryScripts/saveYeastModel.m`
* In Python:
  * Loading: `complementaryScripts/loadYeastModel.py`
  * Saving: currently unavailable

## Model Files

The model is available in `.xml`, `.txt`, `.yml`, `.mat` and `.xlsx` (the last 2 extensions only in `master`). Additionally, the following 2 files are available:
* `dependencies.txt`: Tracks versions of toolboxes & SBML used for saving the model.
* `boundaryMets.txt`: Contains a list of all boundary metabolites in model, listing the id and name.

## Citation

* If you use yeast-GEM please cite the yeast8 paper:
  > Lu, H. et al. _A consensus S. cerevisiae metabolic model Yeast8 and its ecosystem for comprehensively probing cellular metabolism._ Nature Communications 10, 3586 (2019). https://doi.org/10.1038/s41467-019-11581-3.

* Additionally, all yeast-GEM releases are archived in [Zenodo](https://zenodo.org/badge/latestdoi/52777598), for you to cite the specific version of yeast-GEM that you used in your study, to ensure reproducibility. You should always cite the original publication + the specific version, for instance:
  > _The yeast consensus genome-scale model [Lu et al. 2019], version 8.3.4 [Sánchez et al. 2019], was used._

  Find the citation details for your specific version [here](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%221494182%22&sort=-publication_date&all_versions=True).

## Contributing

Contributions are always welcome! Please read the [contributions guideline](https://github.com/SysBioChalmers/yeast-GEM/blob/master/.github/CONTRIBUTING.md) to get started.

## Contributors

* [Mihail Anton](https://www.chalmers.se/en/staff/Pages/mihail-anton.aspx) ([@mihai-sysbio](https://github.com/mihai-sysbio)), Chalmers University of Technology, Sweden
* [Moritz Beber](https://www.dtu.dk/english/service/phonebook/person?id=121829&tab=2&qt=dtupublicationquery) ([@Midnighter](https://github.com/Midnighter)), Danish Technical University, Denmark
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden
* [Dimitra Lappa](https://www.chalmers.se/en/Staff/Pages/lappa.aspx) ([@demilappa](https://github.com/demilappa)), Chalmers University of Technology, Sweden
* [Feiran Li](https://www.chalmers.se/en/staff/Pages/feiranl.aspx) ([@feiranl](https://github.com/feiranl)), Chalmers University of Technology, Sweden
* [Christian Lieven](https://www.dtu.dk/english/service/phonebook/person?id=103199&tab=2&qt=dtupublicationquery) ([@ChristianLieven](https://github.com/ChristianLieven)), Danish Technical University, Denmark
* [Hongzhong Lu](https://www.chalmers.se/en/Staff/Pages/luho.aspx) ([@hongzhonglu](https://github.com/hongzhonglu)), Chalmers University of Technology, Sweden
* [Simonas Marcišauskas](https://www.chalmers.se/en/Staff/Pages/simmarc.aspx) ([@simas232](https://github.com/simas232)), Chalmers University of Technology, Sweden
* [Thomas Pfau](https://wwwen.uni.lu/research/fstc/life_sciences_research_unit/research_areas/systems_biology/people/thomas_pfau) ([@tpfau](https://github.com/tpfau)), University of Luxembourg, Luxembourg
* [Benjamín J. Sánchez](https://www.chalmers.se/en/staff/Pages/bensan.aspx) ([@BenjaSanchez](https://github.com/benjasanchez)), Chalmers University of Technology, Sweden
