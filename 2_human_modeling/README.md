This guidance shows how to reproduce the enhanced flux potential analysis (eFPA) results on the human tissues (Figures 4 and corresponding supplementary figures). 

This module contains the following folders:
- [1_iMAT++](1_iMAT++): the codes to build the tissue-specific networks
- [2_FPA](2_FPA): the codes to perform the tissue eFPA analysis
- [bins](bins): auxillary functions for network analysis
- [input/GEMs](input/GEMs): [DISCONTINUED - TO BE REMOVED]

### Reproducing the eFPA results for human tissue analysis

The human tissue eFPA composes of two sections: building the tissue networks and running the eFPA on the tissue networks. Building context-specific networks is not required for running eFPA, however, our previous study (Yilmaz et al, Mol. Sys. Bio., 2020) on tissue metabolism in C. elegans showed that running FPA on tissue networks better reveals tissue metabolism. Therefore, we built the tissue networks for humans following the previous methods. 

To reproduce the tissue networks, please go to the [working directory](1_iMAT++/).

To reproduce the eFPA analysis, please go to the [working directory](2_FPA/).
