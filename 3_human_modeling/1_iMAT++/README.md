This guidance shows how to reproduce the flux prediction in Figures 2C, 2D, 3A and supp TableS5 by iMAT++. Additionally, we provided a walkthrough guidance for running iMAT++ on the generic model of [iCEL1314](http://wormflux.umassmed.edu/index.html), as well as any other metabolic network model, with any given expression dataset.

### Reproducing the integration for the <i>C. elegans</i> dual tissue model

iMAT++ is run for the seven <i>C. elegans</i> tissues for which aggregated gene expression data is available from [Cao <i>et al</i>., 2017](https://pubmed.ncbi.nlm.nih.gov/28818938/). Input files used are indicated in the walkthrough script [here](TissueOFD.m). In short, these include [the dual-tissue version of iCEL1314](./../input/Tissue.mat), [premade gene categories](./../input/geneCategories.json) and [premade epsilon values](./../input/epsilon.json). 

For converting the excel-formatted model (<i>i.e.</i>, the table downloaded from supplemental material) to COBRA-version model, see [here](loadModelFromExcelTable.m). For categorizing gene expression data as highly, lowly, rarely, and moderately expressed genes, please see [CatExp](../CatExp). The categories used in this study are provided in json format, and can be loaded into MATLAB following instructions in the walkthrough script. Users can view the categories after loading the file. Note that gene expression can be categorized using other methods as well (<i>e.g.</i>, with thresholds on percentiles, such as selecting top 10% percentile as highly expressed genes), depending on the dataset and user preferences, although we strongly recommend avoiding arbitrary expression level cutoffs for gene categorization and using CatExp instead. 

Finally, for defining the epsilon values (flux thresholds for each reaction that discriminates a real flux from zero flux during the integration), please refer to [epsilon generator](./../bins/makeEpsilonSeq.m). The method for defining epsilons can also vary from user to user. For example, a standard epsilon value of 0.001 may be used for all reactions. But we do not recommend this approach as some reactions tend to carry low flux and forcing values of 0.001 on these reactions cause excessive flux values (>1000) in other parts of the network due to mass balance rules.

To reproduce the flux distribution prediction, simply run:
```
matlab < TissueOFD.m
```
Output (flux distribution and other metrics) will be written in output/wormTissue/TissueName.mat

Please refer to the [integration function](scripts/IMATplusplus.m) for the output format.

To reproduce the Flux Variability Analysis (FVA), simply run:
```
matlab < TissueFVA.m
```
Output (upper and lower bounds) will be written in output/wormTissue/FVA/

Additionally, please be advised that the FVA analysis of tissue network is computationally intensive. One should expect the program to take couple of hours to finish in a regular 20-core lab server.

### Running iMAT++ analysis on generic <i>C. elegans</i> model with a given expression dataset

We provided a step-by-step walkthrough script for running iMAT++ on the generic [iCEL1314](http://wormflux.umassmed.edu/index.html) model with any given expression dataset. Users can start with this script to perform iMAT++ analysis on their own expression data. 

We recommend to go over this script interactively in MATLAB graphic interface. Please launch a MATLAB instance and set your current directory in this folder (1_iMAT++). Then, open and go through [the walkthrough script](walkthrough_generic.m) according to the instructions therein.

### Running FVA analysis on generic <i>C. elegans</i> model with a given expression dataset

We provided a step-by-step walkthrough script for running FVA conjoined with iMAT++ on the generic [iCEL1314](http://wormflux.umassmed.edu/index.html) model with any given expression dataset. Users can start with this script to perform iMAT++/FVA/network-building analysis on their own expression data. 

We recommend to first go through the above tutorial for IMAT++, to understand the framework of IMAT++ analysis. Then, the user can follow [the walkthrough script of large-scale FVA](walkthrough_large_scale_FVA.m) to build the context-specific model by iMAT++(and FVA).

### Running iMAT++ and FVA analysis on another metabolic model with a given expression dataset

We provided a walkthrough guidance of how to run the iMAT++ analysis on human metabolic model [Recon 2.2](https://pubmed.ncbi.nlm.nih.gov/27358602/). Publicly available consensus transcriptomes of human tissues  [RNA HPA tissue gene data](https://www.proteinatlas.org/about/download) was used as an example.

As with the tutorial for using generic <i>C. elegans</i> model, we provide a step-by-step [walkthrough script](walkthrough_generic.m) for running iMAT++ on Recon 2.2 model with any given expression dataset. This is in the second section right after the guidance for generic <i>C. elegans</i> model. Users can start with this script to develop iMAT++ analysis for their own model/data. Similarly, the FVA guidance can be found in the second section of [the walkthrough script of large-scale FVA](walkthrough_large_scale_FVA.m).

Please be advised that we provided easy speed tuning of iMAT++ and FVA, which we mentioned as the three speed levels in the paper. For Recon2.2, speed level 2 already provides acceptable computational speed, <I>i.e</I>, ~10 hours running time on a lab server to integrate 17 tissues (iMAT++ and FVA together). However, this could be even much faster (within few hours on a laptop) with speed level 3. Please test and select an appropriate speed level for your input model. For instructions on setting the speed level, please refer to the two walkthrough guidance scripts mentioned above.
