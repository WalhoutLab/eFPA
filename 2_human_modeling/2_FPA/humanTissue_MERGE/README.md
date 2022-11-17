
The entire modeling and interpretation process was splitted into 6 steps. Each step is aimed at a specific analysis and is encoded in one or multiple '.m' file. If you want to reproduce an analysis, please run the corresponding '.m' in MATLAB. 

It is recommended that one runs these scripts in the MATLAB graphical interface, instead of command lines. Many of the analyses are interactive so it is best to open the scripts in MATLAB and follow the comments in the scripts.

### Computation time
Reproducing the entire modeling can be computationally challenging due to the large size of human metabolic network model. We performed our computation in a computing cluster with 512 cores, which completes the eFPA within ~30 minutes. We don't recommend developing your analysis on human model from the scripts here, as they were customized to run eFPA on our cluster (we have custom cluster monitoring functions built in). Please refer to the tutorial section if you are interested in running human model eFPA (but we recommend running on at least a multi-core server!).
