This guidance shows how to reproduce the tissue network building of the 32 human tissues by iMAT++ algorithm.


### Running the scripts

The entire iMAT++ modeling was splitted into 4 steps, from building the human serum model that better handles nutrient input compared with original human 1 model, to the iMAT++ network building. 

Compared with the original iMAT++ algorithm, two new functionalities were added, including (1) allowing user to specific a 'core' set of reactions that would be turned on regardless of their associated gene expression and (2) a mathmetically equivalent but computationally more efficient (~2-fold increase in speed) algorithm to build the networks based on FVA. The algorithm is named Flux Accessibility Analysis (FAA) and built in the provided pipeline.

You may want to try the new iMAT++ functionalities if you are using our iMAT++ algorithm previously published!

### Computation time
We recommend running the scripts in a lab server. The average speed for iMAT++ pipeline is ~30min per tissue in a 20-core regular lab server.
