module load gurobi/810
module load matlab/R2019a
module load curl/7.60.0
module load git/2.9.5
git config --global http.sslverify "false"
matlab < HumanTissue_pipeline_FAA.m > FAA1.log

