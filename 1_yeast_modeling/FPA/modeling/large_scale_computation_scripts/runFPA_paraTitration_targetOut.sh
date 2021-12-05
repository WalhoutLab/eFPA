module load gurobi/911
module load matlab/R2019a
module load perl/5.10.1
module load git/2.9.5
git config --global http.sslverify "false"
matlab < FPA_Yeast_titration_targetOut.m > FPA_paraTitr_targetOut.log
