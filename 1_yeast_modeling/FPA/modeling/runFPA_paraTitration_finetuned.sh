module load gurobi/900
module load matlab/R2019a
module load perl/5.10.1
module load git/2.9.5
git config --global http.sslverify "false"
matlab < FPA_Yeast_titration_fineTuned.m > FPA_paraTitr_finetuned.log
