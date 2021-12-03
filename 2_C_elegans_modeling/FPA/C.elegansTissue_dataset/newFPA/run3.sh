module load gurobi/900
module load matlab/R2019a
module load curl/7.60.0
module load git/2.9.5
git config --global http.sslverify "false"
matlab <./step1_TissueFPA_wtdDist_exp_decay_base2.m> run3.log
wait
