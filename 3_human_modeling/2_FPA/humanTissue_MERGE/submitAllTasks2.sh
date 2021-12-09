bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro4.log -o cLogs/pro4.log -J pro4 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue transporter\" > logs/pro4.log && wait"
sleep 30m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna4.log -o cLogs/rna4.log -J rna4 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_all 6 tissue transporter\" > logs/rna4.log && wait"
sleep 30m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna5.log -o cLogs/rna5.log -J rna5 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_all 6 tissue demand\" > logs/rna5.log && wait"
sleep 30m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna6.log -o cLogs/rna6.log -J rna6 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue allMetDemand\" > logs/rna6.log && wait"
sleep 30m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna7.log -o cLogs/rna7.log -J rna7 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_all 6 tissue allMetDemand\" > logs/rna7.log && wait"
sleep 30m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna8.log -o cLogs/rna8.log -J rna8 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_all 6 tissue rxn\" > logs/rna8.log && wait"
sleep 30m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro7.log -o cLogs/pro7.log -J pro7 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 100 naive transporter\" > logs/pro7.log && wait"
