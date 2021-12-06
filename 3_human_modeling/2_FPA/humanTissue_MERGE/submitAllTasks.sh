# for benchmarking the use of TS weighting
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro1.log -o cLogs/pro1.log -J pro1 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_all 100 naive rxn\" > logs/pro1.log && wait"
sleep 1.5h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro2.log -o cLogs/pro2.log -J pro2 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_raw_all 100 naive rxn\" > logs/pro2.log && wait"
sleep 1.5h
# for tissue-metabolism landscape by the RNA and protein data
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro3.log -o cLogs/pro3.log -J pro3 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue rxn\" > logs/pro3.log && wait"
sleep 1.5h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna1.log -o cLogs/rna1.log -J rna1 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue rxn\" > logs/rna1.log && wait"
sleep 1.5h

# HMDB metabolite FPA, protein
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro4.log -o cLogs/pro4.log -J pro4 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue transporter\" > logs/pro4.log && wait"
sleep 1h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro5.log -o cLogs/pro5.log -J pro5 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue demand\" > logs/pro5.log && wait"
sleep 1h
# HMDB metabolite FPA RNA
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna2.log -o cLogs/rna2.log -J rna2 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue transporter\" > logs/rna2.log && wait"
sleep 1h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna3.log -o cLogs/rna3.log -J rna3 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue demand\" > logs/rna3.log && wait"
sleep 1h
# HMDB metabolite FPA RNA - all genes
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna4.log -o cLogs/rna4.log -J rna4 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_all 6 tissue transporter\" > logs/rna4.log && wait"
sleep 1h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna5.log -o cLogs/rna5.log -J rna5 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_all 6 tissue demand\" > logs/rna5.log && wait"
sleep 1h

# all metabolite FPA - all metabolite demand
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro6.log -o cLogs/pro6.log -J pro6 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue allMetDemand\" > logs/pro6.log && wait"
sleep 1m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna6.log -o cLogs/rna6.log -J rna6 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue allMetDemand\" > logs/rna6.log && wait"
sleep 1m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna7.log -o cLogs/rna7.log -J rna7 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_all 6 tissue allMetDemand\" > logs/rna7.log && wait"
sleep 1m

# metabolite FPA - no network integration (through very local original FPA, dorder = 100)
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro7.log -o cLogs/pro7.log -J pro7 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 100 naive transporter\" > logs/pro7.log && wait"
sleep 1m
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro8.log -o cLogs/pro8.log -J pro8 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 100 naive rxn\" > logs/pro8.log && wait"
