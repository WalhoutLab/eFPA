# for benchmarking the use of TS weighting
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro1.log -o cLogs/pro1.log -J pro1 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_all 100 naive rxn\" > logs/pro1.log && wait"
sleep 1.5h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro2.log -o cLogs/pro2.log -J pro2 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_raw_all 100 naive rxn\" > logs/pro2.log && wait"
sleep 1.5h
# for benchmarking the tissue network
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro3.log -o cLogs/pro3.log -J pro3 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_all 6 naive rxn\" > logs/pro3.log && wait"
sleep 1.5h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro4.log -o cLogs/pro4.log -J pro4 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_all 6 tissue rxn\" > logs/pro4.log && wait"
sleep 1.5h
# for benchmarking the RNA and protein (biological)
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro5.log -o cLogs/pro5.log -J pro5 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue rxn\" > logs/pro5.log && wait"
sleep 1.5h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna1.log -o cLogs/rna1.log -J rna1 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue rxn\" > logs/rna1.log && wait"
sleep 1.5h
# for benchmarking the RNA and protein (methodological), will use TS in stead of raw when full RNA TS is available
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro6.log -o cLogs/pro6.log -J pro6 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_raw_all 6 tissue rxn\" > logs/pro6.log && wait"
sleep 1.5h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna2.log -o cLogs/rna2.log -J rna2 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_raw_all 6 tissue rxn\" > logs/rna2.log && wait"
sleep 1.5h
# for benchmarking the original FPA, new FPA and new FPA with original distance
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro7.log -o cLogs/pro7.log -J pro7 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_all 1.5 tissue rxn\" > logs/pro7.log && wait"
sleep 1.5h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro8.log -o cLogs/pro8.log -J pro8 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original new protein_TS_all 6 tissue rxn\" > logs/pro8.log && wait"
sleep 1.5h
# metabolite FPA, protein
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro9.log -o cLogs/pro9.log -J pro9 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue transporter\" > logs/pro9.log && wait"
sleep 1h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro10.log -o cLogs/pro10.log -J pro10 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue demand\" > logs/pro10.log && wait"
sleep 1h
# metabolite FPA RNA
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna3.log -o cLogs/rna3.log -J rna3 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue transporter\" > logs/rna3.log && wait"
sleep 1h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna4.log -o cLogs/rna4.log -J rna4 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue demand\" > logs/rna4.log && wait"
sleep 1h
# metabolite FPA RNA - all
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna5.log -o cLogs/rna5.log -J rna5 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_raw_all 6 tissue transporter\" > logs/rna5.log && wait"
sleep 1h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna6.log -o cLogs/rna6.log -J rna6 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_raw_all 6 tissue demand\" > logs/rna6.log && wait"
sleep 1h
# metabolite FPA protein - all

# metabolite FPA - all metabolite demand 
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro11.log -o cLogs/pro11.log -J pro11 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new protein_TS_common 6 tissue allMetDemand\" > logs/pro11.log && wait"
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna7.log -o cLogs/rna7.log -J rna7 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_TS_common 6 tissue allMetDemand\" > logs/rna7.log && wait"
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_rna8.log -o cLogs/rna8.log -J rna8 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic weighted new RNA_raw_all 6 tissue allMetDemand\" > logs/rna8.log && wait"

# metabolite FPA - no network integration
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro12.log -o cLogs/pro12.log -J pro12 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 100 naive transporter\" > logs/pro12.log && wait"

# metabolite original MERGE - protein 
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro13.log -o cLogs/pro13.log -J pro13 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 1.5 tissue transporter\" > logs/pro13.log && wait"
sleep 1h
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro14.log -o cLogs/pro14.log -J pro14 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 1.5 tissue demand\" > logs/pro14.log && wait"
sleep 1
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro15.log -o cLogs/pro15.log -J pro15 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 1.5 tissue allMetDemand\" > logs/pro15.log && wait"
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro16.log -o cLogs/pro16.log -J pro16 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 1.5 tissue rxn\" > logs/pro16.log && wait"
# some others 
bsub -q long -W 12:00 -n 1 -R rusage[mem=30000] -e errs/err_pro17.log -o cLogs/pro17.log -J pro17 "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"TissueFPA_generic original original protein_TS_common 100 naive rxn\" > logs/pro17.log && wait"
