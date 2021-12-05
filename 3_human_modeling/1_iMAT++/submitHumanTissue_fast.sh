bsub -q large -W 24:00 -n 20 -R rusage[mem=4000] -e err_iHumanPip_p1f.log -J iHuman_iMAT_p1f "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"HumanTissue_pipeline_fast 1 8\" > iHumanPip_p1f.log && wait"
sleep 30
bsub -q large -W 48:00 -n 20 -R rusage[mem=4000] -e err_iHumanPip_p2f.log -J iHuman_iMAT_p2f "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"HumanTissue_pipeline_fast 9 16\" > iHumanPip_p2f.log && wait"
sleep 30
bsub -q large -W 24:00 -n 20 -R rusage[mem=4000] -e err_iHumanPip_p3f.log -J iHuman_iMAT_p3f "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"HumanTissue_pipeline_fast 17 24\" > iHumanPip_p3f.log && wait"
sleep 30
bsub -q large -W 24:00 -n 20 -R rusage[mem=4000] -e err_iHumanPip_p4f.log -J iHuman_iMAT_p4f "module load gurobi/900 && module load matlab/R2019a && module load git/2.9.5 && git config --global http.sslverify "false" && matlab -nodisplay -nosplash -r \"HumanTissue_pipeline_fast 25 32\" > iHumanPip_p4f.log && wait"
