bsub -q long -W 880 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err.log sh runFPA_paraTitration_finetuned.sh
sleep 1m
bsub -q long -W 880 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err.log sh runFPA_paraTitration_targetOut.sh
sleep 1m
bsub -q long -W 880 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err.log sh runFPA_paraTitration_p1.sh
sleep 1m
bsub -q long -W 880 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err.log sh runFPA_paraTitration_p2.sh
sleep 1m
bsub -q long -W 880 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err.log sh runFPA_paraTitration_p3.sh
