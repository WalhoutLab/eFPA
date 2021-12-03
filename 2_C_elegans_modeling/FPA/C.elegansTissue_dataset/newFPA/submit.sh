bsub -q short -W 120 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err1.log sh run1.sh
bsub -q short -W 120 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err2.log sh run2.sh
bsub -q short -W 120 -n 40 -R rusage[mem=2500] -R span[hosts=1] -o err3.log sh run3.sh
