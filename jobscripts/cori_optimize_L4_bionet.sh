#!/bin/bash -l
export DATE=$(date +%Y%m%d_%H%M%S)
export JOB_NAME=optimize_L4_bionet_"$DATE"
sbatch <<EOT
#!/bin/bash -l
#SBATCH -J $JOB_NAME
#SBATCH -o /global/cscratch1/sd/aaronmil/optimize_L4_bionet/logs/"$JOB_NAME".%j.o
#SBATCH -e /global/cscratch1/sd/aaronmil/optimize_L4_bionet/logs/"$JOB_NAME".%j.e
#SBATCH -q premium
#SBATCH -N 800
#SBATCH -L SCRATCH
#SBATCH -C haswell
#SBATCH -t 03:00:00
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -x

cd $HOME/mouse_l4_arkhipov

srun -N 800 -n 25600 -c 2 --cpu-bind=cores python -m nested.optimize \
    --config-file-path=config/cori_optimize_L4_bionet_config.yaml --disp --output-dir=$SCRATCH/optimize_L4_bionet \
    --pop_size=200 --max_iter=2 --path_length=1 --framework=pc --procs_per_worker=128
EOT
