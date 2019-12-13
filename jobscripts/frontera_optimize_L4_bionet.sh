#!/bin/bash -l
export DATE=$(date +%Y%m%d_%H%M%S)
export JOB_NAME=optimize_L4_bionet_"$DATE"
sbatch <<EOT
#!/bin/bash -l
#SBATCH -J $JOB_NAME
#SBATCH -o /scratch1/06441/aaronmil/src/mouse_l4_arkhipov/logs/"$JOB_NAME".%j.o
#SBATCH -e /scratch1/06441/aaronmil/src/mouse_l4_arkhipov/logs/"$JOB_NAME".%j.e
#SBATCH -p normal
#SBATCH -N 512
#SBATCH -n 14336
#SBATCH -t 12:00:00
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -x

cd $SCRATCH/src/mouse_l4_arkhipov

ibrun -n 14336 python3 -m nested.optimize \
    --config-file-path=config/optimize_L4_bionet_config.yaml --disp --output-dir=data \
    --pop_size=64 --max_iter=2 --path_length=1 --framework=pc --procs_per_worker=224
EOT
