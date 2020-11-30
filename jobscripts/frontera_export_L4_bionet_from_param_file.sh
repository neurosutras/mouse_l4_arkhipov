#!/bin/bash -l
export DATE=$(date +%Y%m%d_%H%M%S)
export JOB_NAME=export_L4_bionet_"$DATE"
sbatch <<EOT
#!/bin/bash -l
#SBATCH -J $JOB_NAME
#SBATCH -o /scratch1/06441/aaronmil/logs/mouse_l4_arkhipov/"$JOB_NAME".%j.o
#SBATCH -e /scratch1/06441/aaronmil/logs/mouse_l4_arkhipov/"$JOB_NAME".%j.e
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 224
#SBATCH -t 2:00:00
#SBATCH --mail-user=neurosutras@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -x

cd $WORK/mouse_l4_arkhipov

ibrun -n 224 python3 -m nested.analyze \
    --config-file-path=config/optimize_L4_bionet_config.yaml --disp \
    --output-dir=$SCRATCH/data/mouse_l4_arkhipov/output --export \
    --bionet_config_file_path=frontera_config.json --framework=pc --procs_per_worker=224 \
    --param-file-path=$1 --model-key=$2
EOT
