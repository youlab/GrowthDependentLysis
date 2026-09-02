#!/usr/bin/env bash
#SBATCH --job-name=amplicon
#SBATCH --partition=youlab-gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/amplicon_%j.out
#SBATCH --error=logs/amplicon_%j.err

mkdir -p logs

echo "Job ID : $SLURM_JOB_ID"
echo "Node   : $SLURMD_NODENAME"
echo "Start  : $(date)"

source /opt/apps/rhel9/Anaconda3-2024.02/etc/profile.d/conda.sh
conda activate /hpc/group/youlab/jlei912/database/amplicon

cd /hpc/group/youlab/jlei912/zhengqing
bash amplicon_pipeline.sh all

echo "End: $(date)"
EOF