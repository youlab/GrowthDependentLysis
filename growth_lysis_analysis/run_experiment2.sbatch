#!/bin/bash
#SBATCH --job-name=growthlysis_exp1
#SBATCH --output=logs/growthlysis_exp1.%j.out
#SBATCH --error=logs/growthlysis_exp1.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=common

# 进入项目目录
cd /hpc/group/youlab/xc232/GrowthLysis-main

# 加载 Anaconda
module load Anaconda3/2024.02

# 激活 conda（适用于 batch 模式）
eval "$(conda shell.bash hook)"
conda activate growthlysis

echo "Python path:"
which python
python --version


# ===== Model 1 =====
#python strain_results/fit_all_strains.py --data data/first_data-lnPTR-growth_rate-plot_57strains.csv --name test_run
#python strain_results/analyze_strains.py --data data/first_data-lnPTR-growth_rate-plot_57strains.csv --name test_run

# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin50_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin50_37_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin50_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin50_37_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin50_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin50_37_1 \
    --model1-run test_run --model2-run 250511_M9CA_sensitive_Carbenicillin50_37_1 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250511_M9CA_sensitive_Carbenicillin50_37_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250511_M9CA_sensitive_Carbenicillin50_37_1



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin50_37_1.txt --name 250511_LB_sensitive_Carbenicillin50_37_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin50_37_1.txt --name 250511_LB_sensitive_Carbenicillin50_37_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin50_37_1.txt --name 250511_LB_sensitive_Carbenicillin50_37_1 \
    --model1-run test_run --model2-run 250511_LB_sensitive_Carbenicillin50_37_1 --time 5.0

python stacked_model/analyze_model4_remove_k30.py --name 250511_LB_sensitive_Carbenicillin50_37_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250511_LB_sensitive_Carbenicillin50_37_1


# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin70_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin70_37_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin70_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin70_37_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin70_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin70_37_1 \
    --model1-run test_run --model2-run 250511_M9CA_sensitive_Carbenicillin70_37_1 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250511_M9CA_sensitive_Carbenicillin70_37_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250511_M9CA_sensitive_Carbenicillin70_37_1



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin70_37_1.txt --name 250511_LB_sensitive_Carbenicillin70_37_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin70_37_1.txt --name 250511_LB_sensitive_Carbenicillin70_37_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin70_37_1.txt --name 250511_LB_sensitive_Carbenicillin70_37_1 \
    --model1-run test_run --model2-run 250511_LB_sensitive_Carbenicillin70_37_1 --time 5.0

python stacked_model/analyze_model4_remove_k30.py --name 250511_LB_sensitive_Carbenicillin70_37_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250511_LB_sensitive_Carbenicillin70_37_1


# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin90_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin90_37_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin90_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin90_37_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250511_M9CA_sensitive_Carbenicillin90_37_1.txt --name 250511_M9CA_sensitive_Carbenicillin90_37_1 \
    --model1-run test_run --model2-run 250511_M9CA_sensitive_Carbenicillin90_37_1 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250511_M9CA_sensitive_Carbenicillin90_37_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250511_M9CA_sensitive_Carbenicillin90_37_1



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin90_37_1.txt --name 250511_LB_sensitive_Carbenicillin90_37_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin90_37_1.txt --name 250511_LB_sensitive_Carbenicillin90_37_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250511_LB_sensitive_Carbenicillin90_37_1.txt --name 250511_LB_sensitive_Carbenicillin90_37_1 \
    --model1-run test_run --model2-run 250511_LB_sensitive_Carbenicillin90_37_1 --time 5.0

python stacked_model/analyze_model4_remove_k30.py --name 250511_LB_sensitive_Carbenicillin90_37_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250511_LB_sensitive_Carbenicillin90_37_1

