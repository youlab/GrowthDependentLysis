#!/bin/bash
#SBATCH --job-name=growthlysis_exp1
#SBATCH --output=logs/growthlysis_exp1.%j.out
#SBATCH --error=logs/growthlysis_exp1.%j.err
#SBATCH --time=30-00:00:00
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
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_30_1.txt --name 250804_M9CA_sensitive_Carbenicillin70_30_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_30_1.txt --name 250804_M9CA_sensitive_Carbenicillin70_30_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_30_1.txt --name 250804_M9CA_sensitive_Carbenicillin70_30_1 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_30_1 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_30_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_30_1



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_LB_sensitive_Carbenicillin70_30_1.txt --name 250804_LB_sensitive_Carbenicillin70_30_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_LB_sensitive_Carbenicillin70_30_1.txt --name 250804_LB_sensitive_Carbenicillin70_30_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_LB_sensitive_Carbenicillin70_30_1.txt --name 250804_LB_sensitive_Carbenicillin70_30_1 \
    --model1-run test_run --model2-run 250804_LB_sensitive_Carbenicillin70_30_1 --time 5.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_LB_sensitive_Carbenicillin70_30_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_LB_sensitive_Carbenicillin70_30_1


# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_1.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_1
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_1.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_1

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_1.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_1 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_1 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_1

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_1



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_2.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_2
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_2.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_2

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_2.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_2 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_2 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_2

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_2


# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_3.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_3
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_3.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_3

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_3.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_3 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_3 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_3

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_3



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_4.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_4
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_4.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_4

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_4.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_4 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_4 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_4

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_4


# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_5.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_5
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_5.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_5

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_5.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_5 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_5 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_5

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_5



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_6.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_6
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_6.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_6

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_6.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_6 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_6 --time 6.0 

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_6

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_6


# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_7.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_7
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_7.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_7

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_7.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_7 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_7 --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_7

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_7



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_8.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_8
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_8.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_8

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_8.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_8 \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_8 --time 6.0 

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_8

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_8



# ===== Model 2 =====
python cubic_growth_rate/fit_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_mono.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_mono
python cubic_growth_rate/analyze_all.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_mono.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_mono

# ===== Model 4 =====
python stacked_model/fit_model4_remove_k30_truncate_to_0.py --data data/od_3rep_250804_M9CA_sensitive_Carbenicillin70_37_mono.txt --name 250804_M9CA_sensitive_Carbenicillin70_37_mono \
    --model1-run test_run --model2-run 250804_M9CA_sensitive_Carbenicillin70_37_mono --time 6.0

python stacked_model/analyze_model4_remove_k30.py --name 250804_M9CA_sensitive_Carbenicillin70_37_mono

# ===== Convergence diagnostics =====
python diagnostics/check_convergence.py --name 250804_M9CA_sensitive_Carbenicillin70_37_mono


