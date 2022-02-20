#!/bin/bash

################################################################################
# List scripts
################################################################################

SCRIPTDIR="/path/to/src/"

SCRIPT1="01_bigWigSummary.R"
# SCRIPT2="04a_discreteBayesianNetworks_mergedRep.R"
# SCRIPT3="plotMethAcceHeatmap.R"

################################################################################
# Summarize bigWig files over enhancer and promotors
################################################################################

# Initialize conda
source ~/.bashrc

# Activate a certain environment
conda activate bayesianNetEnv

# Submit jobs
for i in {1..6}; do
for j in {1..2}; do

#SBATCH --job-name=bigWigSummary $i $j
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -p gpu-legacy
#SBATCH -N 1
#SBATCH --gres=gpu:4
#SBATCH --output=bigWigSumm_output $i $j.out
#SBATCH --mail-user=QianWu.Liao@stud.uni-heidelberg.de
#SBATCH --mail-type=ALL

Rscript ${SCRIPTDIR}${SCRIPT1} $i $j

done
done

################################################################################
# Bayesian network of validation data
################################################################################

# Initialize conda
# source ~/.bashrc

# Activate a certain environment
# conda activate bayesianNetEnv

# Submit jobs

#SBATCH --job-name=BayesianNetworks
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -p gpu-legacy
#SBATCH -N 1
#SBATCH --gres=gpu:4
#SBATCH --output=BN_output.out
#SBATCH --mail-user=QianWu.Liao@stud.uni-heidelberg.de
#SBATCH --mail-type=ALL

# Rscript ${SCRIPTDIR}${SCRIPT2}

################################################################################
# Visualize methylation and accessibility around TSSs
################################################################################

# Initialize conda
# source ~/.bashrc

# Activate a certain environment
# conda activate bayesianNetEnv

# Submit jobs

#SBATCH --job-name=plotHeatmap
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -p gpu-legacy
#SBATCH -N 1
#SBATCH --gres=gpu:4
#SBATCH --output=plotHeatmap.out
#SBATCH --mail-user=QianWu.Liao@stud.uni-heidelberg.de
#SBATCH --mail-type=ALL

# Rscript ${SCRIPTDIR}${SCRIPT3}

## Fin
