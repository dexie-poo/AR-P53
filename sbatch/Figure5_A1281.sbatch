#!/bin/bash
#SBATCH --job-name=Fig5_A1281
#SBATCH --partition=cpu_long
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=/gpfs/scratch/leungd02/ARP53/logs/%x_%j.out
#SBATCH --error=/gpfs/scratch/leungd02/ARP53/logs/%x_%j.err

module purge
module load gcc/15.2.0
module load r/4.4.2

source ~/.bashrc

cd /gpfs/scratch/leungd02/ARP53

mkdir -p logs

Rscript scripts/Figure5_A1281/run_pipeline.R
