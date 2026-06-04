#!/bin/bash
#SBATCH --job-name=ARP53_test
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

echo "Block 2: full annotation"
Rscript scripts/02_full_annotation/02_full_annotation.R || exit 1

echo "Block 3: epithelial subset reconstruction and reclustering"
Rscript scripts/03_epi_subset/03_epi_subset_recluster.R || exit 1

echo "Block 4: Epi/PIN reannotation"
Rscript scripts/04_epi_annotation/04_reannotate_epi_PIN.R || exit 1

echo "DONE through Block 4"
