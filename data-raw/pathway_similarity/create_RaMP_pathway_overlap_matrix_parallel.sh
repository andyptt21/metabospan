#!/bin/bash
#SBATCH --job-name=create_pathway_overlap_matrix_parallel
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --gres=lscratch:500
#SBATCH --mem=24g
#SBATCH --cpus-per-task='16'

module load R
Rscript create_RaMP_pathway_overlap_matrix_parallel.R
