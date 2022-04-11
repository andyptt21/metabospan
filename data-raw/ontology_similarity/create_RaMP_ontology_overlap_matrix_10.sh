#!/bin/bash
#SBATCH --job-name=create_RaMP_pathway_overlap_matrix
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --gres=lscratch:500

module load R
Rscript create_RaMP_ontology_overlap_matrix_10.R
