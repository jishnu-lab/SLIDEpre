#!/bin/bash
#SBATCH -t 5-00:00

#SBATCH --job-name=Step2
# partition (queue) declaration

#SBATCH --mail-user=user@pitt.edu

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --mem=300g

#SBATCH --cpus-per-task=16

module load gcc/8.2.0
module load r/4.1.0

Rscript run_ERPipeS1.R

