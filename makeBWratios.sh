#! /bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="bwRatios"
#SBATCH --time=0-8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=all
#SBATCH --mem-per-cpu=16G

Rscript makeBWratios.R
