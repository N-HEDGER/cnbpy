#!/bin/bash

#SBATCH -J ---J---
#SBATCH -N ---N---
#SBATCH --cpus-per-task=---cpus-per-task---
#SBATCH --mem-per-cpu=---mem-per-cpu---
#SBATCH --mail-user=---mail-user---
#SBATCH --error=---error---
#SBATCH --output=---output---
#SBATCH --mail-type=ALL

singularity run --cleanenv ---mounts--- ---imagepath--- ---datapath--- ---outpath--- ---analysis_level--- --no-sub --participant-label ---pid--- -w ---workpath--- --mem_gb ---mem_gb--- --omp-nthreads ---ot--- ---optionals---

