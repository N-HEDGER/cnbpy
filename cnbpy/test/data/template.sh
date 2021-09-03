#!/bin/bash

#SBATCH -J ---J---
#SBATCH -N ---N---
#SBATCH --cpus-per-task=---cpus-per-task---
#SBATCH --mem-per-cpu=---mem-per-cpu---
#SBATCH --mail-user=---mail-user---
#SBATCH --error=---error---

#SBATCH --mail-type=ALL




singularity run --cleanenv ---mounts--- ---imagepath--- ---datapath--- ---outpath--- --participant-label ---pid--- -w ---workpath--- --output-spaces ---output_spaces--- --mem_mb ---mem_mb--- --fs-subjects-dir ---fsdir--- --omp-nthreads ---ot--- --nthreads ---nt--- --fs-license-file ---fsli--- ---optionals---

