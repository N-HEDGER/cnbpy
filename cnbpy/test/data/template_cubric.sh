#!/bin/bash

#SBATCH -p cubric-centos7
#SBATCH -J ---J---
#SBATCH -N ---N---
#SBATCH --cpus-per-task=---cpus-per-task---
#SBATCH --mem-per-cpu=---mem-per-cpu---
#SBATCH --mail-user=---mail-user---
#SBATCH --error=---error---
#SBATCH --output=---output---
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00

export SINGULARITYENV_TEMPLATEFLOW_HOME=---tfpath---



singularity run --cleanenv ---mounts--- ---imagepath--- ---datapath--- ---outpath--- ---analysis_level--- --participant-label ---pid--- -w ---workpath--- --output-spaces ---output_spaces--- --mem_mb ---mem_mb--- --omp-nthreads ---ot--- --nthreads ---nt--- --fs-license-file $HOME/license.txt ---optionals---

