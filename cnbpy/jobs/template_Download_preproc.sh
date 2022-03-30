#!/bin/bash

#SBATCH -J ---J---
#SBATCH -N ---N---
#SBATCH --cpus-per-task=---cpus-per-task---
#SBATCH --mem-per-cpu=---mem-per-cpu---
#SBATCH --mail-user=---mail-user---
#SBATCH --error=---error---
#SBATCH --output=---output---
#SBATCH --mail-type=ALL


module load python3/anaconda/5.1.0
source activate p3env


python /home/users/yg916972/Software/cnbpy/cnbpy/jobs/Download_preproc.py ---subno--- ---local_path---