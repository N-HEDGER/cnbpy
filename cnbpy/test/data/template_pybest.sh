#!/bin/bash

#SBATCH -J ---J---
#SBATCH -N ---N---
#SBATCH --cpus-per-task=---cpus-per-task---
#SBATCH --mem-per-cpu=---mem-per-cpu---
#SBATCH --mail-user=---mail-user---
#SBATCH --error=---error---
#SBATCH --mail-type=ALL


pybest ---bids_dir--- --space ---space--- --save-all --noise-source ---noise-source--- --hemi L

pybest ---bids_dir--- --space ---space--- --save-all --noise-source ---noise-source--- --hemi R
