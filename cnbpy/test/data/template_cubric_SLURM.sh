#!/bin/bash

#SBATCH -p cubric-centos7
#SBATCH -J ---J---
#SBATCH -N ---N---
#SBATCH --cpus-per-task=---cpus-per-task---
#SBATCH --mem-per-cpu=---mem-per-cpu---
#SBATCH --mail-user=---mail-user---
#SBATCH --error=/cubric/collab/314_wand/bids/derivatives/fmriprep/jobs/WAND_fMRIprep_%j.err
#SBATCH --output=/cubric/collab/314_wand/bids/derivatives/fmriprep/jobs/WAND_fMRIprep_%j.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --array=1-183

export SINGULARITYENV_TEMPLATEFLOW_HOME=---tfpath---

# use the SLURM_ARRAY_TASK_ID to define the subject ID (sub) from the subject_list.txt file
sub=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/sapcm15/gitlab/cnbpy/cnbpy/test/data/subject_list_cubric.txt)

echo "Running fMRIprep for ${sub}"
singularity run --cleanenv ---mounts--- ---imagepath--- ---datapath--- ---outpath--- ---analysis_level--- --participant-label ${sub} -w ---workpath--- --output-spaces ---output_spaces--- --mem_mb ---mem_mb--- --omp-nthreads ---ot--- --nthreads ---nt--- --fs-license-file $HOME/license.txt ---optionals---

