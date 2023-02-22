"""
This script will prepare your Slurm script so that you can batch all your fMRIprep jobs together ready to run on the CUBRIC cluster computer.

This script was created by Carolyn McNabb as part of the cnbpy project.
"""
#Import os
import os

#change into the directory where this repository is stored
os.chdir('/home/sapcm15/gitlab/cnbpy')

#import the FMRIPREP submodule
from cnbpy.fmriprep_cubric import FMRIPREP

#create an object called myfmriprep that loads information about the fmriprep job
myfmriprep=FMRIPREP(jobname='WAND_fmriprep')

#populate fields described in the template.sh file using information from the dictionaries in fmriprep_config.yml
myfmriprep.populate_all()

#write the code to a file that gets stored in your job folder
myfmriprep.writeout()

print("Check that you created a file in your job folder called e.g., myjobfile.sh and if it exists, run it using sbatch myjobfile.sh")