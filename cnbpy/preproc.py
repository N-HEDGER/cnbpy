import os

import nest_asyncio
nest_asyncio.apply()

import nibabel as nb
import datalad.api as dl
import pandas as pd
import pkg_resources
import numpy as np
import s3fs
from IPython.display import Javascript
import datalad.api as dl
from .bids import BIDS
import pkg_resources
from copy import deepcopy

job_path = pkg_resources.resource_filename('cnbpy', 'jobs')

import yaml
from .utils import *


pybest_template=os.path.join(job_path,'template_pybest.sh')
pybest_yaml=os.path.join(job_path,'pybest_config.yml')

class Preprocessor:
    
    """Dataset
    Class for initialising a dataset, using datalad.
    """
        
    def __init__(self,subject):
        self.subject=subject
        self.yaml=self.subject.yaml
        self.internalize_yaml('pybest')
        
    def internalize_yaml(self,subdict):

        """internalize_config_yaml

        """
        with open(self.yaml, 'r') as f:
            self.y = yaml.safe_load(f)

        cdict = self.y[subdict]

        for key in cdict.keys():
            setattr(self, key, cdict[key])
            
            
            
    def make_pybest_script(self,yaml=pybest_yaml,template=pybest_template,jobstr='pybest'):
        
        """ make_pybest_script
        Makes a pybest script for a given subject.
        
          Parameters
        ----------
        yaml: yaml file that contains information about the fmriprep job
        template: A template script file to be populated with the information in pybest_yaml
        """
        
        
        self.jobname=jobstr+'_'+self.subject.subid
        
        
        # We create an additional dictionary to populate with subject-relevant information.
        supdict={'---J---':self.jobname,'---bids_dir---':self.subject.local_fmriprep_path_base,'---error---':os.path.join(self.denoised_out_dir,self.jobname+'.err'),'---output---':os.path.join(self.denoised_out_dir,self.jobname+'.out'),'---out_dir---':self.denoised_out_dir}
        
        
        self.pybest_job=Script_Populator(yaml_file=yaml,template_file=template,out_dir=self.denoised_out_dir,jobname=self.jobname,suffix='.sh',supdict=supdict)
        self.pybest_job.populate()
        self.pybest_job.writeout()
        
    def run_pybest(self):
        """ make_pybest_script
        Runs the pybest job.
        """
        
        self.make_pybest_script()
        self.pybest_job.execute()
        
        
    
        
        