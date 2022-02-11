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
data_path = pkg_resources.resource_filename('cnbpy', 'test/data')

import yaml
from .utils import *


pybest_template=os.path.join(data_path,'template_pybest.sh')
pybest_yaml=os.path.join(data_path,'pybest_config.yml')

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
        self.jobname=jobstr+'_'+self.subject.subid
        supdict={'---J---':self.jobname,'---bids_dir---':self.subject.local_fmriprep_path_base,'---error---':os.path.join(self.subject.local_subject_path,self.jobname+'.err')}
        self.pybest_job=Script_Populator(yaml_file=yaml,template_file=template,out_dir=self.subject.local_subject_path,jobname=self.jobname,suffix='.sh',supdict=supdict)
        self.pybest_job.populate()
        self.pybest_job.writeout()
        
    def run_pybest(self):
        self.make_pybest_script()
        self.pybest_job.execute()
        
        
    
        
        