import os

import yaml
from .utils import *
from sklearn.decomposition import PCA

import pandas as pd
import pkg_resources
import numpy as np
from .bids import BIDS
import pkg_resources
from scipy.signal import savgol_filter

job_path = pkg_resources.resource_filename('cnbpy', 'jobs')
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
        
        
    
        
class PCA_denoiser():
    
    def __init__(self,confound_frame='/content/data/confounds.tsv',ts_data='/content/data/data.gii'):
        self.confound_frame=confound_frame
        self.data=ts_data
        
        #self.csvpath=csvpath
        #self.datapath=datapath
        #self.load_csv()
        #self.load_data()

    def load_csv(self):
        self.confound_frame = pd.read_csv(self.csvpath, sep='\t', header=0,index_col=None)

    def load_data(self):
        data=nib.load(self.datapath)
        self.data=data.agg_data()

    def subset_frame(self,vars):
        self.subsetted_frame=self.confound_frame[vars]

    def prepare_frame(self):
        self.nuissance_array=np.array(self.subsetted_frame)
        medians=np.nanmedian(self.nuissance_array,axis=0)
        for c,v in enumerate(medians):
            self.nuissance_array[:,c][np.isnan(self.nuissance_array[:,c])]=medians[c]
        self.prepared_array=scipy.stats.zscore(self.nuissance_array)
    
    def PCA_regression(self,ncomps):

      # Fit PCA to the prepared array
        self.pca = PCA(n_components=ncomps)  
        self.pca_comps = self.pca.fit_transform(self.prepared_array)
 
          # Add row of 1s to the design matrix.
        self.dm = np.vstack([np.ones(self.pca_comps.shape[0]), self.pca_comps.T]).T

          # Do OLS
        self.betas = np.linalg.lstsq(self.dm, self.data.T)[0]

      # Predictions are dot product of dm and betas.
        self.yhat = np.dot(self.dm, self.betas).T

      # Get model rsq
        self.rsq = 1-(self.data-self.yhat).var(-1)/self.data.var(-1)

      # Get residuals
        self.resid=self.data-self.yhat

        self.resid+= np.nanmean(self.data,axis=-1)[:,np.newaxis]
        self.denoised_data=self.resid
        
        
def sg_filter(data, polyorder=3, deriv=0, window_length = 120,tr=1):
    """ Applies a savitsky-golay filter to a nifti-file.
    Fits a savitsky-golay filter to a 4D fMRI nifti-file and subtracts the
    fitted data from the original data to effectively remove low-frequency
    signals.
    Parameters
    ----------
    in_file : str
    Absolute path to nifti-file.
    polyorder : int (default: 3)
    Order of polynomials to use in filter.
    deriv : int (default: 0)
    Number of derivatives to use in filter.
    window_length : int (default: 120)
    Window length in seconds.
    Returns
    -------
    out_file : str
        Absolute path to filtered nifti-file.
    """


    # TR must be in seconds
    if tr < 0.01:
        tr = np.round(tr * 1000, decimals=3)
    if tr > 20:
        tr = tr / 1000.0

    window = np.int(window_length / tr)

    # Window must be odd
    if window % 2 == 0:
        window += 1

    data_filt = savgol_filter(data, window_length=window, polyorder=polyorder,
                              deriv=deriv, axis=1, mode='nearest')

    data_filt = data - data_filt + data_filt.mean(axis=-1)[:, np.newaxis]

    return data_filt