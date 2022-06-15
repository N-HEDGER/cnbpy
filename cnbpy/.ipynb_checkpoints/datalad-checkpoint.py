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
import re

data_path = pkg_resources.resource_filename('cnbpy', 'test/data')


import yaml
from .utils import *


class Dataset:
    
    """Dataset
    Class for initialising a dataset, using datalad.
    """
        
    def __init__(self,local_base,source):
        
        """
        Initialise the Dataset class.
        
        Parameters
        ----------
        local_base: Path to the local directory where the dataset is to be stored (string)
        source: path to the datalad BIDS repository, as listed on http://datasets.datalad.org/ (string)
        isabide: If this is the ABIDE dataset, will load in extra info.
        
        Returns
        -------
        self.source: source
        self.local_base: local base
        self.dset: datalad class
        
        Notes
        -------
        1. The idea is to set the 'source' to a BIDS directory, that is populated by 'sub-' directories.
        2. This doesn't seem to behave well on my removable storage.  
        
        """
        
        self.source=source
        self.local_base=local_base
        self.dset=dl.Dataset(self.local_base)
        
        
    def install(self):
        
        """
        Installs the dataset using http://docs.datalad.org/en/stable/generated/man/datalad-install.html
        It will create a bunch of symbolic links to the dataset.
        The subject data itself can be retrieved using the 'subject' class.

        """
        
        dl.install(source=self.source,path=self.local_base)
        
        
    def get_all(self):
        
        """
        Gets (downloads) all the files in the dataset.
        """
        
        self.dset.get()
        
    
    def remove(self):
        
        """
        Removes the dataset using http://docs.datalad.org/en/stable/generated/man/datalad-remove.html
        """
        
        self.dset.remove()
        
        
    def get_size(self):
        
        """
        Reports the size of the dataset, in terms of that stored locally and that which is annexed.
        """
        
        a=self.dset.status(annex='all')
        
    def make_flist(self):
        
        """
        Lists all of the files in the dataset.
        """
        
        self.flist=os.listdir(self.local_base)
        
    def get_info_files(self,info_suffixes=['.json','.tsv']):
        self.make_flist()
        self.infolist=[x for x in self.flist if any(y in x for y in info_suffixes)]
        
        self.dset.get([i for i in self.infolist],get_data=True)
        
        
        
                
class Subject:
    
    """DATASET
    Class for a subject.
    """
    
    def __init__(self,Dataset,subid,subprefix='sub-'):
        
        """
        Initialise the SUBJECT class.
        
        Parameters
        ----------
        DATASET: Path to the local directory where the dataset is to be stored (string).
        subid: Whatever identifies the subject after the 'sub-' prefix (string or numeric).
        subprefix: Whatever the subprefix is (the default should be ok for BIDS) (string).
        
        Returns
        -------
        self.subprefix: subprefix.
        self.DATASET: DATASET.
        self.subject_path: String that identifies subject directory within the dataset.

        """
        self.subprefix=subprefix
        self.Dataset=Dataset
        self.subid=subid
        self.subject_path=self.subprefix+str(self.subid)
        self.local_subject_path=os.path.join(self.Dataset.local_base,self.subject_path)
        
    def get(self,flist=None):
        
        """
        Downloads the subject data using http://docs.datalad.org/en/stable/generated/man/datalad-get.html

        """
        self.Dataset.dset.get(self.subject_path,get_data=True)
        
        if flist != None:
            self.Dataset.dset.get(os.path.join(self.subject_path,[f[i] for i in flist]),get_data=True)
            
        
    def drop(self):
        
        """
        Removes the subject data using http://docs.datalad.org/en/stable/generated/man/datalad-drop.html

        """
        
        self.Dataset.dset.drop(self.subject_path)


data_path = pkg_resources.resource_filename('cnbpy', 'test/data')
ABIDE_yaml=os.path.join(data_path,'ABIDE_config.yml')        


class Abide_Dataset(object):
    
    """Abide_Dataset
    Class for an ABIDE dataset.

    """
        
    def __init__(self,local_base,yaml=ABIDE_yaml,data_path=data_path,dicts=['rpaths','csvs']):
        
        """ init
        
        Parameters
        ----------
        local_base: Path to the local directory where the dataset is to be stored (string)
        yaml: yaml file that defines the remote paths to the dataset and information about the csvs.
        data_path: path to the anatomical and functional quality files and phenotypic info.
        
        Returns
        -------        

        """
        self.local_base=local_base
        self.yaml=ABIDE_yaml
        self.data_path=data_path
        self.dicts=dicts
        
        for sdict in self.dicts:
            self.internalize_yaml(sdict)
        
        self.get_abide_info()
        self.make_connection()
        
    
    def list_fmriprepped_subjects(self):
        
        """ list_fmriprepped_subjects 
        Lists the subjects that have been fmriprepped
        Returns
        -------
        self.fmriprepped_list: Full list of files in remote fmriprep directory.
        self.fmriprepped_subs: List of subjects with an fmriprep directory.
        self.giftisubs: List of subjects with gifti files in their fmriprep output.

        """
        
        
        self.fmriprepped_list=self.connection.ls(self.fmriprep_path)
        self.fmriprepped_subs=[x.split('sub-', 1)[1] for x in self.fmriprepped_list if 'html' not in x and 'dataset_description' not in x and 'logs' not in x]
        
        Lgiftifiles=self.connection.glob(os.path.join(self.fmriprep_path,self.gifti_glob))

        self.giftisubs=[[int(s) for s in re.findall(r'\b\d+\b', gfile)][0] for gfile in Lgiftifiles]
        
    
    def make_connection(self):
        
        """ make_connection
        Uses s3fs to make a connection to the S3 services.
        
        Returns
        -------
        self.connection: s3fs.core.S3FileSystem
        """
        
        
        self.connection=s3fs.S3FileSystem(anon=True)
    
    def get_abide_info(self):
        
        """get_abide_info
        Reads in the CSV files associated with the ABIDE dataset.
        
        Returns
        -------
        self.demofile: Demographic info on subjects.
        self.anat_q:  Anatomical quality for subjects.
        self.func_q: Functional quality for subjects.
        """
        
        
        self.demofile=pd.read_csv(os.path.join(self.data_path,self.demofile_name), encoding = self.encoding)
        self.anat_q=pd.read_csv(os.path.join(self.data_path,self.anatfile_name), encoding = self.encoding)
        self.func_q=pd.read_csv(os.path.join(self.data_path,self.funcfile_name), encoding = self.encoding)
        
    
    def internalize_yaml(self,subdict):

        """ internalize_config_yaml
        internalises a subdictionary of a yaml file - makes the elements a property of the class.

        """
        with open(self.yaml, 'r') as f:
            self.y = yaml.safe_load(f)

        cdict = self.y[subdict]

        for key in cdict.keys():
            setattr(self, key, cdict[key])
            
                
        
class Abide_Subject(Subject):
    
    """Abide_Subject
    Class for an ABIDE subject.

    """
        
    def __init__(self,Dataset,subid,subprefix='sub-',dicts=['data','crawler','singularity','pybest']):
            
        """ init
        
        Parameters
        ----------
        Dataset: An AbideDataset object
        subid: The subjects id number (string)
        subprefix: The prefix before the subject ID (string).

        """
        

        super().__init__(Dataset,subid,subprefix)
        
        
        
        self.Dataset=Dataset
        self.yaml=self.Dataset.yaml
        self.dicts=dicts
        
        for sdict in self.dicts:
            self.internalize_yaml(sdict)
                
        self.get_subject_info()
        self.setup_paths()
        
        
    def setup_paths(self):
        
        """ setup_paths
        internalises a subdictionary of a yaml file - makes the elements a property of the class.
        
        Returns
        -------
        self.local_fmriprep_path_base: The base directory to store fmriprepped subjects.
        self.local_fmriprep_path: The 
        """
        
        self.local_fmriprep_path_base=os.path.join(self.local_subject_path,'fmriprep')
        self.local_fmriprep_path=os.path.join(self.local_fmriprep_path_base,self.subject_path)
        
    
    def get_subject_info(self):
        
        """ get_subject_info
        Gets information about the subject based on dataset csv files.
        
        Returns
        -------
        self.anatrow: Row of anatomical quality file.
        self.funcrow: Row of functional quality file.
        self.fmriprep_path: fmriprep path
        self.fmriprep_report: fmriprep report location.
        """
        
        
        self.anatrow=self.Dataset.anat_q.iloc[np.where(self.Dataset.anat_q['Sub_ID']==int(self.subid))[0][0]]
        self.funcrow=self.Dataset.func_q.iloc[np.where(self.Dataset.func_q['Sub_ID']==int(self.subid))[0][0]]
        self.fmriprep_path=os.path.join(self.Dataset.fmriprep_path,self.subject_path)
        self.fmriprep_report=os.path.join(self.Dataset.aws_url,self.Dataset.fmriprep_path,self.subject_path+'.html')
        
        
    def internalize_yaml(self,subdict):

        """internalize_config_yaml

        """
        
        with open(self.yaml, 'r') as f:
            self.y = yaml.safe_load(f)

        cdict = self.y[subdict]
        
        setattr(self, subdict, cdict)

        for key in cdict.keys():
            setattr(self, key, cdict[key])
            
            
    def show_report(self,url):
        
        """ show_report
        Shows the fmriprep report for the subject. Opens up a new browser window.

        """
        
        display(Javascript('window.open("{url}");'.format(url=url)))
        
    def make_fmriprep_dataset(self):
        
        """ make_fmriprep_dataset
        Creates a datalad dataset for the subjects fmriprep data.
        
        Returns
        -------
        self.fmriprep_dset: The fmriprep dataset for the subject (empty until using 'get').

        """
        
        dstring='prefix='+self.fmriprep_path[9:]
        if os.path.exists(self.local_fmriprep_path): # If it exists, then just connect to it.
            self.fmriprep_dset=dl.Dataset(self.local_fmriprep_path)
        else: # Else, initialize a crawler.
            os.makedirs(self.local_fmriprep_path,exist_ok=True)
            self.fmriprep_dset=dl.Dataset.create(self.local_fmriprep_path)
            os.chdir(self.fmriprep_dset.path)
            dl.crawl_init(template=self.template,save=True,args=[self.bucketstring,dstring,'drop_immediately=True',self.fmriprep_exclude])

         
        
    def crawl_fmriprep_dataset(self):
        
        """ crawl_fmriprep_dataset
        Creates a local symbolic link to the dataset on the system.
        
        Returns
        -------
        self.fmriprep_dset: The fmriprep dataset for the subject (empty until using 'get').

        """
        
        os.chdir(self.fmriprep_dset.path)
        dl.crawl()
        self.make_bids()
        
    
    
    def make_bids(self):
        """ make_bids
        Makes a BIDS class for the subject fmriprep directory.
        
        Returns
        -------
        self.bids: bids class.

        """
        
        self.bids=BIDS(os.path.join(self.local_subject_path,'fmriprep'))
        
        
    def get_functional_data(self):
        
        """ get_functional_data
        Gets the functional data for each session.
        
        Returns
        -------
        self.funcdirs: Functional directories for each subject.

        """
        
        
        self.funcdirs=['ses-{ses}/func'.format(ses=ses) for ses in self.bids.sessions[0]]
        
        self.fmriprep_dset.get([funcdir for funcdir in self.funcdirs])
        
        
    
    def get_fmriprep_outcomes(self,ses,run,task='rest'):
        
        """ get_fmriprep_outcomes
        Gets the surface data for the subject given a session, run and task
        
        Parameters
        ----------
        ses: session number (int)
        run: run number (int)
        
        Returns
        -------
        paths for the data for the right and left hemisphere for a given session and run.
        
        
        """
        
        Lfile=self.func_wildcard.format(subject=self.subid,ses=ses,task=task,run=run,hem='L')
        Rfile=self.func_wildcard.format(subject=self.subid,ses=ses,task=task,run=run,hem='R')
        
        self.fmriprep_dset.get(Lfile)
        self.fmriprep_dset.get(Rfile)
        return os.path.join(self.local_fmriprep_path,Lfile),os.path.join(self.local_fmriprep_path,Rfile)
    
    
    def load_fmriprep_outcomes(self,ses,run,task='rest'):
        
        """ load_fmriprep_outcomes
        Loads the surface data for a given session, run and task
        
        Parameters
        ----------
        ses: session number (int)
        run: run number (int)
        
        Returns
        -------
        The surface data (numpy array (vertices, time))
        
        
        """
                
        L,R=self.get_fmriprep_outcomes(ses,run,task='rest')
        mygiftL=nb.gifti.giftiio.read(L)
        mygiftR=nb.gifti.giftiio.read(R)
        ts_dataL=mygiftL.agg_data()
        ts_dataR=mygiftR.agg_data()
        ts_data=np.vstack([ts_dataL,ts_dataR])
        return ts_data
    
    def get_pybest_outcomes(self,ses,run,task='rest'):
        
        """ get_pybest_outcomes
        Gets the pybest-denoised surface data for the subject given a session, run and task
        
        Parameters
        ----------
        ses: session number (int)
        run: run number (int)
        
        Returns
        -------
        paths for the denoised data for the right and left hemisphere for a given session and run.
        
        
        """
        
        
        Lfile=self.denoised_wildcard.format(subject=self.subject_path,ses=ses,task='rest',run=run,hem='L')
        Rfile=self.denoised_wildcard.format(subject=self.subject_path,ses=ses,task='rest',run=run,hem='R')
        
        return os.path.join(self.denoised_out_dir,Lfile),os.path.join(self.denoised_out_dir,Rfile)
    
    def load_pybest_outcomes(self,ses,run,task='rest'):
        
        """ load_pybest_outcomes
        Loads the pybest-denoised surface data for the subject given a session, run and task
        
        Parameters
        ----------
        ses: session number (int)
        run: run number (int)
        
        Returns
        -------
        The denoised surface data (numpy array (vertices, time))
        
        
        """
        
        
        L,R=self.get_pybest_outcomes(ses,run,task='rest')
        ts_data=np.hstack([np.load(L),np.load(R)])
        return ts_data.T
    
    def get_sess_run_combs(self,task='rest',ext='gii'):
        
        """ get_sess_run_combs
        Gets the available runs per session.
        
          Parameters
        ----------
        task: session number (int)
        ext: file extension to search for
        
        Returns
        -------
        self.files_per_session: Number of 
        self.sess_run_combs: A list of session and run combinations (session, run number)
        
        
        """
        
        # Get the valid gifti files for the given task per session.
        self.files_per_session=[self.bids.find_funcs(sub=self.subid,task=task,ses=int(ses),ext=ext) for ses in self.bids.sessions[0]]
        self.runs_per_session=[int(len(i)/2) for i in self.files_per_session] # Get the number of runs per session (number of files /2).
        
        self.sess_run_combs=[]
        # For each session
        for i,v in enumerate(self.bids.sessions[0]):
            # For the range of runs
            for run in range(self.runs_per_session[i]):
                # Append the current session to each of the runs within the session
                self.sess_run_combs.append([v,str(run+1)])
        
    def drop_fmriprep_dataset(self):
        """ drop_fmriprep_dataset
        Drops the subjects fmriprep dataset (can be re-downloaded with get_functional_data)
        """

        
        self.fmriprep_dset.drop()
        
    
    

    
        
    
        
        
        
        

        

        
        
        
    
    
    
    
    
