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
DATA_PATH = pkg_resources.resource_filename('cnbpy', 'test/data')



class DATASET:
    
    """DATASET
    Class for initialising a dataset, using datalad.
    """
        
    def __init__(self,local_base,source):
        
        """
        Initialise the DATASET class.
        
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
        self.flist=os.listdir(self.local_base)
        
    def get_info_files(self,info_suffixes=['.json','.tsv']):
        self.make_flist()
        self.infolist=[x for x in self.flist if any(y in x for y in info_suffixes)]
        
        self.dset.get([i for i in self.infolist],get_data=True)
        
        
        
        
        

        
                
        
class SUBJECT:
    
    """DATASET
    Class for a subject.
    """
    
    def __init__(self,DATASET,subid,subprefix='sub-'):
        
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
        self.DATASET=DATASET
        self.subid=subid
        self.subject_path=self.subprefix+str(self.subid)
        self.local_subject_path=os.path.join(self.DATASET.local_base,self.subject_path)
        
    def get(self,flist=None):
        
        """
        Downloads the subject data using http://docs.datalad.org/en/stable/generated/man/datalad-get.html

        """
        self.DATASET.dset.get(self.subject_path,get_data=True)
        
        if flist != None:
            self.DATASET.dset.get(os.path.join(self.subject_path,[f[i] for i in flist]),get_data=True)
            
        
    def drop(self):
        
        """
        Removes the subject data using http://docs.datalad.org/en/stable/generated/man/datalad-drop.html

        """
        
        self.DATASET.dset.drop(self.subject_path)



class ABIDE_DATASET(DATASET):
    
    """DATASET
    Class for an ABIDE dataset.
    
    Additional methods for the ABIDE dataset go here.
    Will inherit the properties of a regular DATASET class.
    """
        
    def __init__(self,local_base,source):
            
        """
        Retrieves ABIDE information.
        
        Returns
        -------
        self.demofile: The demographic files (pandas dataframe).
        self.anat_q: The anatomical quality file (pandas dataframe).
        self.func_q: The functional quality file (pandas dataframe).
        """
        
        super().__init__(local_base,source)
        
        self.demofile=pd.read_csv(os.path.join(DATA_PATH,'ABIDEII_Composite_Phenotypic.csv'), encoding = "ISO-8859-1")
        self.anat_q=pd.read_csv(os.path.join(DATA_PATH,'anat_qap.csv'), encoding = "ISO-8859-1")
        self.func_q=pd.read_csv(os.path.join(DATA_PATH,'functional_qap.csv'), encoding = "ISO-8859-1")
        self.remote_base='fcp-indi/data/Projects/ABIDE2'
        
        self.fmriprep_path=os.path.join(self.remote_base,'Outputs/fmriprep/fmriprep')
        self.freesurfer_path=os.path.join(self.remote_base,'Outputs/fmriprep/freesurfer')
        self.mriqc_path=os.path.join(self.remote_base,'Outputs/fmriprep/mriqc')
        self.connection=s3fs.S3FileSystem(anon=True)
        self.awsurl='https://s3.amazonaws.com/'
        
        fmriprepped_list=self.connection.ls(self.fmriprep_path)
        self.fmriprepped_subs=[x.split('sub-', 1)[1] for x in fmriprepped_list if 'html' not in x and 'dataset_description' not in x and 'logs' not in x]
        
        
        
        
        
class ABIDE_SUBJECT(SUBJECT):
    
    """DATASET
    Class for an ABIDE dataset.
    
    Additional methods for the ABIDE dataset go here.
    Will inherit the properties of a regular DATASET class.
    """
        
    def __init__(self,DATASET,subid,outpath,subprefix='sub-'):
            
        """
        Retrieves ABIDE information.
        
        Returns
        -------
        self.demofile: The demographic files (pandas dataframe).
        self.anat_q: The anatomical quality file (pandas dataframe).
        self.func_q: The functional quality file (pandas dataframe).
        """
        super().__init__(DATASET,subid,subprefix='sub-')
    
        #self.demorow=DATASET.demofile.iloc[np.where(DATASET.demofile['SUB_ID']==int(subid))[0][0]]
        self.anatrow=DATASET.anat_q.iloc[np.where(DATASET.anat_q['Sub_ID']==int(subid))[0][0]]
        self.funcrow=DATASET.func_q.iloc[np.where(DATASET.func_q['Sub_ID']==int(subid))[0][0]]
        self.fmriprep_path=os.path.join(DATASET.fmriprep_path,self.subject_path)
        self.freesurfer_path=os.path.join(DATASET.freesurfer_path,self.subject_path)
        self.mriqc_path=os.path.join(DATASET.mriqc_path,self.subject_path)
        self.fmriprep_report=os.path.join(DATASET.awsurl,DATASET.fmriprep_path,self.subject_path+'.html')
        self.outpath=outpath
        self.local_subject_path=os.path.join(self.outpath,self.subject_path)
        
        
        
    def show_report(self,url):
        
        """
        Shows the fmriprep report for the subject. Opens up a new browser window.

        """
        
        display(Javascript('window.open("{url}");'.format(url=url)))
        
    def get_fmriprep_data(self):
        
        """
        Creates a datlad dataset for the subjects fmriprep data.
        
        Returns
        -------
        self.fmriprep_dset: The fmriprep dataset for the subject (empty until using 'get').

        """
        
        self.local_fmriprep_path=os.path.join(self.local_subject_path,'fmriprep',self.subject_path)
        
        
        if os.path.exists(self.local_fmriprep_path): # If it exists, then just connect to it.
            self.fmriprep_dset=dl.Dataset(self.local_fmriprep_path)
            os.chdir(self.fmriprep_dset.path)
        else: # Else, initialize a crawler.
            os.makedirs(self.local_fmriprep_path,exist_ok=True)
            self.fmriprep_dset=dl.Dataset.create(self.local_fmriprep_path)
            os.chdir(self.fmriprep_dset.path)
            dl.crawl_init(template='simple_s3',save=True,args=['bucket=fcp-indi','prefix=data/Projects/ABIDE2/Outputs/fmriprep/fmriprep/{subject}'.format(subject=self.subject_path),'drop_immediately=True'])
        dl.crawl()
        self.make_bids()
    
    
    def make_bids(self):
        self.bids=BIDS(os.path.join(self.local_subject_path,'fmriprep'))
    
    def get_fmriprep_outcomes(self,ses,run,task='rest'):
        
        """
        Gets the surface data for the subject given a session, run and task
        
        """
        
        wcard='ses-{ses}/func/sub-{subject}_ses-{ses}_task-{task}_run-{run}_space-fsaverage5_hemi-{hem}.func.gii'
        Lfile=wcard.format(subject=self.subid,ses=ses,task=task,run=run,hem='L')
        Rfile=wcard.format(subject=self.subid,ses=ses,task=task,run=run,hem='R')
        
        print(Lfile,Rfile)
        self.fmriprep_dset.get(Lfile)
        self.fmriprep_dset.get(Rfile)
        return os.path.join(self.local_fmriprep_path,Lfile),os.path.join(self.local_fmriprep_path,Rfile)
    
    
    def get_sess_run_combs(self,task='rest',ext='gii'):
        
        """
        Gets the available runs per session.
        
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
        
        
    def load_fmriprep_outcomes(self,ses,run,task='rest'):
        
        """
        Loads the surface data for a given session, run and task
        
        """
        
        L,R=self.get_fmriprep_outcomes(ses,run,task='rest')
        mygiftL=nb.gifti.giftiio.read(L)
        mygiftR=nb.gifti.giftiio.read(R)
        ts_dataL=mygiftL.agg_data()
        ts_dataR=mygiftR.agg_data()
        ts_data=np.vstack([ts_dataL,ts_dataR])
        return ts_data
        
    
        
        
        
        

        

        
        
        
    
    
    
    
    
