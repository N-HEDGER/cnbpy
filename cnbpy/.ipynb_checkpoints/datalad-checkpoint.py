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



class Abide_Dataset(Dataset):
    
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
        
        self.data_path = pkg_resources.resource_filename('cnbpy', 'test/data')
        self.yaml=os.path.join(DATA_PATH,'ABIDE_config.yml')
        self.dicts=['rpaths','csvs']
        
        for sdict in self.dicts:
            self.internalize_yaml(sdict)
        
        self.get_abide_info()
        self.make_connection()
        self.list_fmriprepped_subjects()
        
    
    def list_fmriprepped_subjects(self):
        self.fmriprepped_list=self.connection.ls(self.fmriprep_path)
        self.fmriprepped_subs=[x.split('sub-', 1)[1] for x in self.fmriprepped_list if 'html' not in x and 'dataset_description' not in x and 'logs' not in x]
        
        #self.freesurfered_list=self.connection.ls(self.freesurfer_path)
        #self.freesurfered_subs=[x.split('sub-', 1)[1] for x in self.freesurfered_list if 'html' not in x and 'dataset_description' not in x and 'logs' not in x]
    
    def make_connection(self):
        """make_connectionto the S3 bucket
        """
        
        
        self.connection=s3fs.S3FileSystem(anon=True)
    
    def get_abide_info(self):
        """ read the ABIDE information.
        """
        
        self.demofile=pd.read_csv(os.path.join(self.data_path,self.demofile_name), encoding = self.encoding)
        self.anat_q=pd.read_csv(os.path.join(self.data_path,self.anatfile_name), encoding = self.encoding)
        self.func_q=pd.read_csv(os.path.join(self.data_path,self.funcfile_name), encoding = self.encoding)
        
    
    def internalize_yaml(self,subdict):

        """internalize_config_yaml

        """
        with open(self.yaml, 'r') as f:
            self.y = yaml.safe_load(f)

        cdict = self.y[subdict]

        for key in cdict.keys():
            setattr(self, key, cdict[key])
            
                
        
class Abide_Subject(Subject):
    
    """Abide_Subject
    Class for an ABIDE subject.
    
    Additional methods for the ABIDE Subject go here.
    Will inherit the properties of a regular Subject class.
    """
        
    def __init__(self,Dataset,subid,subprefix='sub-'):
            
        """
        
        Returns
        -------

        """
        super().__init__(Dataset,subid,subprefix='sub-')
        
        self.yaml=self.Dataset.yaml
        
        self.dicts=['data','crawler','neuropythy','singularity','pybest']
        
        for sdict in self.dicts:
            self.internalize_yaml(sdict)
                
        self.get_subject_info()
        self.setup_paths()
        
        
    def setup_paths(self):
        self.local_fmriprep_path_base=os.path.join(self.local_subject_path,'fmriprep')
        
        self.local_freesurfer_path=os.path.join(self.freesurfer_subject_dir,self.subject_path)
        self.local_fmriprep_path=os.path.join(self.local_fmriprep_path_base,self.subject_path)
        #self.local_freesurfer_path=os.path.join(self.local_freesurfer_path_base,self.subject_path)
        self.pybestdir=os.path.join(self.local_subject_path,'pybest')
        
    
    def get_subject_info(self):
        self.anatrow=self.Dataset.anat_q.iloc[np.where(self.Dataset.anat_q['Sub_ID']==int(self.subid))[0][0]]
        self.funcrow=self.Dataset.func_q.iloc[np.where(self.Dataset.func_q['Sub_ID']==int(self.subid))[0][0]]
        self.fmriprep_path=os.path.join(self.Dataset.fmriprep_path,self.subject_path)
        self.freesurfer_path=os.path.join(self.Dataset.freesurfer_path,self.subject_path)
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
        
        """
        Shows the fmriprep report for the subject. Opens up a new browser window.

        """
        
        display(Javascript('window.open("{url}");'.format(url=url)))
        
    def make_fmriprep_dataset(self):
        
        """
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

         
    def make_freesurfer_dataset(self):
        
        """
        Creates a datlad dataset for the subjects freesurfer data.
        
        Returns
        -------
        self.freesurfer_dset: The freesurfer dataset for the subject (empty until using 'get').

        """
        dstring='prefix='+self.freesurfer_path[9:]
        
        if os.path.exists(self.local_freesurfer_path): # If it exists, then just connect to it.
            self.freesurfer_dset=dl.Dataset(self.local_freesurfer_path)
        else: # Else, initialize a crawler.
            os.makedirs(self.local_freesurfer_path,exist_ok=True)
            self.freesurfer_dset=dl.Dataset.create(self.local_freesurfer_path)
            os.chdir(self.freesurfer_dset.path)
            dl.crawl_init(template=self.template,save=True,args=[self.bucketstring,dstring,'drop_immediately=True',self.freesurfer_exclude])
            
                        
    def crawl_freesurfer_dataset(self):
        os.chdir(self.freesurfer_dset.path)
        dl.crawl()
        
    def crawl_fmriprep_dataset(self):
        os.chdir(self.fmriprep_dset.path)
        dl.crawl()
        self.make_bids()
        
        
    def get_freesurfer_outcomes(self):
        
        """
        Gets the required freesurfer folders
        
        Returns
        -------
        self.freesurfer_dset: The fmriprep dataset for the subject (empty until using 'get').

        """
        
        
        for folder in self.required_freesurfer_subfolders:
            self.freesurfer_dset.get(folder)
    
    
    def make_bids(self):
        self.bids=BIDS(os.path.join(self.local_subject_path,'fmriprep'))
        
        
    def get_functional_data(self):
        
        self.funcdirs=['ses-{ses}/func'.format(ses=ses) for ses in self.bids.sessions[0]]
        
        self.fmriprep_dset.get([funcdir for funcdir in self.funcdirs])
    
    def get_fmriprep_outcomes(self,ses,run,task='rest'):
        
        """
        Gets the surface data for the subject given a session, run and task
        
        """
        
        Lfile=self.func_wildcard.format(subject=self.subid,ses=ses,task=task,run=run,hem='L')
        Rfile=self.func_wildcard.format(subject=self.subid,ses=ses,task=task,run=run,hem='R')
        
        self.fmriprep_dset.get(Lfile)
        self.fmriprep_dset.get(Rfile)
        return os.path.join(self.local_fmriprep_path,Lfile),os.path.join(self.local_fmriprep_path,Rfile)
    
    def get_pybest_outcomes(self,ses,run,task='rest'):
        
        Lfile=self.denoised_wildcard.format(subject=self.subject_path,ses=ses,task='rest',run=run,hem='L')
        Rfile=self.denoised_wildcard.format(subject=self.subject_path,ses=ses,task='rest',run=run,hem='R')
        
        return os.path.join(self.pybestdir,Lfile),os.path.join(self.pybestdir,Rfile)
    
    def load_pybest_outcomes(self,ses,run,task='rest'):
        L,R=self.get_pybest_outcomes(ses,run,task='rest')
        ts_data=np.hstack([np.load(L),np.load(R)])
        return ts_data.T
    
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
    
    
    
    def drop_all_manual(self,dset):
        things2drop=listdir_nohidden(dset.path)
        dset.drop([i for i in things2drop])
        
    
    
    def make_benson_outputs(self):
        
        """
        Produces the list of files produced by the benson retinotopy command.
        
        Returns
        -------
        self.benson_outputs: Filenames for the files produced for the benson retinotopy command.

        """
        
        
        self.benson_outputs=[]
        self.resampled_benson_outputs=[]
        for a in ['lh','rh']:
            for b in self.produced_params:
                self.benson_outputs.append(self.outsurf_wildcard.format(hem=a,param=b,fmt=self.output_fmt))
               
    
    
    def load_benson_output(self,hem,param):
         
        surf_file=os.path.join(self.local_freesurfer_path,'surf',self.resampled_outsurf_wildcard.format(hem=hem,param=param,fmt=self.output_fmt))
        surfobj=nib.load(surf_file)
        surfdat=surfobj.get_data()
        return surfdat.flatten()[:self.resamp_verts]
    
    def load_all_benson_outputs(self):
        
        benson_outputs=[]
        for b in self.produced_params:
            benson_outputs.append(np.concatenate([self.load_benson_output('lh',b),self.load_benson_output('rh',b)]))
        benson_frame=pd.DataFrame(np.array(benson_outputs).T)
        benson_frame.columns=self.produced_params
        benson_frame['hem']=np.repeat([1,2],self.resamp_verts)
        benson_frame['angle'][benson_frame['hem']==2]=-benson_frame['angle'][benson_frame['hem']==2]
        benson_frame['mask']=benson_frame['varea']==1

        self.benson_frame=benson_frame
        self.V1maskL,self.V1maskR=np.array(benson_frame[benson_frame['hem']==1]['mask']),np.array(benson_frame[benson_frame['hem']==2]['mask'])
    
    
    def do_benson_retinotopy(self):
        
        self.make_benson_outputs()
        
        """
        Performs the benson retinotopy command on the subjects freesurfer directory.
        
        """
        
        self.benson_cmd=self.benson_wildcard.format(sub=self.subject_path,subdir=self.freesurfer_subject_dir,fmt=self.output_fmt)
        print(self.benson_cmd)
        os.system(self.benson_cmd)
        
        
    def run_singularity(self,scriptname): 
    
        self.singularity_cmd=self.singularity_cmd_wildcard.format(lmount=self.local_freesurfer_path_base,smount=self.smount,simage=self.singularity_img,scriptname=scriptname)
        
        
    def make_scripts(self):
        
        """
        Makes a set of scripts to resample the benson outputs to fsaverage 5 space.
        
        """
        
        
        self.script_dicts = [{'---sub---':self.subject_path,'---src---':os.path.join('/fs',self.subject_path,'surf',item),'---trg---': 'fsaverage',
               '---trgsurf---':os.path.join('/fs',self.subject_path,'surf',item.replace('.mgz','_resamp.mgz')),'---hemi---':item[:2]} for item in self.benson_outputs]
        self.scriptfiles=[]
        supdict={'---fsli---':'/fs/license.txt','---sdir---':'/fs'}
        for i,v in enumerate(self.script_dicts):
            mp=Script_Populator(self.script_yml,self.script_template,self.script_outdir,self.script_prefix+'_'+str(i))
            #mp.yaml=v
            mp.yaml={**v, **supdict}
            mp.populate()
            mp.writeout()            
            self.scriptfiles.append(os.path.join('/fs','scripts',os.path.split(mp.outfile)[-1]))
            
    
    def make_singularity_commands(self):
        
        """
        Makes a set of singularity commands.
        
        """
        
        self.singularity_commands=[]
        for script,num in enumerate(self.scriptfiles):
            self.singularity_commands.append( self.singularity_cmd_wildcard.format(lfspath=self.freesurfer_subject_dir,sfspath=self.smount,script_prefix=self.script_prefix,simage=self.singularity_img,scriptloc=num))
    
    
    def execute_singularity_command(self,cmd):
        os.system(cmd)
        
    def execute_singularity_commands(self):
        
        """
        Executes the singularity commands.
        
        """
        
        for i,cmd in enumerate(self.singularity_commands):
            self.execute_singularity_command(cmd)
    
    def resample_benson_retinotopy(self):
        self.make_scripts()
        self.make_singularity_commands()
        self.execute_singularity_commands()
        
    
        
    
        
        
        
        

        

        
        
        
    
    
    
    
    
