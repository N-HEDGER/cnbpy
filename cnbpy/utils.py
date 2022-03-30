import pkg_resources
import yaml,os
import re
import nibabel as nib
import neuropythy as npy
import numpy as np

from scipy import stats
from sklearn.model_selection import KFold
import random
import os


DATA_PATH = pkg_resources.resource_filename('cnbpy', 'test/data')

MMP_PATH = pkg_resources.resource_filename('cnbpy', 'test/data/HCP-MMP')


    



class Script_Populator:
    
    """Script_Populator
    Class for populating a script.
    
    The idea is to take a template file and populate it with information contained in a yaml file.
    
    
    """
    
    
    def __init__(self,yaml_file,template_file,out_dir,jobname='myfmriprep',suffix='.sh',**kwargs):
        
        """
        Initialise the FMRIPREP class.
        
        Parameters
        ----------
        yaml_file: yaml file containing the information to be fed into the jobscript.
        template_file: template file for the jobscript 
        out_dir: Where to save the populated script.
        jobname: The name given to the job.
        
        
        An additional 'supdict' dictionary can be provided in kwargs to populate additional information.
        This is useful in the case where the script needs to be populated on the fly.
        
        Parameters
        ----------
        
        self.outfile: Script output location.
        self.working_string: The unpopulated template script.
        
        """
        self.jobname=jobname # Name to be given to job.
        
        self.yaml_file=yaml_file
        
        
        with open(yaml_file, 'r') as f:
            self.yaml = yaml.safe_load(f)
        
        # Append the supdict if it exists.
        if 'supdict' in kwargs:
            supdict = kwargs.get('supdict') 
            self.yaml={**self.yaml, **supdict}
        
        # Read the jobscript template into memory.
        self.jobscript = open(template_file)
        self.working_string = self.jobscript.read()
        self.jobscript.close()
        
        self.outfile=os.path.join(out_dir,jobname+suffix)
        
        
    def populate(self):
        
        """ populate
        
        Populates the jobscript with items from the yaml
        
        Returns
        ----------
        self.working_string: populated script.
        
        """
        
        
        for e in self.yaml:
            rS = re.compile(e)
            self.working_string = re.sub(rS, self.yaml[e], self.working_string)
            
    def writeout(self):
        
        """  writeout
       
        Writes out the jobscript file to the outfile location.
        
        """
        
        
        of = open(self.outfile,'w')
        of.write(self.working_string)
        of.close()
        
    def execute(self,execute_type='sbatch'):
        
        """ execute
        Executes the script.
        """
        
        os.system(execute_type + " " + self.outfile)
        print('{job} sent to SLURM'.format(job=self.jobname))



def send_slurm_job(template,yaml,supdict,jobname,jobout):
    
    exdict={'---J---':jobname,'---error---':os.path.join(jobout,jobname+'.err'),'---output---':os.path.join(jobout,jobname+'.out')}

    supdict={**supdict, **exdict}

    job=Script_Populator(yaml_file=yaml,template_file=template,out_dir=jobout,jobname=jobname,suffix='.sh',supdict=supdict)
    job.populate()
    job.writeout()
    job.execute()
        
        
class MMP_masker:
    def __init__(self,MMPloc=MMP_PATH,space='fsaverage'):
        self.MMPloc=MMPloc
        self.space=space
        self.load()
        
    def load(self):    
        self.annotfile_L = os.path.join(self.MMPloc,'{space}_lh.HCP-MMP1.annot'.format(space=self.space))
        self.annotfile_R = os.path.join(self.MMPloc,'{space}_rh.HCP-MMP1.annot'.format(space=self.space))
        
        self.lh_labels, self.lh_ctab, self.lh_names = nib.freesurfer.io.read_annot(self.annotfile_L)
        self.rh_labels, self.rh_ctab, self.rh_names = nib.freesurfer.io.read_annot(self.annotfile_R)
        self.lh_names=self.decode_list(self.lh_names)
        self.rh_names=self.decode_list(self.rh_names)
        
    def decode_list(self,inlist):
        outlist=[x.decode() for x in inlist]
        return outlist
        
    def get_roi_index(self,label,hem='L'):
        idx=self.lh_names.index('{hem}_{label}_ROI'.format(label=label,hem=hem))
        
        return idx
    
    def get_roi_verts(self,label):
        Lverts,Rverts=np.where(self.lh_labels==self.get_roi_index(label))[0],np.where(self.rh_labels==self.get_roi_index(label))[0]
        return Lverts, Rverts
    
    def downsample(self,inarray):
        
        outarray=inarray[:vertsperhem]
        return outarray
        
    def make_roi_mask(self,label,downsample=False,boolean=True,vertsperhem=10242):
        L_empty,R_empty=np.zeros(len(self.lh_labels)),np.zeros(len(self.rh_labels))
        Lverts,Rverts=self.get_roi_verts(label)
        L_empty[Lverts]=1
        R_empty[Rverts]=1
        
        if downsample==True:
            L_empty,R_empty=self.downsample(L_empty,vertsperhem=10242),self.downsample(R_empty,vertsperhem=10242)
        
        combined_mask=np.concatenate([L_empty,R_empty])
        
        if boolean==True:
             L_empty,R_empty,combined_mask=L_empty.astype(bool),R_empty.astype(bool),combined_mask.astype(bool)
            
        return L_empty, R_empty, combined_mask
    
    def make_composite_mask(self,labels):
        roimasks=np.sum([self.make_roi_mask(label) for label in labels],axis=0)        
        return roimasks

    
class retinotopy_prior:
    def __init__(self,sub='fsaverage'):
        self.sub=sub
        self.subs = npy.data['benson_winawer_2018'].subjects
    
    def get_param(self,param,downsample=True):
        LDAT=self.subs[self.sub].lh.prop(param)
        RDAT=self.subs[self.sub].rh.prop(param)
        
        if downsample==True:
            LDAT,RDAT=self.downsample(LDAT),self.downsample(RDAT)
        combined_DAT=np.concatenate([LDAT,RDAT])
        return combined_DAT
        
    def downsample(self,inarray,vertsperhem=10242):
        outarray=inarray[:vertsperhem]
        return outarray  
        

def listdir_nohidden(path):
    return [f for f in os.listdir(path) if not f.startswith('.')]



def get_k_folds(data,nsplits=5):
    
    import random
    kf = KFold(n_splits=nsplits,shuffle=True,random_state=1)
    trainfolds=list()
    testfolds=list()
    i=0
    for train_index, test_index in kf.split(data):
        
        trainfolds.append(train_index)
        testfolds.append(test_index)
        
        random.Random(i).shuffle(trainfolds[i])
        
        random.Random(i).shuffle(testfolds[i])
        i=i+1
    return trainfolds, testfolds


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)