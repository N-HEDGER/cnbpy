import os

import nest_asyncio
nest_asyncio.apply()

import datalad.api as dl
import pandas as pd
import pkg_resources


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
        
        
    def remove(self):
        
        """
        Removes the dataset using http://docs.datalad.org/en/stable/generated/man/datalad-remove.html
        """
        
        self.dset.remove()
        

        
                
        
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
        self.subject_path=self.subprefix+str(subid)
        
    def get(self):
        
        """
        Downloads the subject data using http://docs.datalad.org/en/stable/generated/man/datalad-get.html

        """
        self.DATASET.dset.get(self.subject_path,get_data=True)
        
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
        
    
        self.demofile=pd.read_csv(os.path.join(DATA_PATH,'ABIDEII_Composite_Phenotypic.csv'), encoding = "ISO-8859-1")
        self.anat_q=pd.read_csv(os.path.join(DATA_PATH,'anat_qap.csv'), encoding = "ISO-8859-1")
        self.func_q=pd.read_csv(os.path.join(DATA_PATH,'functional_qap.csv'), encoding = "ISO-8859-1")
        
        super().__init__(local_base,source)
    
