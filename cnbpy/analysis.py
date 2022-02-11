import os
import yaml
import glob 
import pickle
import numpy as np

from glob import glob
import os.path as op
from prfpy.stimulus import CFStimulus
from prfpy.utils import subsurface
from prfpy.model import CFGaussianModel
from prfpy.fit import CFFitter
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import pkg_resources


from .utils import *

data_path = pkg_resources.resource_filename('cnbpy', 'test/data')
MMP_PATH = pkg_resources.resource_filename('cnbpy', 'test/data/HCP-MMP')



class CF_analysis(object):
    """AnalysisBase
    
    Base class for preparing everything for analysis. 
    Reads in the config specified on a yaml file
    Sets up all the relevant paths, given a condition,subject, task and subtask.
    Creates a vertex mask for the prf fitting.
    Reads in the data and forms cross-validation folds.
    
    """
    
    def __init__(self,subject,analysis_name='Test_analysis'):
            
        """
        
        Returns
        -------

        """
        
        self.yaml=os.path.join(data_path,'ABIDE_config.yml')
        self.analysis_name=analysis_name
        
        self.subject=subject
        
        self.dicts=['cf_analysis']
        
        for sdict in self.dicts:
            self.internalize_yaml(sdict)
                
        
    def internalize_yaml(self,subdict):

        """internalize_config_yaml

        """
        with open(self.yaml, 'r') as f:
            self.y = yaml.safe_load(f)

        cdict = self.y[subdict]
        
        setattr(self, subdict, cdict)

        for key in cdict.keys():
            setattr(self, key, cdict[key])
            
    def startup(self,session,run,roilab):
        self.session=session
        self.run=run
        self.roilab=roilab

    
    def load_data(self):
        self.data2fit=self.subject.load_pybest_outcomes(*[self.session,self.run])
        self.datadims=self.data2fit.shape
        self.TRs=self.datadims[-1]
        
    def get_mask(self):
        masker=MMP_masker()
        self.Lmask,self.Rmask,self.Fullmask=masker.make_roi_mask(self.roilab)
    
    def make_subsurf(self):
        self.subsurf=subsurface(self.pcx_subject,[self.Lmask,self.Rmask])
        self.subsurf.create()
    
    def load_lookup(self):
        self.lookup=pd.read_csv(os.path.join(MMP_PATH,self.lookup_wildcard.format(ROI=self.roilab)))
    
    def prepare_analysis(self):
        self.load_data()
        self.make_train_test_folds()
        self.get_mask()
        self.make_subsurf()
        self.load_lookup()
        
        
    def make_train_test_folds(self):

        """ make_train_test_folds
        Organises the data into LOO train and test folds.
        ----------

        """
        data_train,data_test=[],[]

                
        self.trainfolds,self.testfolds=get_k_folds(self.data2fit.T,self.ksplits)

        for fold in range(len(self.trainfolds)):

            data_train.append(self.data2fit[:,self.trainfolds[fold]]) #Thus concatenate this.
            data_test.append(self.data2fit[:,self.testfolds[fold]]) # And put the left out data in the corresponding place for the test set.
        
        self.data_train=data_train
        self.data_test=data_test

    
    
    def prepare_out_dirs(self):
        self.out_csv=os.path.join(self.out_base,self.analysis_name,'csvs')
        self.out_flat=os.path.join(self.out_base,self.analysis_name,'flatmaps')
        self.out_webGL=os.path.join(self.out_base,self.analysis_name,'webGL',self.webGL_wildcard.format(subject=self.subject,session=self.session,run=self.run,ROI=self.roilab))

        if not os.path.isdir(self.out_flat):
            os.makedirs(self.out_flat,exist_ok = True)

        if not os.path.isdir(self.out_csv):

            os.makedirs(self.out_csv,exist_ok = True)
            
        if not os.path.isdir(self.out_webGL):

            os.makedirs(self.out_webGL,exist_ok = True)
            
        self.avcsv=os.path.join(self.out_csv,self.csv_wildcard.format(subject=self.subject.subid,session=self.session,run=self.run,ROI=self.roilab))
        self.longcsv=os.path.join(self.out_csv,self.long_csv_wildcard.format(subject=self.subject.subid,session=self.session,run=self.run,ROI=self.roilab))
        
        
    def fit(self,fold):

        """ fit
        Performs the fitting for a given fold
        Outputs a pandas dataframe.
        ----------

        """

        # Define train and test folds
        mydat_train=self.data_train[fold]
        mydat_test=self.data_test[fold]

        # Define train and test data
        train_stim=CFStimulus(mydat_train,self.subsurf.subsurface_verts,self.subsurf.distance_matrix)
        test_stim=CFStimulus(mydat_test,self.subsurf.subsurface_verts,self.subsurf.distance_matrix)

        # Make model and fitter.
        model=CFGaussianModel(train_stim)
        gf=CFFitter(data=mydat_train,model=model)

        # Perform fitting
        gf.quick_grid_fit(np.array(self.cf_sizes))

        # Perform cross validation
        gf.quick_xval(test_data=mydat_test,test_stimulus=test_stim)

        # Splice the best fitting CF centres with the lookup table.
        spliced_lookup=self.lookup.iloc[gf.quick_gridsearch_params[:,0].astype(int),:]
        spliced_lookup.columns = ['spliced_' + s for s in spliced_lookup.columns]

        # Make a dataframe of the fits.
        inframe=pd.DataFrame(gf.quick_gridsearch_params,columns=self.CF_fit_columns)
        inframe['xval_R2']=gf.xval_R2
        inframe['V0']=gf.quick_vertex_centres
        inframe['hempref']=(inframe['V0']>self.subsurf.leftlim).astype(int)
        


        # We need to ensure that weights are not negative for the weighted averaging, due to the way numpy works
        xval4weights=np.copy(gf.xval_R2)
        xval4weights[xval4weights<0]=0.00001
        inframe['weights']=xval4weights

        # Concatenate fits and spliced parameters.
        foldframe=pd.concat([inframe,spliced_lookup.set_index(inframe.index)], axis=1)
        foldframe['rad']=np.deg2rad(foldframe['spliced_prf_polar_angle'])
        foldframe['x'],foldframe['y']=pol2cart(foldframe['spliced_prf_eccentricity'],foldframe['rad'])
        
        foldframe.assign(fold=fold)

        return foldframe
                                
    def fit_all_folds(self):

        """ fit_all_folds

        Performs fitting on all folds.

        ----------

        """

        self.all_fits=[self.fit(i) for i in tqdm(range(len(self.data_train)))]
                                
    def summarise_fits(self):

        """ Summarise_fits

        Performs averaging, or weighted averaging to make a summary dataframe across folds.

        ----------

        """

        self.fits_concat = pd.concat(self.all_fits)
        self.fits_concat['vert']=self.fits_concat.index
        
        self.wavs=pd.DataFrame(np.array([self.weighted_mean(self.fits_concat,var,'weights','vert') for var in self.vars2wav]).T)
        self.wavs.columns=self.vars2wav

        # Regular averaging of variance explained.
        self.avs=pd.concat([self.fits_concat.groupby('vert', as_index=False)[var].mean() for var in self.vars2av],axis=1)

        self.av_frame=pd.concat([self.wavs,self.avs],axis=1)

        # Set the data from the source region to np.nan.
        self.av_frame.loc[self.subsurf.subsurface_verts]=np.nan

        self.av_frame=self.av_frame.assign(session=self.session,run=self.run,roi=self.roilab,TRs=self.TRs)
        
    def weighted_mean(self,df, values, weights, groupby):
        df = df.copy()
        grouped = df.groupby(groupby)
        df['weighted_average'] = df[values] / grouped[weights].transform('sum') * df[weights]

        return grouped['weighted_average'].sum(min_count=1) #min_count is required for Grouper objects
        
        
        
        