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

from cnbpy.pycortex_vis import webplotter

from .utils import *

import os
from .datalad import Dataset,Subject,Abide_Subject,Abide_Dataset
from .bids import BIDS
from .fmriprep import FMRIPREP
import cortex
from .preproc import Preprocessor



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
        
        """startup
        
        Provides basic info on the part of data to perform the analysis on.
        
        Parameters
        ----------
        session: session number (int)
        run: session number (int)
        roilab: ROI to serve as the CF source region (string)
        """
        
        
        self.session=session
        self.run=run
        self.roilab=roilab

    
    def load_data(self):
        
        """load_data
        
        Loads the relevant data into memory
        
        Returns
        ----------
        self.data2fit: The data to perform fitting on (numpy array: vertices X time)
        self.datadims: The dimensions of the data (tuple)
        self.TRs: The number of time points
        """
        
        self.data2fit=self.subject.load_pybest_outcomes(*[self.session,self.run])
        self.datadims=self.data2fit.shape
        self.TRs=self.datadims[-1]
        
    def get_mask(self):
        
        """get_mask
        
        Makes a boolean mas for the connective field source region.
        
        Returns
        ----------
        self.Lmask: ROI mask, left hemisphere
        self.Rmask: ROI mask Right hemisphere
        self.Fullmask: ROI mask across hemispheres.
        """
        
        
        masker=MMP_masker()
        self.Lmask,self.Rmask,self.Fullmask=masker.make_roi_mask(self.roilab)
    
    def make_subsurf(self):
        
        """make_subsurf
        
        Makes a suburface of the region defined by the mask.
        
        Returns
        ----------
        self.subsurf: Subsurface of source region.
        """
        
        
        self.subsurf=subsurface(self.pcx_subject,[self.Lmask,self.Rmask])
        self.subsurf.create()
    
    def load_lookup(self):
        
        """load_lookup
        
        Loads in the lookup table to be spliced with the conective field centres.
        
        Returns
        ----------
        self.lookup: Pandas dataframe of parameters to be spliced.
        """
        
        
        self.lookup=pd.read_csv(os.path.join(MMP_PATH,self.lookup_wildcard.format(ROI=self.roilab)))
        
        
    def make_train_test_folds(self):

        """ make_train_test_folds
        Organises the data into LOO train and test folds.
        
        Returns
        ----------
        self.data_train: Training data (list of numpy arrays).
        self.data_test Test data (list of numpy arrays).

        """
        data_train,data_test=[],[]

                
        self.trainfolds,self.testfolds=get_k_folds(self.data2fit.T,self.ksplits)

        for fold in range(len(self.trainfolds)):

            data_train.append(self.data2fit[:,self.trainfolds[fold]]) #Thus concatenate this.
            data_test.append(self.data2fit[:,self.testfolds[fold]]) # And put the left out data in the corresponding place for the test set.
        
        self.data_train=data_train
        self.data_test=data_test    
    
    
    def prepare_analysis(self):
        
        """prepare_analysis
        
        Prepares the data for analysis, using the above stages.
        
        Returns
        ----------
        self.lookup: Pandas dataframe of parameters to be spliced.
        """
        
        self.load_data()
        self.make_train_test_folds()
        self.get_mask()
        self.make_subsurf()
        self.load_lookup()
        
    
    
    def prepare_out_dirs(self):
        
        """prepare_out_dirs
        
        Prepares the output directories for storing the results.
        
        Returns
        ----------
        self.out_csv: Directory for the csv files to be output
        self.out_flat: Directory for flatmaps to be output
        self.avcsv: Name of the csv file containing fold-averaged outcomes.
        self.longcsv: Name of the csv file containing fold-wise outcomes.
        self.out_webGL: Path to the webGL viewer.
        """
        
        
        self.out_csv=os.path.join(self.out_base,self.analysis_name,'csvs',self.subject.subid)
        self.out_flat=os.path.join(self.out_base,self.analysis_name,'flatmaps',self.subject.subid)
        self.out_webGL=os.path.join(self.out_base,self.analysis_name,'webGL',self.webGL_wildcard.format(subject=self.subject.subid,session=self.session,run=self.run,ROI=self.roilab))

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
        
        Returns
        ----------
        Pandas Dataframe containing foldwise outcomes

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
        foldframe['rad']=np.radians(foldframe['spliced_prf_polar_angle'])
        foldframe['x'],foldframe['y']=pol2cart(foldframe['spliced_prf_eccentricity'],foldframe['rad'])
        
        foldframe.assign(fold=fold)

        return foldframe
                                
    def fit_all_folds(self):

        """ fit_all_folds

        Performs fitting on all folds.
        
        Returns 
        ----------
        self.fits_concat: Dataframe of fits across all folds of data.

        """

        self.all_fits=[self.fit(i) for i in tqdm(range(len(self.data_train)))]
        
        self.fits_concat = pd.concat(self.all_fits)
        self.fits_concat['vert']=self.fits_concat.index
                                
    def summarise_fits(self):

        """ Summarise_fits

        Performs averaging, or weighted averaging to make a summary dataframe across folds.
        
        Returns
        ----------
        self.av_frame: Averaged outcomes across folds.

        """

        
        self.wavs=pd.DataFrame(np.array([self.weighted_mean(self.fits_concat,var,'weights','vert') for var in self.vars2wav]).T)
        self.wavs.columns=self.vars2wav

        # Regular averaging of variance explained.
        self.avs=pd.concat([self.fits_concat.groupby('vert', as_index=False)[var].mean() for var in self.vars2av],axis=1)

        self.av_frame=pd.concat([self.wavs,self.avs],axis=1)

        # Set the data from the source region to np.nan.
        self.av_frame.loc[self.subsurf.subsurface_verts]=np.nan

        self.av_frame=self.av_frame.assign(session=self.session,run=self.run,roi=self.roilab,TRs=self.TRs)
        
        self.av_frame['spliced_prf_eccentricity'],self.av_frame['spliced_prf_polar_angle']=cart2pol(self.av_frame['x'],self.av_frame['y'])
        
    def weighted_mean(self,df, values, weights, groupby):
        df = df.copy()
        grouped = df.groupby(groupby)
        df['weighted_average'] = df[values] / grouped[weights].transform('sum') * df[weights]

        return grouped['weighted_average'].sum(min_count=1) #min_count is required for Grouper objects
    
    def saveout(self):
        
        """ saveout
         
        Saves out the averaged and foldwise outcomes.
         
        """
        
        self.av_frame.to_csv(self.avcsv)
        self.fits_concat.to_csv(self.longcsv)
        
        
        
class Plotter:

    """Plotter
    Class for plotting the analysis.
    Reads in the parts of the yaml relating to visualisation
    Reads in the csv relating to the saved outcomes.
    Saves out the flatmaps.
    """

    def __init__(self,base):
        self.yaml=os.path.join(data_path,'ABIDE_config.yml')
        self.base=base
        self.get_outcomes()
        
        
        self.dicts=['visualisation']
        
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
            


    def get_outcomes(self):

        """ get_outcomes
            Gets the saved outcomes in the form of the CSV file.
        """
        self.frame=pd.read_csv(self.base.avcsv)


    def saveout(self):

        """ saveout
        Saves out all the plots
        """

        for cv in enumerate(self.vars2plot):
            self.make_plot(i)
            self.make_plot(i,alpha=True)


    def make_plot(self,plotnum,alpha=False,save=True):

        """ make_plot
        Make a plot (not finished).
        """
        figure= plt.figure(figsize=(self.x,self.y))


        if alpha:
            filename=os.path.join(self.base.out_base,self.alpha_plot_wildcard.format(subject=self.base.subject,session=self.base.session,run=self.base.run,ROI=self.base.roilab,Outcome=self.vars2plot[plotnum]))


            #filename=os.path.join(self.base.out_flat,'alpha_'+self.vars2plot[plotnum])
            plot=alpha_plot(self.frame[self.vars2plot[plotnum]],self.frame[self.alphavar],subject=self.pycortex_subject,cmap=self.cmaps[plotnum]+'_alpha',vmin=self.lims[plotnum][0],vmax=self.lims[plotnum][1],vmin2=np.nanmin(self.frame[self.alphavar]),vmax2=self.lims2[-1],ax=figure,colorbar=self.colorbar,rois=self.rois,labels=self.labels)

            figure.savefig(filename,dpi=self.dpi)


        else:
            filename=os.path.join(self.base.out_base,self.plot_wildcard.format(subject=self.base.subject,session=self.base.session,run=self.base.run,ROI=self.base.roilab,Outcome=self.vars2plot[plotnum]))


            #filename=os.path.join(self.base.out_flat,self.vars2plot[plotnum])
            plot=basic_plot(self.frame[self.vars2plot[plotnum]],subject=self.pycortex_subject,cmap=self.cmaps[plotnum],vmin=self.lims[plotnum][0],vmax=self.lims[plotnum][1],ax=figure,colorbar=self.colorbar,rois=self.rois,labels=self.labels)
            figure.savefig(filename,dpi=self.dpi)

        plt.cla()   # Clear axis
        plt.clf()   # Clear figure
        plt.close()
        
        
    def make_webplotter(self):
        
        """ make_webplotter
        Saves out surface based outcomes
        """ 
                
        labels=self.vars2plot
        data=[np.array(self.frame[var]) for var in self.vars2plot]
        alphadat=[np.array(self.frame[self.alphavar])]
        lims=self.lims
        cmaps=self.cmaps
        
        
        mp=webplotter(data,lims,cmaps,labels,subject=self.pcx_sub,outpath=self.base.out_webGL,port=2245,pause=1,alpha=True,alphadat=alphadat,lims2=self.lims2)
        mp.prep_all_data()
        mp.make_static()

        
        
        
import time
def Download_preproc_sub(subno,local_path='/storage/basic/nh_leverhulme/scratchpad'):
    
    """Download_preproc_sub
    A set of functions for downloading and preprocessing a subject.
    """
    
    Dset=Abide_Dataset(local_path) # Make ABIDE dataset
    msub=Abide_Subject(Dset,str(subno)) # Make ABIDE subject

    msub.make_fmriprep_dataset() # Make the fmriprep datset
    msub.crawl_fmriprep_dataset() # Crawl
    msub.get_functional_data() # Get the functional data.
    msub.make_bids()
    msub.get_sess_run_combs()
    myprepoc=Preprocessor(msub) # Preprocessing
    myprepoc.make_pybest_script()
    myprepoc.pybest_job.execute(execute_type='sh') # 
    
    while not os.path.isfile(os.path.join(msub.denoised_out_dir,msub.denoised_wildcard.format(subject=msub.subject_path,ses=msub.sess_run_combs[0][0],task='rest',run=msub.sess_run_combs[0][1],hem='L'))):
        time.sleep(2)
        print('pybest not finished')
    print('completed pybest')
    
    msub.drop_fmriprep_dataset() # Drop the raw dataset.
    
    
def CF_fit_subject(subno,analysis_name='Test_analysis',local_path='/storage/basic/nh_leverhulme/scratchpad',roilab='V1'):
    
    """Download_preproc_sub
    A set of functions for performing the CF fitting on the subject.
    """
    
    # Re-establish The ABIDE subject - get the relevant info.
    Dset=Abide_Dataset(local_path)
    msub=Abide_Subject(Dset,str(subno))
    msub.make_fmriprep_dataset()
    msub.make_bids()
    msub.get_sess_run_combs()
    
    # Do the fitting and save out.
    for counter,comb in enumerate(msub.sess_run_combs):
        myan=CF_analysis(msub,analysis_name=analysis_name)
        myan.startup(*comb,roilab)
        myan.prepare_analysis()
        myan.fit_all_folds()
        myan.summarise_fits()
        myan.prepare_out_dirs()
        myan.saveout()
        mp=Plotter(myan)
        mp.make_webplotter()