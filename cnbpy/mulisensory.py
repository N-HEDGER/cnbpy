import yaml
#from fracridge import fracridge
from himalaya.scoring import r2_score_split
from himalaya.kernel_ridge import MultipleKernelRidgeCV
from sklearn.utils.validation import check_is_fitted
import matplotlib.pyplot as plt
from himalaya.kernel_ridge import Kernelizer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import check_cv

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from himalaya.backend import set_backend
from .utils import generate_leave_one_run_out
from sklearn import set_config
from himalaya.kernel_ridge import ColumnKernelizer

import pkg_resources
import pandas as pd
import os
from copy import deepcopy,copy
from cnbpy.pycortex_vis import webplotter
import pickle

import matplotlib
import matplotlib.pyplot as plt
import PyQt5
from .pycortex_vis import alpha_plot
from .utils import pol2cart,cart2pol
from .surface import Subsurface
import cortex
import psutil

from himalaya.kernel_ridge._predictions import predict_weighted_kernel_ridge
from himalaya.validation import check_array
import numpy as np


DATA_PATH = pkg_resources.resource_filename('cnbpy', 'test/data')

pkg_yaml=os.path.join(DATA_PATH,'multisensory.yml')

class sensory_dm():
    
    """ sensory_dm
    Class for performing mutlisensory CF modeling.
    """

    def __init__(self,subject,yaml=pkg_yaml):
        self.subject=subject
        self.yaml=yaml
        self.load_yaml()
        self.internalize_config(self.y,'dm')
        self.internalize_config(self.y,'extra')
        self.internalize_config(self.y,'output')
        self.internalize_config(self.y,'subsurfaces')
        
                
        
        # Set up info about modality indices in the design matrix.
        self.comps_per_hem=[self.y['patches'][v]['laplacians'] for c,v in enumerate(self.modalities)]
        vals=np.insert(np.cumsum(self.comps_per_hem)*2, 0, 0)
            
        self.modality_idxs=[]    
        for v in range(len(vals)-1):
            self.modality_idxs.append(np.array(range(vals[v],vals[v+1])))  
        
    
    def load_yaml(self):
        """ load_yaml
        Loads the yaml file into memory.
        """
        
        with open(self.yaml, 'r') as f:
            self.y = yaml.safe_load(f)
            
            
    def prepare_out_dirs(self,analysis_name):
        
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
        
        self.analysis_name=analysis_name
        self.out_csv=os.path.join(self.out_base,self.analysis_name,'csvs',self.subject.subid)
        self.out_flat=os.path.join(self.out_base,self.analysis_name,'flatmaps',self.subject.subid)
        self.out_webGL=os.path.join(self.out_base,self.analysis_name,'webGL',self.webGL_wildcard.format(subject=self.subject.subid))

        if not os.path.isdir(self.out_flat):
            os.makedirs(self.out_flat,exist_ok = True)

        if not os.path.isdir(self.out_csv):

            os.makedirs(self.out_csv,exist_ok = True)
            
        if not os.path.isdir(self.out_webGL):

            os.makedirs(self.out_webGL,exist_ok = True)
            
        
            
    def internalize_config(self,y,subdict):
        """ internalize_config
        Internalises a subdictionary of the yaml file.
        """
        
        subdict = y[subdict]

        for key in subdict.keys():
            setattr(self, key, subdict[key])
            
    def make_subsurface_from_mask(self,modality,laplacians):
        """ make_subsurface
        Makes a subsurface based on mask csv file.
        """
        maskdatL=pd.read_csv(os.path.join(self.maskdir,self.maskwcard.format(hem='L',modality=modality)))
        maskdatR=pd.read_csv(os.path.join(self.maskdir,self.maskwcard.format(hem='R',modality=modality)))
        
        maskL,maskR=np.array(maskdatL['mask']),np.array(maskdatR['mask'])
        
        surf=Subsurface(self.pcx_sub,[maskL,maskR],surftype=self.surftype)
        surf.create()
        print('Making {ncomponents} laplacian components'.format(ncomponents=laplacians))
        surf.make_laplacians(laplacians)
        return surf
        
        
            
    def make_subsurface(self,labels,laplacians):
        
        """ make_subsurface
        Makes a subsurface based on the labels in the MMP parcellation.
        """
        
        print('Making subsurface for ROI {ROI}'.format(ROI='_'.join(labels)))
        mymasker=MMP_masker(self.MMPloc)
        mask=mymasker.make_composite_mask(labels) # Makes a composite mask, therefore handles things if two ROIs are provided.
        
        # If V1, we additionally curtail the mask to include only those vertices within a certain R2 and eccentricity bound.
        if labels[-1]=='V1':
            roimask=mask[-1].astype(bool)
            R2mask=np.array(self.retprior['R2'])>self.minimum_R2
            eccmask=np.array(self.retprior['eccentricity'])<self.maximum_eccentricity
            combmask=roimask * R2mask * eccmask
            mask[0]=combmask[:self.vertsperhem]
            mask[1]=combmask[self.vertsperhem:]
            
        surf=Subsurface('hcp_999999',[mask[0],mask[1]],surftype='sphere')
        surf.create()
        print('Making {ncomponents} laplacian components'.format(ncomponents=laplacians))
        surf.make_laplacians(laplacians)
        return surf
    
    
    
    def load_lookups(self):
        
        """ load_lookups
        Loads in a lookup table for each modality, these are a pandas dataframe with the number of rows equal to the number of vertcies
        """
        
        self.lookups=[pd.read_csv(os.path.join(self.lookupdir,self.lookupwcard.format(modality=v))) for c,v in enumerate(self.modalities)]
        
        #self.lookups=[pd.read_csv(self.y['patches'][v]['lookup']) for c,v in enumerate(self.modalities)]
        # Subset to the specific variables we want to splice.
        self.lookups=[self.lookups[c][self.y['patches'][v]['vars2splice']] for c,v in enumerate(self.modalities)]
        
        
    def make_subsurfaces(self,force_new=False):
        
        """ make_subsurfaces
        Makes subsurfaces based on the information given in the yaml file.
        """
        
        if force_new: # Make a new subsurface
        
            
            self.subsurfaces=[self.make_subsurface_from_mask(v,self.y['patches'][v]['laplacians']) for c,v in enumerate(self.modalities)]
            
            
        else: # Load a pre-made subsurface.
            self.subsurfaces=[]
            
            for c,v in enumerate(self.modalities):
            
                with open(os.path.join(self.surfdir,self.surfwcard.format(modality=v)), 'rb') as handle:
                    self.subsurfaces.append(pickle.load(handle))
                
    def make_roi_data(self):
        
        """ make_roi_data
        Makes ROI data seperately for each subsurface and hemisphere. 
        """
        
        self.roidat=[]

        for c,v in enumerate(self.subsurfaces):
            self.roidat.append(self.data[getattr(v,'subsurface_verts_L'),:])
            self.roidat.append(self.data[getattr(v,'subsurface_verts_R'),:])
            
    def make_eigs(self):
        
        """ make_eigs
        Gets the eigenvectors for seperately for each subsurface and hemisphere. 
        """
        
        self.eig=[]
        for c,v in enumerate(self.subsurfaces):
            self.eig.append(getattr(v,'L_eigenvectors'))
            self.eig.append(getattr(v,'R_eigenvectors'))
    
    def make_dm(self,data):
        
        """ make_dm
        Takes the dot product of each of the data sets of data in each subsurface and the respective eigenvector. 
        """ 
        
        self.data=data
        self.make_roi_data()
        self.make_eigs()
        
        self.dms=[]
        for c,v in enumerate(self.roidat):
            self.dms.append(np.dot(v.T, self.eig[c].real))
        self.dm=np.hstack(self.dms)
        
    
    def set_backend(self):
        
        """ set_backend
        Sets the backend for the modeling. 
        """ 

        
        self.backend = set_backend(self.backend_engine, on_error="warn")
        

    
    def make_cv(self,run_durations=None,ksplits=None):
        
        """ make_cv
        Make the cross validation scheme.
        """ 
        
        self.n_samples_train = self.dm.shape[0]
        
        if run_durations!=None:
            
            self.cv = generate_leave_one_run_out(self.n_samples_train, run_durations)
            self.cv = check_cv(self.cv)  
        
        if ksplits!=None:
            self.cv=int(ksplits)
        
        
    def make_preproc(self):
        
        """ make_preproc
        Make the pre-processing pipeline.
        """ 
        
        set_config(display='diagram')  # requires scikit-learn 0.23

        self.preprocess_pipeline = make_pipeline(
        StandardScaler(with_mean=self.with_mean, with_std=self.with_std),
        Kernelizer(kernel="linear")
        )
        
    def kernelize(self):
        
        """ kernelize
        Define the different feature spaces.
        """ 
    
        start_and_end = np.concatenate([[0], np.cumsum(np.array(self.comps_per_hem)*2)])
        slices = [
        slice(start, end)
        for start, end in zip(start_and_end[:-1], start_and_end[1:])]

        kernelizers_tuples = [(name, self.preprocess_pipeline, slice_)
                      for name, slice_ in zip(self.modalities, slices)]
        self.column_kernelizer = ColumnKernelizer(kernelizers_tuples)
                
        
    def setup_model(self):
        
        """ setup_model
        Setup the model object.
        """ 
        
        self.alphas=np.logspace(self.alpha_min,self.alpha_max,self.alpha_vals)
        
        solver_params = dict(n_iter=self.n_iter, alphas=self.alphas,
                     n_targets_batch=self.n_targets_batch,
                     n_alphas_batch=self.n_alphas_batch,
                     n_targets_batch_refit=self.n_targets_batch_refit)

        self.model = MultipleKernelRidgeCV(kernels="precomputed", solver=self.solver,
                                  solver_params=solver_params, cv=self.cv)
        
    def complete_pipeline(self):
        
        """ complete_pipeline
        Put together the preprocessing, kernelizing and modeling.
        """ 
        
        self.pipeline = make_pipeline(
        self.column_kernelizer,
        self.model,)
        
        self.pipeline2=deepcopy(self.pipeline)
        
    def prep_pipeline(self,run_durations=None,ksplits=None):
        
        """ prep_pipeline
        Puts the entire pipeline together.
        """ 
        
        self.internalize_config(self.y,'modeling')
        self.set_backend()
        self.make_cv(run_durations,ksplits)
        self.setup_model()
        self.make_preproc()
        self.kernelize()
        self.complete_pipeline()
    
    

    
    def fit_by_chunk_on_gpu(self,fit_data,test_data,GPU_num):
        
        """ fit_by_chunk_on_gpu
        
        Partititions the fitting into chunks. Performs the fitting and testing simultaneously since the load on the GPU is reduced by deleting the pipeline object on each iteration.
        """ 
        
        torch.cuda.set_device(GPU_num) 
        
        self.fit_data=fit_data
        self.test_data=test_data
        
        
        
        self.split_indices = np.array_split(np.arange(fit_data.shape[-1]), self.gpu_batches)
        
        self.d2fit=[]
        self.test_d2fit=[]
        
        self.pipeline_chunk=[]
        self.xfit_chunk=[]
        self.betas_chunk=[]
        self.Y_hat_split_chunk=[]
        self.split_scores_chunk=[]
        
        self.xval_score_chunk=[]
        self.Y_hat_split_test_chunk=[]
        self.test_split_scores_chunk=[]
        self.best_alphas_chunk=[]
        
       
        
        for c,v in enumerate(self.split_indices):
            self.d2fit.append(self.fit_data[:,v])
            self.test_d2fit.append(self.test_data[:,v])
            
            self.pipeline_chunk.append(deepcopy(self.pipeline))
            self.pipeline_chunk[c].fit(self.dm, self.d2fit[c])
            self.xfit_chunk.append(self.pipeline_chunk[c][-2].get_X_fit())
            self.betas_chunk.append(self.be_to_npy(self.pipeline_chunk[c][-1].get_primal_coef(self.xfit_chunk[c])))
            
            self.Y_hat_split_chunk.append(self.be_to_npy(self.pipeline_chunk[c].predict(self.dm, split=True)))
            self.split_scores_chunk.append(self.be_to_npy(r2_score_split(self.d2fit[c], self.Y_hat_split_chunk[c])))
            
            self.best_alphas_chunk.append(self.be_to_npy(self.pipeline_chunk[c][-1].best_alphas_))
            
            
            # Cross validation.
            self.xval_score_chunk.append(self.be_to_npy(self.pipeline_chunk[c].score(self.test_dm,self.test_d2fit[c])))
            self.Y_hat_split_test_chunk.append(self.be_to_npy(self.pipeline_chunk[c].predict(self.test_dm, split=True))) 
            self.test_split_scores_chunk.append(self.be_to_npy(r2_score_split(self.test_d2fit[c], self.Y_hat_split_test_chunk[c])))
            
            
            # Removing pipeline should reduce the load on the GPU
            self.pipeline_chunk[c]=None
            
        self.Y_hat_split=np.concatenate(self.Y_hat_split_chunk,axis=-1)
        self.betas=np.concatenate(self.betas_chunk,axis=-1)
        self.split_scores=np.hstack(self.split_scores_chunk)
        self.best_alphas=np.concatenate(self.best_alphas_chunk)

        self.test_split_scores=np.hstack(self.test_split_scores_chunk)
        self.Y_hat_split_test=np.concatenate(self.Y_hat_split_test_chunk,axis=-1)
            
            
    def make_test_dm(self,surf_data):
        
        """ make_test_dm
        tests the model on out of sample data.
        """ 
        
        train_roidat=copy(self.roidat) # Copy the training roidata.
        train_data=copy(self.data)
        
        self.data=surf_data
        self.make_roi_data() # Make new ROI data
        
        # Make a new design matrix, ensuring that the current eigenvectors are used.
        test_dms=[]
        for c,v in enumerate(self.roidat):
            test_dms.append(np.dot(v.T, self.eig[c].real))
        self.test_dm=np.hstack(test_dms)
        
        self.roidat=train_roidat # Relace the roidat.
        self.data=train_data
        
            
        
    
    def fit(self,fit_data):
        
        """ fit
        Performs the fitting.
        """ 
        
        self.fit_data=fit_data
        self.pipeline.fit(self.dm, self.fit_data)
        
        
        
    def get_params(self):
        """ get_params
        Gets some params from the fitting
        """ 
        
        self.xfit=self.column_kernelizer.get_X_fit()
        self.betas=self.be_to_npy(self.pipeline[-1].get_primal_coef(self.xfit)) # Get the beta weights.
        
        
        self.Y_hat_split = self.be_to_npy(self.pipeline.predict(self.dm, split=True)) # Get the split predictions (training data)
        self.split_scores = self.be_to_npy(r2_score_split(self.fit_data, self.Y_hat_split)) # Get the split R2 (training data).
        self.best_alphas=self.be_to_npy(self.pipeline[-1].best_alphas_)
        
    def test_xval(self,surf_data,data):
        
        """ test_xval
        tests the model on out of sample data.
        """ 
        
        train_roidat=copy(self.roidat) # Copy the training roidata.
        train_data=copy(self.data)
        
        self.data=surf_data
        self.make_roi_data() # Make new ROI data
        
        # Make a new design matrix, ensuring that the current eigenvectors are used.
        test_dms=[]
        for c,v in enumerate(self.roidat):
            test_dms.append(np.dot(v.T, self.eig[c].real))
        test_dm=np.hstack(test_dms)
        
        self.xval_score= self.be_to_npy(self.pipeline.score(test_dm,data))
        

        self.Y_hat_split_test=self.be_to_npy(self.pipeline.predict(test_dm, split=True)) 
        self.test_split_scores = self.be_to_npy(r2_score_split(data, self.Y_hat_split_test)) 
        
        self.roidat=train_roidat # Relace the roidat.
        self.data=train_data
        
        
    def reconstruct_profiles(self,zscore=False):
        
        """ reconstruct_profiles
        Reconstructs the CF profiles by dotting the eigenvectors with the betas.
        """ 
      
            
        self.modality_profiles=[np.vstack([np.dot(self.subsurfaces[c].L_eigenvectors.real,self.betas[c][:self.comps_per_hem[c]]),np.dot(self.subsurfaces[c].R_eigenvectors.real,self.betas[c][self.comps_per_hem[c]:])]) for c,v in enumerate(self.modalities)]
        
        
        
    def get_hempref(self):
        
        """ get_hempref
            Gets the peak hemisphere from the modality profiles.
            Defined in 2 ways i) which hemisphere the peak is in ii) The difference of the sum of each profile.
        """ 
        
        # Binary of whether the peaks is in the left or right subsurface
        hempref_peak=[(np.nanargmax(self.modality_profiles[c],axis=0)>self.subsurfaces[c].subsurface_verts_L.shape[-1]).astype(int) for c,v in enumerate(self.modalities)]
        
        # The difference in the sum of each modality profile.
        hempref_diff=[np.sum(self.modality_profiles[c][self.subsurfaces[c].subsurface_verts_L.shape[-1]:,:],axis=0)-np.sum(self.modality_profiles[c][:self.subsurfaces[c].subsurface_verts_L.shape[-1],:],axis=0) for c,v in enumerate(self.modalities)]
        
        self.hempref_frames=[pd.DataFrame(np.vstack([hempref_diff[c],hempref_peak[c]]).T,columns=['hempref_'+ p + v for p in ['diff_','peak_']]) for c,v in enumerate(self.modalities)]

        self.spliced_lookups=[pd.concat([self.spliced_lookups[c].reset_index(drop=True),self.hempref_frames[c]],axis=1) for c,v in enumerate(self.modalities)]
        
    def quantify_size(self,profiles,distance_matrix):
        
        """ quantify_size
        Quantifies the sizes of a set of modality profiles
        This is defined as the maximum distance from the peak at which the profile is above half maximum
        """ 
        
        # Find the  peak of each profile.
        peaks=np.nanargmax(profiles,axis=0)

        # Get peak value of each profile
        peakvals=np.nanmax(profiles,axis=0)

        # Create a mask corresponding to parts of the profile that are larger than the HM.
        ismore=profiles>(peakvals/2)[np.newaxis,:]
        abovehm=np.sum(ismore,axis=0)

        # Index distances matrix relative to the peak of each profile. 
        peakdistmat=distance_matrix[peaks,:]

        # Distances less than the HM are nan.
        peakdistmat[np.invert(ismore.T)]=np.nan

        # 0 distances are also nan. 
        peakdistmat[peakdistmat==0]=np.nan

        # Distances in the opposite hemisphere are nan.
        peakdistmat[peakdistmat==np.inf]=np.nan

        # what is the maximum distance where the profile is above its half maximum?
        maxdist_above_hm=np.nanmax(peakdistmat,axis=1)
        return maxdist_above_hm
    
    
    def quantify_sizes(self):
        
        """ quantify_sizes
        Quantifies the sizes of the modality profiles
        """ 
        
    
        self.size_frames=[pd.DataFrame(self.quantify_size(self.modality_profiles[c],self.subsurfaces[c].distance_matrix),columns=['size_' + v]) for c,v in enumerate(self.modalities)]
        
        self.spliced_lookups=[pd.concat([self.spliced_lookups[c].reset_index(drop=True),self.size_frames[c]],axis=1) for c,v in enumerate(self.modalities)]
    
        
    def save_outcomes(self):
        
        """ save_outcomes
        Saves out beta, alpha, split scores for test and train
        """ 
        
        self.beta_names=[]
        self.score_names=[]
        self.Yhat_names=[]
        
        for e,v in enumerate(self.modalities):
            self.beta_names.append(['{modality}_{x}_beta'.format(modality=v,x=str(x).zfill(3)) for x in range(len(self.betas[e]))])
            
            #self.Yhat_names.append(['{modality}_{x}_Yhat'.format(modality=v,x=str(x).zfill(3)) for x in range(len(self.Y_hat_split_test[e]))])
            
            self.score_names.append('{modality}_score'.format(modality=v))
            self.saveout(self.betas[e],'betas_{modality}'.format(modality=v),self.beta_names[e])
            if self.save_yhat:
                self.saveout(self.Y_hat_split_test[e],'Yhat_{modality}'.format(modality=v),self.Yhat_names[e])
                
                    
        self.saveout(self.split_scores,self.train_scorename,self.score_names)
        #self.saveout(self.test_split_scores,self.test_scorename,self.score_names)
        self.saveout([self.best_alphas],self.alphaname,[self.alphaname])
        
    def be_to_npy(self,var):
        
        """ be_to_npy
        Converts the backend format to numpy - but also handles lists.
        """ 
        
        if type(var)==list:
            res=[self.backend.to_numpy(el) for el in var]
        else:
            res=self.backend.to_numpy(var)
        return res
    
    def saveout(self,dlist,pname,nlist=None):
        
        """ saveout
        Saves out parameters to a cifti file.
        """ 
        

        fname=os.path.join(self.out_csv,self.out_npy_wildcard.format(param=pname))
        np.save(fname,dlist)

    
    
    def make_mean_profiles(self):
        
        """ make_mean_profiles
        Makes the mean profile across eigenvectors.
        """ 
            
        self.mean_modality_profiles=[np.concatenate([np.mean(self.subsurfaces[c].L_eigenvectors.real,axis=1),np.mean(self.subsurfaces[c].R_eigenvectors.real,axis=1)]) for c,v in enumerate(self.modalities)]
        
        
    def regress_out(self,profiles, mean_profile):
        
        """ regress_out
        Regresses out a pattern from the profiles. 
        """
        # Make design matrix with intercept.
        dm = np.vstack([np.ones(profiles.shape[-1]), mean_profile]).T
        # Get betas.
        betas = np.linalg.lstsq(dm, profiles.T)[0]
    
        # Get prediction
        yhat = np.dot(dm, betas).T
    
        # Get residuals
        resid=profiles-yhat
    
        # Add back intercept
        resid+= np.mean(profiles,axis=-1)[:,np.newaxis]
    
        return resid
        
        
    def regress_out_mean_profiles(self):
        
        """ regress_out_mean_profile
        Regresses out the mean modality profiles.
        """ 
        
        
        if not hasattr(self,'mean_modality_profiles'):
            self.make_mean_profiles()
            
        
        self.modality_profiles=[self.regress_out(self.modality_profiles[c].T,self.mean_modality_profiles[c]).T for c,v in enumerate(self.modalities)]
        
    
    
    def splice_lookup(self,lookup_wb,subsurface,modality_profiles,label,dot_product=False,pos_only=True):
        
        """ splice_lookup
        Splices the modality profiles with the lookup table.
        """ 
        
        lookup=lookup_wb.loc[np.concatenate([subsurface.subsurface_verts_L,subsurface.subsurface_verts_R])]
        
        colnames=lookup.columns
        
        
        if dot_product==True:
            
            weights2=np.copy(modality_profiles) # Copy the modality profiles
    
            if pos_only==True:
                weights2[weights2<0]=0 # Set all negative weights to zero.
        
            else: # Else, raise them all so that they are >0.
                weights2=weights2+np.abs(np.min(weights2,axis=0))[np.newaxis,:]
        
            summed_weights=np.sum(weights2,axis=0)
            dp=np.dot(np.array(lookup.T),weights2)
            wav=dp/summed_weights
            spliced_lookup=pd.DataFrame(wav.T,columns=[s +'_'+label for s in colnames])
            
        else:
            peaks=np.nanargmax(modality_profiles,axis=0)
            spliced_lookup=lookup.iloc[np.array(peaks),:]
            spliced_lookup.columns = [s +'_'+label for s in colnames]
                                                                     
        return spliced_lookup
        
    def splice_lookups(self):
        
        self.internalize_config(self.y,'splicing')
        
        self.calc_size=False
        
        """ splice_lookups
        Splices all modality profiles with all lookup tables.
        """ 
        if self.regress_out_mean==True:
            self.regress_out_mean_profiles()
        
        self.spliced_lookups=[self.splice_lookup(self.lookups[c],self.subsurfaces[c],self.modality_profiles[c],v,dot_product=self.dot_product[c],pos_only=self.pos_only[c]) for c,v in enumerate(self.modalities)]
        
        if self.calc_hempref:
            self.get_hempref()
        if self.calc_size:
            self.quantify_sizes()
        
        
    def save_spliced_lookups(self):
        
        """ save_spliced_lookups
        Saves out the spliced lookup tables.
        """ 
        
        self.subject.brainmodel=Ciftihandler()
        
        joint_names=[self.spliced_lookups[c].columns.tolist() for c,v in enumerate(self.spliced_lookups)]
        joint_names = [val for sublist in joint_names for val in sublist]
        joint_array=np.concatenate(self.spliced_lookups,axis=1).T
            
        self.saveout(joint_array,self.spliced_paramname,joint_names)
            
        
    def load_param(self,pname,split=False,npy=False):
        
        """ load_param
        Loads a parameter saved out.
        """ 
        
        if self.npy==True:
            fname=os.path.join(self.out_csv,self.out_npy_wildcard.format(param=pname))
            return np.load(fname)
        
        else:
            fname=os.path.join(self.out_csv,self.out_cifti_wildcard.format(param=pname))
            bm=Ciftihandler(fname)
            bm.get_data()
            if split==True:
                split_dat=bm.decompose_data(bm.data)
                return split_dat
            else:
                return bm.data
        
        
    
    def load_precomputed_betas(self):
        
        """ load_precomputed_betas
        Loads precomputed betas into memory.
        Useful for when we want to splice a lookup table after the fitting has been performed.
        """ 
        
        self.betas=[self.load_param('betas_{modality}'.format(modality=v)) for e,v in enumerate(self.modalities)]
    
    def make_flat_mapping(self):
        
        """ make_flat_mapping
            Makes a flat mapping for the subject.
        """ 
        
        (linfpts, linfpolys), (rinfpts, rinfpolys) = cortex.db.get_surf(self.pcx_sub, "flat",
                                                                nudge=True)
        mapping_flat=np.vstack([linfpts,rinfpts])
        self.flat_mapping=pd.DataFrame(mapping_flat,columns=['flat_dim1','flat_dim2','flat_dim3'])
    
    def make_subsurface_flat_mappings(self):
        
        """ make_subsurface_flat_mappings
            Subsets the flat mapping to the subsurface.
        """ 
        
        self.make_flat_mapping()

        self.subsurface_flat_mappings=[self.flat_mapping.loc[self.subsurfaces[c].subsurface_verts].reset_index() for c,v in enumerate(self.subsurfaces)]
        
        
        
 
    def prepare_frame(self,with_spliced=False):
    
        self.internalize_config(self.y,'plotting')
        
        print('yo')
        
        data=[]
        labels=[]
        for param in self.vars2plot:
                
            dat=self.load_param(param)
    
            if 'scores' in param:
                nlist=['{var}_{mod}_score'.format(var=param,mod=mod) for mod in self.modalities]
                for n in nlist:
                    labels.append(n)
            for d in range(dat.shape[0]):
                # Need to log alphas
                if param=='best_alphas':
                    data.append(np.log10(dat[d,:]))
                    labels.append(param+'_'+param)
                else:
                    data.append(dat[d,:])
            
        frame=pd.DataFrame(np.array(data).T,columns=labels)
        self.frame=frame.assign(analysis_name=self.analysis_name,subject=self.subject.subid)

        if with_spliced==True:
            joint_names=[self.spliced_lookups[c].columns.tolist() for c,v in enumerate(self.spliced_lookups)]
            joint_names = ['spliced_params_'+val for sublist in joint_names for val in sublist]
            joint_array=np.concatenate(self.spliced_lookups,axis=1).T


            self.frame=pd.concat([self.frame,pd.DataFrame(joint_array.T,columns=joint_names)],axis=1) 
            
    
    def save_frame(self,with_spliced=False):
        
        if not hasattr(self,'frame'):
            self.prepare_frame(with_spliced)
            
        self.frame.to_csv(os.path.join(self.out_csv,self.out_csv_wildcard))

    

        