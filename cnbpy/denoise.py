from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import os
import nibabel as nib
import scipy



class PCA_denoiser():
    """
    PCA_denoiser
    Class for preparing nuisance regressors and then running PCA denoising on fMRIprep output using confounds.tsv file  
    
    """
    
    def __init__(self,csvpath='/content/data/confounds.tsv',datapath='/content/data/data.gii'):
        """
        Initialise the PCA_denoiser class.
        
        Input arguments
        ----------
        csv_path: path to the .tsv file containing noise components
        data_path: path to fMRIprep output
        
        Parameters
        ----------
        
        self.csvpath: the path to the noise components tsv file.
        self.datapath: the path to the fMRIprepped data.
        self.load_csv: reads the .tsv file using pandas
        self.load_data: loads the fMRIprepped data using nibabel
        
        """
      self.csvpath=csvpath
      self.datapath=datapath
      self.load_csv()
      self.load_data()

    def load_csv(self):
      self.confound_frame = pd.read_csv(self.csvpath, sep='\t', header=0,index_col=None)

    def load_data(self):
      data=nib.load(self.datapath)
      self.data=data.agg_data()

    def subset_frame(self,vars):
         """ subset_frame
        
        Subsets the noise components to only include variables defined by the user
        
        Returns
        ----------
        self.subsetted_frame: subsetted noise components.
        
        """
      self.subsetted_frame=self.confound_frame[vars]

    def prepare_frame(self):
        """ prepare_frame
        
        Prepares the noise regressors.

        1.   Converts this subsetted frame to a numpy array. 
        2.   Converts NAs for each regressor to the median.
        3. Zscores each regressor over time.
        4. Returns the result in self.prepared_array

     
        Returns
        ----------
        self.nuissance_array: a numpy array of nuisance regressors with NANs replaced with medians.
        self.prepared_array: a z-scored numpy array of the nuisance regressors
        """
        self.nuissance_array=np.array(self.subsetted_frame)
        medians=np.nanmedian(self.nuissance_array,axis=0)
        for c,v in enumerate(medians):
          self.nuissance_array[:,c][np.isnan(self.nuissance_array[:,c])]=medians[c]
        self.prepared_array=scipy.stats.zscore(self.nuissance_array)
    
    def PCA_regression(self,ncomps):
        """ PCA_regression
        
        Runs a PCA on the noise components and removes these from the data.

        1. Fits a PCA of n components to self.prepared_array
        2. Creates a design matrix with these regressors and an intercept.
        3. Fits this to self.data
        4. Takes the residuals of the regression.
        5. Adds back in the intercept.
        6. Returns the result in self.denoised_data
     
        Returns
        ----------
       
        self.pca_comps: PCA output of 
        self.dm: design matrix (including intercept)
        self.betas: Ordinary least squares beta coefficients for respective noise components
        self.yhat: predictions (dot products) 
        self.rsq: root mean squared differences
        self.resid: residuals after removing PCA components
        self.denoised_data: residuals (now cleaned data)
        """
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

      # Add back in intercept
      self.resid+= np.nanmean(self.data,axis=-1)[:,np.newaxis]
      self.denoised_data=self.resid