class OPENNEURO:
    
    """OPENNEURO
    
    A class for interacting with OPENNEURO DATA.
        
    """
    
    
    # Local dir is the location where we want our downloads to be stored.
    def __init__(self,local_base,dset='ds000228/ds000228_R1.0.1/uncompressed'):
        
        """__init__
        constructor for OPENNEURO class, takes stimulus object as argument
        
        Parameters
        ----------
        local_base : Local path in which to sync files. 
        """
        
        
        self.dset=dset
        self.connection=s3fs.S3FileSystem(anon=True)
        self.rpath=os.path.join('openneuro',self.dset)
        self.rlist=self.connection.ls(self.rpath)



class subject_defunct:
    def __init__(self,subno,ABIDE):

        self.sub='sub-'+str(subno)
        self.row=np.where(ABIDE.demofile['SUB_ID']==subno)[0][0]
        self.site=ABIDE.demofile['SITE_ID'][self.row]
        self.raw_path=os.path.join(ABIDE.raw_path,self.site,self.sub)
        self.fmriprep_path=os.path.join(ABIDE.fmriprep_path,self.sub)
        self.freesurfer_path=os.path.join(ABIDE.freesurfer_path,self.sub)
        self.mriqc_path=os.path.join(ABIDE.mriqc_path,self.sub)

        self.rpaths=[self.raw_path,self.fmriprep_path,self.freesurfer_path]

        self.local_dir=os.path.join(ABIDE.local_base,self.sub)
        self.l_raw_path=os.path.join(self.local_dir,'raw_dat')
        self.l_fmriprep_path=os.path.join(self.local_dir,'fmriprep')
        self.l_freesurfer_path=os.path.join(self.local_dir,'freesurfer')

        self.lpaths=[self.l_raw_path,self.l_fmriprep_path,self.l_freesurfer_path]


        self.fmriprep_rep=os.path.join(ABIDE.aws_suffix,ABIDE.fmriprep_path,self.sub+'.html')

        mriqclist=ABIDE.connection.ls(ABIDE.mriqc_path)
        self.mriqc_reps = [s for s in mriqclist if self.sub in s and 'html' in s]

    def show_report(self,ABIDE):

        for html in self.mriqc_reps:
            webbrowser.open_new_tab(os.path.join(ABIDE.aws_suffix,html))
        webbrowser.open_new_tab(self.fmriprep_rep)

        return

    def download(self,ABIDE):

        os.mkdir(self.local_dir)

        for path in self.lpaths:
            os.mkdir(path)

        for i in tqdm(range(len(self.rpaths))):
            ABIDE.connection.get(self.rpaths[i],self.lpaths[i],recursive=True)


        return


    

    

class ABIDE_defunct:
    # Local dir is the location where we want our downloads to be stored.
    def __init__(self,local_base):
        self.demofile=pd.read_csv(os.path.join(DATA_PATH,'ABIDEII_Composite_Phenotypic.csv'), encoding = "ISO-8859-1")
        self.connection=s3fs.S3FileSystem(anon=True)
        self.fmriprep_path='fcp-indi/data/Projects/ABIDE2/Outputs/fmriprep/fmriprep'
        self.freesurfer_path='fcp-indi/data/Projects/ABIDE2/Outputs/fmriprep/freesurfer'
        self.mriqc_path='fcp-indi/data/Projects/ABIDE2/Outputs/mriqc/mriqc'
        self.raw_path='fcp-indi/data/Projects/ABIDE2/RawData'
        self.aws_suffix='https://s3.amazonaws.com/'
        self.anat_q=pd.read_csv(os.path.join(DATA_PATH,'anat_qap.csv'), encoding = "ISO-8859-1")
        self.func_q=pd.read_csv(os.path.join(DATA_PATH,'functional_qap.csv'), encoding = "ISO-8859-1")
        self.local_base=local_base
        
        

        
        
    # Shows a file.
    def show(self,npath):

        niftifiledim=len(image.load_img(npath).shape)
        if niftifiledim == 3:
           display=plotting.view_img(npath,bg_img=False)
        else:
            print ('>1 volume, plotting only the first for convenience')
            firstim=image.index_img(npath, 0)
            display=plotting.view_img(firstim,bg_img=False)
        return display


    def getniftibits(self,npath):
        nifti = nib.load(npath)
        VOXSIZE = nifti.header['pixdim'][1:4]
        SHAPE= (nifti.header['dim'][1:5])
        TR = (nifti.header['pixdim'][4:5])
        VOXFRAME=pd.DataFrame(VOXSIZE)
        VOXFRAME=VOXFRAME.T
        SHAPEFRAME=pd.DataFrame(SHAPE)
        SHAPEFRAME=SHAPEFRAME.T
        VOXFRAME.columns=['VoxsizeX','VoxsizeY','VoxsizeZ']
        SHAPEFRAME.columns=['ShapeX','ShapeY','ShapeZ','Volumes']
        CFRAMEi=pd.concat([VOXFRAME,SHAPEFRAME],axis=1)
        CFRAMEi['TR'] = TR
        return(CFRAMEi)
