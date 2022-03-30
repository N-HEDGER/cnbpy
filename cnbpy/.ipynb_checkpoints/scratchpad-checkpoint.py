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
        
    def get_freesurfer_outcomes(self):
        
        """
        Gets the required freesurfer folders
        
        Returns
        -------
        self.freesurfer_dset: The fmriprep dataset for the subject (empty until using 'get').

        """
        
        
        for folder in self.required_freesurfer_subfolders:
            self.freesurfer_dset.get(folder)
            
            
import time
def CF_fit_subject(subno,analysis_name='Test_analysis',local_path='/storage/basic/nh_leverhulme/scratchpad',roilab='V1'):

    Dset=Abide_Dataset(local_path)
    msub=Abide_Subject(Dset,str(subno))

    msub.make_fmriprep_dataset()
    msub.crawl_fmriprep_dataset()
    msub.get_functional_data()
    msub.make_bids()
    msub.get_sess_run_combs()
    #myprepoc=Preprocessor(msub)
    #myprepoc.make_pybest_script()
    #myprepoc.pybest_job.execute()
    
    while not os.path.isfile(os.path.join(msub.denoised_out_dir,msub.denoised_wildcard.format(subject=msub.subject_path,ses=msub.sess_run_combs[0][0],task='rest',run=msub.sess_run_combs[0][1],hem='L'))):
        time.sleep(2)
        print('not there yet')
            
    
    for counter,comb in enumerate(msub.sess_run_combs):
        myan=CF_analysis(msub)
        myan.startup(*comb,roilab)
        myan.prepare_analysis()
        myan.fit_all_folds()
        myan.summarise_fits()
        myan.prepare_out_dirs()
        myan.saveout()
        mp=Plotter(myan)
        mp.make_webplotter()