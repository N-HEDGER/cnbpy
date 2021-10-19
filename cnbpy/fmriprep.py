import pkg_resources
import yaml,os
import re

DATA_PATH = pkg_resources.resource_filename('cnbpy', 'test/data')

yaml_file=os.path.join(DATA_PATH,'fmriprep_config.yml')
template_file=os.path.join(DATA_PATH,'template.sh')

class FMRIPREP:
    
    """FMRIPREP
    Class for running an fmriprep job using slurm.
    """
    
    
    def __init__(self,yaml_file=yaml_file,template_file=template_file,jobname='myfmriprep'):
        
        """
        Initialise the FMRIPREP class.
        
        Parameters
        ----------
        yaml_file: yaml file containing the information to be fed into the fmriprep jobscript (provided in repo).
        template_file: template file for the jobscript 
        jobname: The name given to the job.
        
        
        """
        self.jobname=jobname # Name to be given to job.
        
        self.yaml_file=yaml_file
        with open(yaml_file, 'r') as f:
            self.yaml = yaml.safe_load(f)
            
        self.slurm_dict = self.yaml['slurm']
        self.l_dict= self.yaml['lpaths']
        self.s_dict= self.yaml['spaths']
        self.m_dict= self.yaml['mounts']
        self.c_dict= self.yaml['call']
        
        # Read the jobscript template into memory.
        self.jobscript = open(template_file)
        self.working_string = self.jobscript.read()
        self.jobscript.close()
        
        # Make the output paths.
        self.make_out_paths()
        self.prep_call_dict()
        
        
    def reload_yaml(self):
        
        """
        reload_yaml
        
        Reloads the yaml file. 
        
        """
        
        with open(self.yaml_file, 'r') as f:
            self.yaml = yaml.safe_load(f)
            
    def prep_call_dict(self):
        
        """
        prep_call_dict
        
        Converts the lists into null-seperated strings so they can be read by the command line
        
        
        """
        
        
        self.c_dict['---output_spaces---']=' '.join(self.c_dict['---output_spaces---'])
        self.c_dict['---optionals---']=' '.join(self.c_dict['---optionals---'])
        
            
    def make_out_paths(self):
        
        """
        make_out_paths
        
        Defines the location for the outputs based on the base paths and the name of the job.
        
        Returns
        ----------
        self.outfile: The location for the jobscript to be output.
        self.errfile: The location for the errorfile to be output.
        self.work_path: The location for the intermediate outputs to be stored 
        self.out_path: The location for the outputs of fmriprep to be stored.
        
        """
        
        
        self.outfile=os.path.join(self.l_dict['job_path'],self.l_dict['job_wildcard'].format(jobname=self.jobname))
        self.errfile=os.path.join(self.l_dict['job_path'],self.l_dict['err_wildcard'].format(jobname=self.jobname))
        self.out_path=os.path.join(self.l_dict['outpath_base'],self.jobname)
        self.work_path=os.path.join(self.l_dict['workpath_base'],self.jobname)
        
        
        if not os.path.isdir(self.out_path):
            print(f'making output directory at {self.out_path}')
            os.mkdir(self.out_path)
            
        if not os.path.isdir(self.work_path):
            print(f'making working directory at {self.work_path}')
            os.mkdir(self.work_path)
        
        self.l_dict['---outpath---']=self.out_path
        self.l_dict['---workpath---']=self.work_path
        self.slurm_dict['---error---']=self.errfile
        
        
    def make_mount(self,lpath,spath):
        """
        make_mounts
        Converts the local and singularity paths into a mount. 
        """
        mount='-B ' + lpath+':'+spath
        return(mount)
    
    def make_mounts(self):
        
        """
        make_mounts
        
        Makes the mounts from the local locations to the singularity image.
        
        Returns
        ----------
        self.mounts: A list of mounts
        self.mount_dict: A dictionary with key ---mounts--- that contains the mounts as a string.
        
        """
        
        
        self.mounts=[self.make_mount(self.l_dict[mount],self.s_dict[mount]) for mount in self.m_dict['st_paths2mount']]
        self.mount_dict = {"---mounts---":(' ').join(self.mounts)}
        
    def populate(self,cdict):
        
        """
        populate
        
        Populates the jobscript with items from the specififed dictionary.
        
        """
        
        
        for e in cdict:
            rS = re.compile(e)
            self.working_string = re.sub(rS, cdict[e], self.working_string)
            
            
    def populate_all(self):
        
        """
        populate_all
        
        Populates the jobscript with all dictionaries.
        
        """
        
        self.make_mounts()
        self.populate(self.c_dict)
        self.populate(self.mount_dict)
        self.populate(self.l_dict)
        self.populate(self.slurm_dict)
    
    def writeout(self):
        
        """
        writeout
        
        Writes out the jobscript file to the outfile location.
        
        """
        
        
        of = open(self.outfile,'w')
        of.write(self.working_string)
        of.close()
        

        
        
        
        
        
        
        
        
    

