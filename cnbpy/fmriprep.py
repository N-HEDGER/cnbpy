import pkg_resources
import yaml,os

DATA_PATH = pkg_resources.resource_filename('cnbpy', 'test/data')

yaml_file=os.path.join(DATA_PATH,'config.yml')
template_file=os.path.join(DATA_PATH,'template.sh')

class FMRIPREP:
    def __init__(self,yaml_file=yaml_file,template_file=template_file,job='myfmriprep'):
        
        self.job=job
        self.yaml_file=yaml_file
        with open(yaml_file, 'r') as f:
            self.yaml = yaml.safe_load(f)
            
        self.slurm_dict = self.yaml['slurm']
        self.l_dict= self.yaml['lpaths']
        self.s_dict= self.yaml['spaths']
        self.m_dict= self.yaml['mounts']
        
        self.jobscript = open(template_file)
        self.working_string = self.jobscript.read()
        self.jobscript.close()
        
    def reload_yaml(self):
        with open(self.yaml_file, 'r') as f:
            self.yaml = yaml.safe_load(f)
            
    def make_mount(self,lpath,spath):
        mount='-B ' + lpath+':'+spath
        return(mount)
    
    def make_mounts(self):
        
        self.mounts=[self.make_mount(self.l_dict[mount],self.s_dict[mount]) for mount in self.m_dict['st_paths2mount']]
        

        
        
        
        
        
        
        
        
    

