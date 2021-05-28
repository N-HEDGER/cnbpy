import pkg_resources
import yaml,os

DATA_PATH = pkg_resources.resource_filename('cnbpy', 'test/data')
yaml_file=os.path.join(DATA_PATH,'config.yml')

class FMRIPREP:
    def __init__(self,yaml_file=yaml_file):
        self.yaml_file=yaml_file
        with open(yaml_file, 'r') as f:
            self.yaml = yaml.safe_load(f)
    

