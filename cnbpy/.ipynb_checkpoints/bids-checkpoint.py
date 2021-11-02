import os.path as op
from glob import glob
import pandas as pd
import json



class BIDS:
    
    """BIDS
    Class for a BIDS dataset.
    """
    
    def __init__(self,path,func_fmt='.nii'):
        
        
        """
        Initialise the BIDS class.
        
        Parameters
        ----------
        path: Path to the directory where the dataset is stored (string).
        func_fmt: What filetype identifies the functional data (optional, default= '.nii' string)
        
        Returns
        -------
        self.subjects: Subjects in the BIDS directory (list).
        self.submessage: Message that prints what subjects were found.
        self.sessions: What tasks were found per session, per subject (list).
        self.sessmessage: Message that prints what sessions were found .
        self.tasks: What tasks were found per session, per subject (list).
        self.taskmessage: Message that prints the tasks found.

        """
        
        
        self.path=path
        self.subjects=[op.basename(s).split('-')[1] for s in
            sorted(glob(op.join(self.path, 'sub-*'))) if op.isdir(s)]

        self.submessage=f"Found {len(self.subjects)} participant(s) {self.subjects}"

        sessions = []
        sessmessage=[]
        for this_sub in self.subjects:
            these_ses = [
                op.basename(s).split('-')[1] for s in
                sorted(
                    glob(op.join(self.path, f'sub-{this_sub}', 'ses-*')))
                if op.isdir(s)
            ]

            sessmessage.append(f"Found {len(these_ses)} session(s) for sub-{this_sub} {these_ses}")
            sessions.append(these_ses)

        self.sessions=sessions
        self.sessmessage=sessmessage

        tasks=[]
        taskmessage=[]
        for this_sub, these_ses in zip(self.subjects, self.sessions):
            these_task = []
            for this_ses in these_ses:
                if this_ses is None:
                    tmp = glob(op.join(
                        self.path,
                        f'sub-{this_sub}',
                        'func',
                        f"*{func_fmt}*"
                    ))
                else:
                    tmp = glob(op.join(
                        self.path,
                        f'sub-{this_sub}',
                        f'ses-{this_ses}',
                        'func',
                        f"*{func_fmt}*"
                    ))

                these_ses_task = list(set(
                    [op.basename(f).split('task-')[1].split('_')[0]
                     for f in tmp]
                ))


                nullstring = "" if this_ses is None else f"and ses-{this_ses}"

                taskmessage.append(f"Found {len(these_ses_task)} task(s) for sub-{this_sub} {nullstring} {these_ses_task}")
                these_task.append(these_ses_task)



            self.taskmessage=taskmessage
            tasks.append(these_task)
        self.tasks=tasks
        
        self.load_all_jsons()
    


        
    def elaborate(self):
        
        """
        Prints out all the messages from the BIDS initialisation.

        """
        
        print(self.submessage)
        print(self.sessmessage)
        print(self.taskmessage)

    def find_funcs(self,*,sub,task,ses=None,ext='.nii.gz'):
        
        """
        Finds functional files for a given subject, session and task.
        Keywords are forced. Not positional.
        
        Parameters
        ----------
        sub: Subject identifier (string)
        ses: Session identifier (string)
        task: Task identifier (string)
        
        Returns
        -------
        funcs: Functional files that meet the criteria (list).
        

        """
    
        if ses != None:

            ffunc_dir = op.join(self.path, f'sub-{sub}', f'ses-{ses}', 'func')
            
        else:
            ffunc_dir = op.join(self.path, f'sub-{sub}', 'func')
            
        print(ffunc_dir)

        funcs = sorted(glob(op.join(ffunc_dir, f'*task-{task}*')))
        funcs =[x for x in funcs if ext in x]
        if not funcs:
            raise ValueError(
                "Could not find functional data with the following parameters:\n"
                f"sub-{sub}, ses-{ses}, task-{task}")
        return funcs

    def find_anats(self,*,sub,ses=None,weight='T1'):
        
        """
        Finds anatomical files for a given subject, session.
        
        Parameters
        ----------
        sub: Subject identifier (string)
        ses: Session identifier (string)
        weight: What type of weight (optional, default ='T1')
        
        Returns
        -------
        anats: Anatomical files that meet the criteria (list).
        
        """
        
        if ses != None:
            
            anat_dir = op.join(self.path, f'sub-{sub}', f'ses-{ses}', 'anat')
            
        else: 
            anat_dir = op.join(self.path, f'sub-{sub}', 'anat')

        anats = sorted(glob(op.join(anat_dir, f'*{weight}*')))

        if not anats:
            raise ValueError(
                "Could not find anatomical data with the following parameters:\n"
                f"sub-{sub}, ses-{ses}, weight-{weight}")

        return anats
    
    
    def find_jsons(self):
        self.jsons=glob(op.join(self.path, '*.json'))
        
    def load_json(self,file):
        with open(file) as f:
            key='json_'+op.split(file)[-1][:-5]
            setattr(self, key, json.load(f))  
                     
    def load_all_jsons(self):
        self.find_jsons()
        [self.load_json(json) for json in self.jsons] 
    


