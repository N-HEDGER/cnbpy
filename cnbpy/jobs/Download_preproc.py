from cnbpy.analysis import Download_preproc_sub
import sys


subj=str(sys.argv[1])
local_path=str(sys.argv[2])

Download_preproc_sub(subno=subj,local_path=local_path)