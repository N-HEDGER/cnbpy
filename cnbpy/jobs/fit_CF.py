from cnbpy.analysis import CF_fit_subject
import sys



subj=str(sys.argv[1])
analysis_name=str(sys.argv[2])
local_path=str(sys.argv[3])
roilab=str(sys.argv[4])


CF_fit_subject(subno=subj,analysis_name=analysis_name,local_path=local_path,roilab=roilab)

