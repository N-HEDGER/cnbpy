import cortex
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import nibabel as nib
import pkg_resources
import yaml
import itertools
import json
import time
from matplotlib.colors import Normalize



DATA_PATH = pkg_resources.resource_filename('cnbpy', 'test/data')

def load_example(demo=1):
    fsaverage_vertsperhem=163842
    fsaverage5_vertsperhem=10242
    
    if demo==1:
        pframe=pd.read_csv(os.path.join(DATA_PATH,'retinotopy_params.csv'))
        pframe

        lhem=pframe.iloc[:fsaverage5_vertsperhem]
        rhem=pframe.iloc[fsaverage_vertsperhem:fsaverage_vertsperhem+fsaverage5_vertsperhem]
        frame=pd.concat([lhem,rhem])
        
    elif demo==2:
    
        pframe=pd.read_csv(os.path.join(DATA_PATH,'retinotopy_params.csv'))
        pframe
        fsaverage_vertsperhem=163842
        fsaverage5_vertsperhem=10242
        lhem=pframe.iloc[:fsaverage_vertsperhem]
        rhem=pframe.iloc[fsaverage_vertsperhem:]
        frame=pd.concat([lhem,rhem])
        
    elif demo==3:
        l,r=nib.load(os.path.join(DATA_PATH,'L.params_59k.func.gii')),nib.load(os.path.join(DATA_PATH,'R.params_59k.func.gii'))
        ldat,rdat=l.agg_data(),r.agg_data()
        alldat=np.array([np.concatenate([ldat[i],rdat[i]]) for i in range(len(ldat))])
        frame=pd.DataFrame(alldat.T)
        frame.columns=['Mu','Sigma','Exponent','R2','FWHM','Speech','Speech_beta', 'PCA_1', 'PCA_2', 'PCA_3']
        frame['Mu']=np.exp(frame['Mu'])
        
    return frame
        


def basic_plot(dat,vmax,subject='fsaverage',vmin=0,rois=False,colorbar=False,cmap='plasma',ax=None,labels=True):
    
    dat=np.array(dat)
    
    light=cortex.Vertex(dat,subject=subject, vmin=vmin, vmax=vmax,cmap=cmap)
    mfig=cortex.quickshow(light,with_curvature=True,with_rois=rois,with_colorbar=colorbar,with_labels=labels,fig=ax)
    return mfig

def alpha_plot(dat,dat2,vmin,vmax,vmin2,vmax2,subject='fsaverage',rois=False,labels=False,colorbar=False,cmap='nipy_spectral_alpha',ax=None):
    light=cortex.Vertex2D(dat,dat2,subject=subject, vmin=vmin, vmax=vmax,vmin2=vmin2,vmax2=vmax2,cmap=cmap)
    mfig=cortex.quickshow(light,with_curvature=True,with_rois=rois,with_colorbar=colorbar,fig=ax,with_labels=labels)
    
    
def zoom_to_roi(subject, roi, hem, margin=15.0):
    roi_verts = cortex.get_roi_verts(subject, roi)[roi]
    roi_map = cortex.Vertex.empty(subject)
    roi_map.data[roi_verts] = 1

    (lflatpts, lpolys), (rflatpts, rpolys) = cortex.db.get_surf(subject, "flat",
                                                                nudge=True)
    sel_pts = dict(left=lflatpts, right=rflatpts)[hem]
    roi_pts = sel_pts[np.nonzero(getattr(roi_map, hem))[0],:2]

    xmin, ymin = roi_pts.min(0) - margin
    xmax, ymax = roi_pts.max(0) + margin
    
    
    plt.axis([xmin, xmax, ymin, ymax])
    print([xmin, xmax, ymin, ymax])
    return


def zoom_to_rect(myrect):
    plt.axis(myrect)


    
def zoomed_plot(dat,vmin,vmax,ROI,hem,subject='fsaverage',rois=False,colorbar=False,cmap='plasma',ax=None,labels=True,alpha=False):
    
    basic_plot(dat,vmax,subject,vmin,rois,colorbar,cmap,ax,labels)
        
    zoom_to_roi(subject,ROI,hem)
    
def zoomed_alpha_plot(dat,dat2,vmin,vmax,vmin2,vmax2,ROI,hem,subject='fsaverage',rois=False,colorbar=False,cmap='plasma',ax=None,labels=True,alpha=False):
    alpha_plot(dat,dat2,vmin,vmax,vmin2,vmax2,cmap=cmap,rois=rois,labels=labels)
    zoom_to_roi(subject,ROI,hem)
    
    
def zoomed_plot2(dat,vmin,vmax,subject='fsaverage',rect=[-229.33542, -121.50809, -117.665405, 28.478895],rois=False,colorbar=False,cmap='plasma',ax=None,labels=True):
    basic_plot(dat,vmax,subject,vmin,rois,colorbar,cmap,ax,labels)
    zoom_to_rect(rect)
    
    
def zoomed_alpha_plot2(dat,dat2,vmin,vmax,vmin2,vmax2,subject='fsaverage',rect=[-229.33542, -121.50809, -117.665405, 28.478895],rois=False,colorbar=False,cmap='plasma',ax=None,labels=True):
    alpha_plot(dat,dat2,vmin,vmax,vmin2,vmax2,cmap=cmap,rois=rois,labels=labels)
    zoom_to_rect(rect)


def NormalizeData(data):
    return (data - np.nanmin(data)) / (np.nanmax(data) - np.nanmin(data))





class webplotter():
    
    
    def __init__(self,data,lims,cmaps,labels,subject,outpath='/Users/nicholashedger/Documents/tmp',port=2245,pause=1,alpha=True,alphadat=[],lims2=[]):
                
        self.port=port
        self.data=data
        self.lims=lims
        self.lims2=lims2
        self.cmaps=cmaps
        self.labels=labels
        self.subject=subject
        self.curv = cortex.db.get_surfinfo(subject)
        self.curv.vmin = np.nanmin(self.curv.data)
        self.curv.vmax = np.nanmax(self.curv.data)
        self.outpath=outpath
        self.pause=pause
        self.alpha=alpha
        self.alphadat=alphadat
        
        self.alphalabs= ['alpha_'+ var for var in self.labels]
        
        if self.alpha==True:
            self.labels=self.labels+self.alphalabs
            self.prep_curv()
            
            if len(self.lims2)==1:
                self.lims2=self.lims2*len(self.data)
                self.alphadat=self.alphadat*len(self.data)
            
        
    def prep_curv(self):
    
        curv = cortex.db.get_surfinfo(self.subject)

        curv.data = curv.data * .5

        curv = cortex.Vertex(curv.data, self.subject, vmin=-.5,vmax=.5,cmap='gray')
    
       
        
        self.curv_rgb = np.vstack([curv.raw.red.data, curv.raw.green.data, curv.raw.blue.data])
        
        
    
    def prep_alpha_data(self,data1, data2, vmin, vmax, vmin2, vmax2,cmap):
    
        vx = cortex.Vertex(data1, self.subject, cmap=cmap, vmin=vmin, vmax=vmax)
    
        # Map to RGB
        vx_rgb = np.vstack([vx.raw.red.data, vx.raw.green.data, vx.raw.blue.data])
        norm2 = Normalize(vmin2, vmax2)   
        # Pick an arbitrary region to mask out
        # (in your case you could use np.isnan on your data in similar fashion)
        alpha = np.clip(norm2(data2), 0, 1)
        
        # Alpha mask
        display_data =  self.curv_rgb * (1-alpha) + vx_rgb * alpha
    
        # Create vertex RGB object out of R, G, B channels
        return cortex.VertexRGB(*display_data, self.subject)   
    
        
    def prep_data(self,data,vmin,vmax,cmap):
            
        vx = cortex.Vertex(np.array(data),vmin=vmin,vmax=vmax,cmap=cmap,subject=self.subject)
        return vx
    

    def prep_all_data(self):
        self.vx_data=[]
        for i in range(len(self.data)):
            self.vx_data.append(self.prep_data(self.data[i],self.lims[i][0],self.lims[i][1],self.cmaps[i]))
            
        
        for i in range(len(self.data)):
            self.vx_data.append(self.prep_alpha_data(self.data[i],self.alphadat[i],self.lims[i][0],self.lims[i][1],self.lims2[i][0],self.lims2[i][1],self.cmaps[i]))
            
        self.data_dict = dict(zip(self.labels, self.vx_data))
                    
    def show(self):
        
        self.prep_all_data()
        self.handle=cortex.webshow(self.data_dict,port=self.port,labels_visible=(),overlays_visible=(),open_browser=False)
        self.handle=self.handle.get_client()
        
    def internalize_plot_yaml(self,myml):

        """internalize_config_yaml

        """
        self.yaml=myml
        with open(self.yaml, 'r') as f:
            self.y = yaml.safe_load(f)

        self.camera_dict = self.y['camera']
        
        self.animation_dict=self.y['animation']

        self.camera_dicts=[self.camera_dict[key] for key in self.camera_dict.keys()]
        
        self.size_dict = self.y['size']

        for key in self.size_dict.keys():
            setattr(self, key, self.size_dict[key])        
    
    def make_animation_dicts(self):
        self.animseq=[]
        
        for key in self.animation_dict.keys():
            vals=np.array([np.linspace(self.animation_dict[key]['start_dict'][skey],self.animation_dict[key]['end_dict'][skey],self.animation_dict[key]['nframes']) for skey in self.animation_dict[key]['start_dict'].keys()]).T

            self.animseq.append([dict(zip(self.animation_dict[key]['start_dict'].keys(),i)) for i in vals])
    
    
    def make_anim_snaps(self):
        
        for idx,cam in enumerate(self.animation_dict.keys()):
        
            for idx,cdict in enumerate(self.animseq[idx]): 
                
                
                
                cdict['camera.azimuth']= int(cdict['camera.azimuth'])
                cdict['camera.altitude']= int(cdict['camera.altitude'])
                cdict['camera.radius']= int(cdict['camera.radius'])
                
                print(cdict)
                self.handle._set_view(**cdict)
                
                for _ in range(100):
                    for k, v in cdict.items():
                        k = k.format(subject=self.subject) if '{subject}' in k else k
                        print(self.handle.ui.get(k)[0])
                        if self.handle.ui.get(k)[0] != v:
                        
                            print('waiting for', k, self.handle.ui.get(k)[0], '->', v)
                            time.sleep(0.1)
                            continue
                    break
                time.sleep(0.1)
                
                file_name=os.path.join(self.outpath,'Anim_'+cam+"_{0:03}.png".format(idx))
                self.handle.getImage(file_name, size=(self.sizex,self.sizey))
                while not os.path.exists(file_name):
                    pass
                time.sleep(0.1)
                
    
    
    def make_snaps(self,label,camera_dict):
        
        
        self.handle._set_view(**camera_dict)
        
        self.handle.draw()
        time.sleep(self.pause)
        outstring=json.dumps(camera_dict)
        figpath=os.path.join(self.outpath,label+'_'+outstring+'.png')
        self.handle.getImage(figpath, size=(self.sizex,self.sizey))
        
        
    def make_data_snaps(self,label):
        for camera_dict in self.camera_dicts:
            self.make_snaps(label,camera_dict)
            
    def make_all_snaps(self):
        for lab in self.labels:
            self.make_data_snaps(lab)
            self.handle.nextData()


    def make_static(self):
    
        cortex.webgl.make_static(outpath=self.outpath, data=self.data_dict, recache=True)








