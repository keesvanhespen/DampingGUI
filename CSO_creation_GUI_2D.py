from mevis import *
import numpy as np
import os
import sys
import matplotlib.pyplot as plt 
import datetime
import pydicom
import subprocess as subp
import scipy.io as sio

from subprocess import Popen

def init():
  ctx.field("XMarkerListContainer.deleteAll").touch()
  ctx.field("ROISelect.startVoxelX").setValue(-1000)
  ctx.field("ROISelect.startVoxelY").setValue(-1000)
  ctx.field("ROISelect.endVoxelX").setValue(-1000)
  ctx.field("ROISelect.endVoxelY").setValue(-1000)
  ctx.field("CSOManager2.removeAllCSOsAndGroups").touch()
  
  
def contourcreation():
  global Results_and_settings
  # Creates iso contours at the FWHM value, and initializes the graphs on which the velocity curves are drawn.
  ctx.field("CSOManager2.removeAllCSOsAndGroups").touch()

  #Take maximum of magnitude image across all time points
  magnitude = ctx.field("ROISelect.outImage").image()
  magnitude = magnitude.getTile( (0,0,0,0,0), (magnitude.UseImageExtent,magnitude.UseImageExtent,magnitude.UseImageExtent,magnitude.UseImageExtent,magnitude.UseImageExtent) )   
  img = ctx.field("ROISelect.outImage").image()
  img = img.getTile( (0,0,0), (img.UseImageExtent,img.UseImageExtent,img.UseImageExtent) )
  cso_areas = np.zeros((magnitude.shape[0],img.shape[0]))
  for timepoint in range(0,magnitude.shape[0]):
    ctx.field("SubImage.t").setValue(timepoint)
    #Generate contours
    ctx.field("CSOIsoGenerator2.timePoint").setValue(timepoint)
    ctx.field("itkOtsuThresholdImageFilter.update").touch()     
    ctx.field("RunPythonScript.execute").touch()
    ctx.field("CSOIsoGenerator2.startPosition").setStringValue(ctx.field("XMarkerAtIndex.position3D").stringValue())
    ctx.field("CSOIsoGenerator2.addCSOToGroupWithLabel").setValue(0)
    ctx.field("CSOIsoGenerator2.apply").touch()
    
    #check contour area
    csolist = ctx.field("CSOManager2.outCSOList").object()
    ctx.field("CSOInfo.csoShowByIndexOrId").setValue(csolist.getCSOIdList()[-1])
    cso_areas[timepoint,0] = ctx.field("CSOInfo.csoArea").value
  #Check if there are outliers with low areas. These could be incorrect contours in flow voids.
  d = np.abs(cso_areas- np.median(cso_areas[~np.all(cso_areas == 0, axis=1),:],axis=0))
  axe=1
  mdev = np.median(d[cso_areas!=0])
  s = d/mdev if mdev else 0.
  outliers = s<4
  np.all(outliers,axis=axe)
  np.argwhere(outliers == 0)
  flow_void_candidates = np.argwhere(outliers == 0)
  if np.squeeze(np.argwhere(np.all(cso_areas == 0,axis=0))).size:
    flow_void_candidates = flow_void_candidates[np.invert(np.isin(flow_void_candidates[:,1],np.squeeze(np.argwhere(np.all(cso_areas == 0,axis=0))))),:]
    

  if flow_void_candidates.size!=0:
    #try to fix flow void candidates
    for candidate in range(0, flow_void_candidates.shape[0]):
      timepoint,slice = flow_void_candidates[candidate,0],flow_void_candidates[candidate,1]
      ctx.field("SubImage.t").setValue(timepoint)
    
      #Generate contours
      
      ctx.field("itkOtsuThresholdImageFilter.update").touch()
      ctx.field("RunPythonScript.execute").touch()
      ctx.field("CSOIsoGenerator2.timePoint").setValue(timepoint)
      
      #get seedpoint contour one slice back
      csolist = ctx.field("CSOManager2.outCSOList").object()
      if timepoint>0:
        temp, previous_index = find_cso_at_time_index(csolist,slice,timepoint-1)
      elif timepoint==0:
        temp, previous_index = find_cso_at_time_index(csolist,slice,timepoint+1)
        
      first_contour_seedpoint = csolist.getCSOById(previous_index).getSeedPointsAsNumPyArray()[0]
      ctx.field("WorldVoxelConvert2.worldPos").setStringValue(np.array2string(first_contour_seedpoint)[1:-1])

      ctx.field("WorldVoxelConvert2.voxelPos").setStringValue(np.array2string(np.array([ctx.field("WorldVoxelConvert2.voxelX").value,ctx.field("WorldVoxelConvert2.voxelY").value,0]))[1:-1])
      ctx.field("CSOIsoGenerator2.startPosition").setStringValue(ctx.field("WorldVoxelConvert2.worldPos").stringValue())
      ctx.field("CSOIsoGenerator2.addCSOToGroupWithLabel").setValue(ctx.field("WorldVoxelConvert2.voxelZ").intValue())
      ctx.field("CSOIsoGenerator2.apply").touch()
      #delete old flow void CSO
      if ctx.field("CSOManager2.numCSOs").value > magnitude.shape[0]:
        temp, deleteindex= find_cso_at_time_index(csolist, slice,timepoint)
        csolist.removeCSO(deleteindex)
  ctx.field("CSOConvertToImage.apply").touch()
  
def csolistener(): 
  #Listens if a new CSO is drawn, or an existing one is edited. Will update the CSOtoimage module, so velocity values across the newly drawn/edited CSO contour are used.
  magnitude = ctx.field("ROISelect.outImage").image()
  magnitude = magnitude.getTile( (0,0,0,0,0), (magnitude.UseImageExtent,magnitude.UseImageExtent,magnitude.UseImageExtent,magnitude.UseImageExtent,magnitude.UseImageExtent) )   

  if ctx.field("CSOManager2.numCSOs").value > magnitude.shape[0]:
    csolist = ctx.field("CSOManager2.outCSOList").object()
    temp, deleteindex= find_cso_at_time_index(csolist, 0,ctx.field("View2D.timePoint").value)
    csolist.removeCSO(deleteindex)
  ctx.field("CSOConvertToImage.apply").touch()      
def find_cso_at_time_index(csolist,slice,timepoint):
  #Returns index in group, and Id of the cso.
  curslicelist=csolist.getGroupByLabel(slice)
  numcsos = curslicelist.getNumCSOs()
  for cso in range(0,numcsos):
    curslice = curslicelist.getCSOAt(cso)
    curtime = curslice.getTimePointIndex()
    
    if curtime == timepoint:
      
      return cso, curslice.getId()          
 
def save():
  mask = ctx.field("CSOConvertToImage.output0").image()
  mask = np.squeeze(mask.getTile( (0,0,0,0,0), (mask.UseImageExtent,mask.UseImageExtent,mask.UseImageExtent,mask.UseImageExtent,mask.UseImageExtent) ) ,axis=1 ) 
  
  mask = np.transpose(mask,(2,3,1,0))
  sio.savemat(ctx.field("StringUtils2.string1").value,{'mask':mask.astype('uint8')})
  print('Saving complete!')

def set_initial_candidate():
  ctx.field("XMarkerListContainer.deleteAll").touch()
  ctx.field("XMarkerListContainer.add").touch()
  ctx.field("IntervalThreshold2.threshMax").setValue(ctx.field("ImageStatistics.totalMaxVal").value*0.9)
  ctx.field("DistanceFromXMarkerList.CalcDistance").touch()
  ctx.field("XMarkerListContainer1.index").setValue(ctx.field("DistanceFromXMarkerList.Idx").value)
  ctx.field("XMarkerListContainer.posXYZ").setValue(ctx.field("XMarkerListContainer1.posXYZ").value)
  if not ctx.field("XMarkerListContainer1.numItems").value<1:
    contourcreation()