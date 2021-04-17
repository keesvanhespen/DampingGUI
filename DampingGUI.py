#!/usr/bin/env python

"""
Script to compute pulsatility values and flow velocity curves from 4D PCA images.

Written by:
Kees van Hespen, UMC Utrecht, The Netherlands
"""
from mevis import *
import numpy as np
import ntpath
import os
import scipy.ndimage 
import sklearn.linear_model
import glob
from pandas import DataFrame
import pyperclip
import shutil
import matplotlib.pyplot as plt 
import matplotlib as mpl
import time
from skimage.restoration import unwrap_phase

def init():
  global Results_and_settings
  #Initialization of the network
  #Closes all opened images
  ctx.field("currenttab").setValue(0) 
  ctx.field("ImageLoad.close").touch()
  ctx.field("ImageLoad1.close").touch()
  ctx.field("ImageLoad2.close").touch()
  ctx.field("ImageLoad5.close").touch()
  ctx.field("ImageLoad5.filename").setValue('')
  ctx.module("Load_split_dicoms").field("BoolInt.boolValue").setValue(False)
  ctx.module("Load_split_dicoms").field("BoolInt1.boolValue").setValue(False)
  ctx.module("Load_split_dicoms").field("BoolInt2.boolValue").setValue(False)
  ctx.module("Load_split_dicoms").field("ImageLoad1.filename").setValue('')
  ctx.module("Load_split_dicoms").field("ImageLoad2.filename").setValue('')
  ctx.module("Load_split_dicoms").field("ImageLoad3.filename").setValue('')
  ctx.module("Load_split_dicoms").field("ImageLoad4.filename").setValue('')
  ctx.module("Load_split_dicoms").field("ImageLoad5.filename").setValue('')
  ctx.module("Load_split_dicoms").field("ImageLoad6.filename").setValue('')
  ctx.module("Load_split_dicoms").field("ImageLoad7.filename").setValue('')
  ctx.field("FileInformation1.path").setValue('')
  ctx.module("Load_split_dicoms").field("ImageLoad1.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad3.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad4.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad5.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad6.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad2.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad7.close").touch()
  ctx.field("ImageLoad.filename").setValue('')
  ctx.field("ImageLoad1.filename").setValue('')
  ctx.field("ImageLoad2.filename").setValue('')
  ctx.module("DTF_path_and_MPR").field("BaseBypass.bypass").setValue(True)
  ctx.field("Bypass1.noBypass").setValue(True)
  #Remove generated centerline
  ctx.module("Render_Centerline").field("Vesselness.update").touch()
  ctx.module("Render_Centerline").field("DtfSkeletonization.update").touch()
  ctx.module("Render_Centerline").field("itkImageFileReader.fileName").setValue('')
  ctx.field("ImageLoad5.filename").setValue("")
  ctx.field("DrawVoxels3D1.clear").touch()
  ctx.module("Render_Centerline").field("OtsuThreshold.apply").touch()
  interface = ctx.module("Load_split_dicoms").module("PythonImage").call("getInterface")
  interface.unsetImage()
  interface = ctx.module("Load_split_dicoms").module("PythonImage1").call("getInterface")
  interface.unsetImage()
  interface = ctx.module("Load_split_dicoms").module("PythonImage2").call("getInterface")
  interface.unsetImage()
  ctx.field("Bypass1.noBypass").setValue(False)
  ctx.module("Render_Centerline").field("IntervalThreshold.threshMin").setValue(70)
  ctx.module("DTF_path_and_MPR").field("MPRPath4.fieldOfView").setValue(15)
  #sets the dicom load progressbar back to 0%
  ctx.field("prbar_loading_images").setValue(0)
  ctx.field("prbar_pulsatility").setValue(0)

  #Resets the loaded dicom trees, used in loading the dicom images
  ctx.field("SoLUTEditor.colorPoints").setStringValue('[ 0 0 0 0, 0.5 255 0 0, 0.8 1 1 1, 1 1 1 1]')
  ctx.field("SoLUTEditor.alphaPoints").setStringValue('[ 0 0, 0.5 0, 0.8 1, 1 1]')
  #Resets values for DTF skeletonization
  ctx.field("SoLUTEditor.currentIndex").setValue(1)
  ctx.module("Start_end_marker_selection").field("So3DMarkerEditor.deleteAll").touch()
  ctx.module("Render_Centerline").field("WEMIsoSurface1.startTaskSynchronous").touch()
  ctx.module("DTF_path_and_MPR").field("PathToKeyFrame.numSmoothes").setValue(5)
  #Deletes all xmarker containers
  ctx.module("DTF_path_and_MPR").field("XMarkerListContainer.deleteAll").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer.deleteAll").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer1.deleteAll").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer2.deleteAll").touch()
  ctx.module("DTF_path_and_MPR").field("XMarkerShortestPath.updateButton").touch()
  ctx.module("DTF_path_and_MPR").field("MPRPath1.upVectorColor").setValue(0, 0.66, 0)
  
  #Remove contour files
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.removeAllCSOsAndGroups").touch()
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.clear").touch()
  #Reinitialize numerical default values
  ctx.module("Magnitude_velocity_viewers").module("Colorbar").field("SoView2DMarkerEditor.textFontSize").setValue(20)

  ctx.module("FWHM_iso_value_calculation").field("Arithmetic06.arg2X").setValue(0.5)
  ctx.module("Measurements_3D_Viewer").field("Bypass_tub.bypass").setBoolValue(False)
  ctx.module("Measurements_3D_Viewer").field("Bypass_dtf.bypass").setBoolValue(True)
  
  ctx.module("Measurements_3D_Viewer").field("WEMIsoSurface.clear").touch()
  ctx.module("Select_Vessels_3D_Viewer").field("WEMIsoSurface.clear").touch()
  ctx.module("SegmentationVisualization").field("SoDrawStyle2.style").setValue('FILLED')
  ctx.module("DTF_path_and_MPR").field("MPRPath1.showUpVectors").setBoolValue(False)
  ctx.module("Render_Centerline").field("IntervalThreshold.threshMin").setValue(70)

  ctx.field("leftright").setValue('Left')
  ctx.field("Cname").setValue(0)
  #Build empty results dictionary
  Results_and_settings = {
    'Centerline Positions':np.zeros((10, 1)),
    'Centerline Vectors':np.zeros((10, 1)),
    'leftright':'left',
    'Anatomical Location':np.zeros((10, 1)),
    'Velocity Profiles':np.zeros((10, 1)),
    'Flow Profiles':np.zeros((10, 1)),
    'Arithmetic':'None',
    'Velocity pulsatility':np.zeros((10, 1)),
    'Path length':0,
    'Area':np.zeros((10, 1)),
    'Slice world matrix':[],
    'MPR world matrix':[],
    'Variable FWHM': True,
    'FWHM value': ctx.module("FWHM_iso_value_calculation").field("Arithmetic06.arg2X").value
    }

def image_corrections():
  # Does a background correction based on the static tissue mask. Used implementation of  Wouter Potters, AMC, originally written for Matlab. Performs phase unwrapping as well.
  # The corrected images will overwrite the split dicoms in the MLoutput subfolder. The original dicoms are unaltered.
  os.chdir(ctx.field("FileInformation3.dirname").value)
  ctx.module("Load_split_dicoms").field("ImageLoad4.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad5.close").touch()
  ctx.module("Load_split_dicoms").field("ImageLoad6.close").touch()
  ppca_filenames = sorted(glob.glob("*_PPCA*"), key = lambda filenm: filenm[-10])
  mpca_filenames = sorted(glob.glob("*_MPCA*"), key = lambda filenm: filenm[-10])
  mffe_filenames = sorted(glob.glob("*_MFFE*"), key = lambda filenm: filenm[-10])
  if len(ppca_filenames)!=3:#probably missing a velocity encoding direction
      print('missing a velocity encoding direction. Reload with override option for dicom splitting.')
  else:
      for scan in range(0,3):
        ctx.field("ImageLoad.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value,mffe_filenames[scan]))
        
        if len(mpca_filenames)!=3: #E.g. manual loading without MPCA images
          print('Bias field correction performed without MPCA. It is recommended to also load in MPCA images under "Manual load images".')
          mpca_available = 0
        else:
          ctx.field("ImageLoad1.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value,mpca_filenames[scan]))
          mpca_available = 1
        ctx.field("ImageLoad2.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value,ppca_filenames[scan]))
        #check if background correction has already been performed for the loaded scan
        if ctx.field("StringUtils3.boolResult").value==False:
          print('Background correction and phase unwrapping for velocity direction: %i' %scan)
          MFFE = ctx.field("ImageLoad.output0").image()
          MFFE = np.array(MFFE.getTile( (0, 0, 0, 0, 0), (MFFE.UseImageExtent, MFFE.UseImageExtent, MFFE.UseImageExtent, MFFE.UseImageExtent, MFFE.UseImageExtent) ))
          MFFE = MFFE.transpose(3,4,2,0,1)
          if mpca_available:
            MPCA = ctx.field("ImageLoad1.output0").image()
            MPCA = np.array(MPCA.getTile( (0, 0, 0, 0, 0), (MPCA.UseImageExtent, MPCA.UseImageExtent, MPCA.UseImageExtent, MPCA.UseImageExtent, MPCA.UseImageExtent) ))
          
            MPCA = MPCA.transpose(3,4,2,0,1)
          else:
            MPCA = np.zeros(MFFE.shape)
            
          PPCA = ctx.field("ImageLoad2.output0").image()
          PPCA = np.array(PPCA.getTile( (0, 0, 0, 0, 0), (PPCA.UseImageExtent,PPCA.UseImageExtent,PPCA.UseImageExtent,PPCA.UseImageExtent,PPCA.UseImageExtent) ))
        
          PPCA = PPCA.transpose(3,4,2,0,1)
          mean_ma = np.mean(MFFE, axis = 3)
          mean_cd = np.mean(MPCA, axis = 3)
          PPCA = PPCA[...,:]
          if PPCA.shape[3]>1:
            pm = np.squeeze(np.power(PPCA, 2))
          else:
            pm = np.squeeze(np.power(PPCA, 2),axis=3)
          pm = scipy.ndimage.filters.convolve(pm, np.ones((3, 3, 1, 1))/9, mode = 'constant')
          mean_pm = np.mean(pm, axis = 3)
          mean_pm = mean_pm/np.median(mean_pm)
          
          strelshape = np.array(
          [[0,0,1,0,0],
          [0,1,1,1,0],
          [1,1,1,1,1],
          [0,1,1,1,0],
          [0,0,1,0,0]])
          strelshape = strelshape[..., None, None]
          static_vessel_mask = scipy.ndimage.morphology.binary_erosion(scipy.ndimage.morphology.binary_dilation(((mean_ma > np.mean(mean_ma) + .9*np.std((mean_ma))-20)==True) & ((mean_ma < np.max(mean_ma) / 2)==True) &  ((mean_ma!=0)==True) &
              ((mean_cd <= np.median(mean_cd) + 1.5* np.std(mean_cd))==True), strelshape), strelshape)
          mean_phase_fh = np.mean(PPCA, axis = 3) 
          
          static_vessel_mask = np.squeeze(static_vessel_mask)
          interface = ctx.module("PythonImage1").call("getInterface")
          interface.setImage(static_vessel_mask.astype('int').transpose(2,0,1), minMaxValues =(np.min(static_vessel_mask),np.max(static_vessel_mask)))
          xx, yy, zz = np.meshgrid(np.linspace(0,PPCA.shape[0]-1, PPCA.shape[0]), np.linspace(0,PPCA.shape[1]-1, PPCA.shape[1]), np.linspace(0, PPCA.shape[2]-1, PPCA.shape[2]))
          k = np.random.permutation(xx[static_vessel_mask].size)
          X = np.array([xx[static_vessel_mask], yy[static_vessel_mask], zz[static_vessel_mask]])
          X = np.transpose(X)
          Y = mean_phase_fh[static_vessel_mask]
          if k.size<100000:
            numberofelements = k.size
          else:
            numberofelements = 100000
          reg = sklearn.linear_model.LinearRegression().fit(X[k[0:numberofelements], :], Y[k[0:numberofelements], :])
          correction_matrix = xx* reg.coef_[0][0] + yy * reg.coef_[0][1]+ zz*reg.coef_[0][2]+ reg.intercept_
          correction_matrix = correction_matrix[..., None]
          correction_matrix = np.tile(correction_matrix, (1, 1, 1, 1, PPCA.shape[3])).transpose(4, 0, 3, 1, 2)
          
          interface = ctx.module("PythonImage").call("getInterface")
          interface.setImage(correction_matrix, minMaxValues = (np.min(correction_matrix), np.max(correction_matrix)))
               
          #Phase unwrapping
          venc = int(ctx.field("DicomRescale5.inputIntercept").value)/np.pi
          if venc != 0:
            corrected_img = ctx.field("DicomRescale1.output0").image()
            PH = np.array(corrected_img.getTile( (0,0,0,0,0), (corrected_img.UseImageExtent,corrected_img.UseImageExtent,corrected_img.UseImageExtent,corrected_img.UseImageExtent,corrected_img.UseImageExtent) ))
            mask_arr = PH==int(ctx.field("DicomRescale5.inputIntercept").value)
            PH = np.ma.array(PH/venc,mask = mask_arr,fill_value = np.pi)
            PH_unwrapped = np.zeros(PH.shape)
            for timepoint in range(0,PH.shape[0]):
              PH_unwrapped[timepoint,0,:,:,:] = unwrap_phase(np.squeeze(PH[timepoint,...]))
            PH_unwrapped[mask_arr==True]=np.pi
            interface = ctx.module("PythonImage4").call("getInterface")  
            interface.setImage(PH_unwrapped*venc)
            
            ctx.field("ImageLoad1.close").touch()
            ctx.field("ImageLoad.close").touch()    
            ctx.field("ImageSave1.save").touch()
            tempname = ctx.field("ImageSave1.filename").value
            ctx.field("ImageLoad2.close").touch()
            shutil.move(tempname, ctx.field("ImageLoad2.filename").value)
          else:
            print('venc is zero')
      #Load background corrected scans  
      print('Loaded corrected scans')
      ctx.module("Load_split_dicoms").field("ImageLoad4.load").touch()
      ctx.module("Load_split_dicoms").field("ImageLoad5.load").touch()
      ctx.module("Load_split_dicoms").field("ImageLoad6.load").touch()
      phase_images = ctx.module("Load_split_dicoms").field("Phase_images.output0").image()
      phase_images = phase_images.getTile( (0, 0, 0, 0, 0, 0), (phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent) )   
      interface = ctx.module("Load_split_dicoms").module("PythonImage2").call("getInterface")
      interface.setImage(phase_images, minMaxValues = (np.min(phase_images), np.max(phase_images)), voxelToWorldMatrix = ctx.module("Load_split_dicoms").field("Info2.worldMatrix").value)    
      
def manual_load():
  #Loads all MFFE, MPCA, and PPCA images if manual loading is chosen.
  if not os.path.exists(ctx.field("FileInformation3.dirname").value):  
    os.mkdir(ctx.field("FileInformation3.dirname").value); 
  if not ctx.field("FileInformation3.exists").value == True or ctx.field("override_splitdicoms").boolValue() == True:
    head = ctx.field("FileInformation3.dirname").value
    temp, tail = ntpath.split(ctx.field("manual_mffe1").value)
    shutil.copyfile(ctx.field("manual_mffe1").value,head+tail.replace('.dcm', '')+'_0_MFFE.dcm')
    temp, tail = ntpath.split(ctx.field("manual_mffe2").value)
    shutil.copyfile(ctx.field("manual_mffe2").value,head+tail.replace('.dcm', '')+'_1_MFFE.dcm')
    temp, tail = ntpath.split(ctx.field("manual_mffe3").value)
    shutil.copyfile(ctx.field("manual_mffe3").value,head+tail.replace('.dcm', '')+'_2_MFFE.dcm')
    if ctx.field("manual_mpca1").value:
      temp, tail = ntpath.split(ctx.field("manual_mpca1").value)
      shutil.copyfile(ctx.field("manual_mpca1").value,head+tail.replace('.dcm', '')+'_0_MPCA.dcm')
      temp, tail = ntpath.split(ctx.field("manual_mpca2").value)
      shutil.copyfile(ctx.field("manual_mpca2").value,head+tail.replace('.dcm', '')+'_1_MPCA.dcm')
      temp, tail = ntpath.split(ctx.field("manual_mpca3").value)
      shutil.copyfile(ctx.field("manual_mpca3").value,head+tail.replace('.dcm', '')+'_2_MPCA.dcm')
    temp, tail = ntpath.split(ctx.field("manual_ppca1").value)
    shutil.copyfile(ctx.field("manual_ppca1").value,head+tail.replace('.dcm', '')+'_0_PPCA.dcm')
    temp, tail = ntpath.split(ctx.field("manual_ppca2").value)
    shutil.copyfile(ctx.field("manual_ppca2").value,head+tail.replace('.dcm', '')+'_1_PPCA.dcm')
    temp, tail = ntpath.split(ctx.field("manual_ppca3").value)
    shutil.copyfile(ctx.field("manual_ppca3").value,head+tail.replace('.dcm', '')+'_2_PPCA.dcm')
    
  else:
    print('Magnitude, modulus, and phase images already split. For override, check the options.')
  load_images()
  
def check_dropped_image():
  # Splits dicoms based on MFFE/MPCA/PPCA imagetype. Also checks for existance of T1W image.
  # Also checks for the existance of the MLOutput folder, and will create one if non-existant.
  ctx.field("prbar_loading_images").setValue(0)
  if "PCA" in ctx.field("DicomFrameTagInfo4.tagValue").value:
    #Creates MLoutput folder if it does not exist
    if not os.path.exists(ctx.field("FileInformation3.dirname").value):
      os.mkdir(ctx.field("FileInformation3.dirname").value)
    # if split dicoms do not exist in the MLoutput folder
    if not ctx.field("FileInformation3.exists").value == True or ctx.field("override_splitdicoms").boolValue() == True:
      print('Loading dicoms...')
      ctx.field("DirectDicomImport.dplImport").touch()
      ctx.field("DirectDicomImport2.dplImport").touch()
      ctx.field("prbar_loading_images").setValue(0.25)
      for volume_nr in range(0, ctx.field("DirectDicomImport.numVolumes").value):
        
        ctx.field("DirectDicomImport.outVolume").setValue(volume_nr)        
        ctx.field("ImageSave.save").touch()

        ctx.field("prbar_loading_images").setValue(0.25 + 0.5 * (volume_nr + 1) / ctx.field("DirectDicomImport.numVolumes").value)
        
      print('Dicoms split into magnitude, modulus and phase images')
      #Checks if a T1W M_FFE is present in the folder. If so, its  
      ctx.field("DirectDicomImport1.dplImport").touch()
      if ctx.field("DirectDicomImport1.numVolumes").value > 0:
        for volume_nr in range(0,ctx.field("DirectDicomImport1.numVolumes").value):
          ctx.field("DirectDicomImport1.outVolume").setValue(volume_nr)
          if 'M_FFE' in ctx.field("DicomFrameTagInfo6.sharedTagValue").stringValue():
            ctx.field("ImageSave2.save").touch()
      else:
        print('No T1w image found in this folder')  
    else:
      print('Magnitude, modulus, and phase images already split. For override, check the options.')
    load_images()   
  #if dropped image does not contains "PCA" in the file description
  else:
    print('Image is not a PCA image')
  ctx.field("ImageLoad5.close").touch()
  ctx.field("ImageLoad4.close").touch() 
def load_images():
  # Split dicom image loading:
  # change current dir to MLoutput subfolder, and find all PPCA images there.
  os.chdir(ctx.field("FileInformation3.dirname").value)
  ctx.field("fileName").setValue(ctx.field("FileInformation1.dirname").value)
  ctx.field("DirectDicomImport2.dplImport").touch()
  ppca_filenames = glob.glob("*_PPCA*")
  mffe_filenames = glob.glob("*_MFFE*")
  mpca_filenames = glob.glob("*_MPCA*")
  t1w_filenames = glob.glob("*T1w*")
  
  if len(ppca_filenames)!=3:#probably missing a velocity encoding direction
    print('probably missing a velocity encoding direction. Reload with override option for dicom splitting.')
  else:
    # load mffe images
    ctx.module("Load_split_dicoms").field("ImageLoad3.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value, mffe_filenames[0]))
    ctx.module("Load_split_dicoms").field("ImageLoad2.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value, mffe_filenames[1]))
    ctx.module("Load_split_dicoms").field("ImageLoad7.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value, mffe_filenames[2]))
    magnitude = ctx.module("Load_split_dicoms").field("PythonArithmetic.output0").image()
    magnitude = magnitude.getTile( (0,0,0,0,0,0), (magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent) )   
    
    interface = ctx.module("Load_split_dicoms").module("PythonImage").call("getInterface")
    interface.setImage(magnitude, minMaxValues = (np.min(magnitude), np.max(magnitude)), voxelToWorldMatrix = ctx.module("Load_split_dicoms").field("Info.worldMatrix").value)    
    # load ppca images
    for scan in range(0,3):

      found_enc_dir = 0
      #AP direction
      if ppca_filenames[scan].find('_1_PPCA')>0:
        ctx.module("Load_split_dicoms").field("ImageLoad6.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value, ppca_filenames[scan]))
        found_enc_dir = 1
      #RL direction  
      if ppca_filenames[scan].find('_0_PPCA')>0:
        ctx.module("Load_split_dicoms").field("ImageLoad5.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value, ppca_filenames[scan]))
        found_enc_dir = 1
      #FH direction  
      if ppca_filenames[scan].find('_2_PPCA')>0:
        ctx.module("Load_split_dicoms").field("ImageLoad4.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value, ppca_filenames[scan]))
        found_enc_dir = 1
      if found_enc_dir != 1:
        print('non orthogonal/no velocity encoding directions found in header')

    # set phase images output
    phase_images = ctx.module("Load_split_dicoms").field("Phase_images.output0").image()
    phase_images = phase_images.getTile( (0, 0, 0, 0, 0, 0), (phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent) )   
    interface = ctx.module("Load_split_dicoms").module("PythonImage2").call("getInterface")
    interface.setImage(phase_images, minMaxValues = (np.min(phase_images),np.max(phase_images)), voxelToWorldMatrix = ctx.module("Load_split_dicoms").field("Info2.worldMatrix").value)    
    
    # set view orientation/world matrix PCA
    ds = ctx.field("DirectDicomImport2.output0").getDicomTree()
    vieworientation = ds.getTag("(2001,105f)").getSequenceItem(0).getTag("(2005,1081)").value()
    Tmic = np.array([[ 0, -1, 0, 0],[ -1, 0, 0, 0],[ 0, 0 ,1, 0],[ 0, 0, 0, 1]])
    if vieworientation == 'AP':
      Tsom = np.array([[0, -1, 0, 0],[ 0, 0, 1, 0],[ 1, 0, 0, 0],[ 0, 0, 0, 1]])
    elif vieworientation == 'RL':
      Tsom = np.array([[0, 0, -1, 0],[ 0, -1, 0, 0],[ 1, 0, 0, 0],[ 0, 0, 0, 1]])
    elif vieworientation == 'FH':
      Tsom = np.array([[0, -1, 0, 0],[ -1, 0, 0, 0],[ 0, 0, 1, 0],[ 0, 0, 0, 1]])
    else:
      print('View angle not found in dicom header. Assuming transversal (FH) view angle')
      Tsom = np.array([[0 ,-1 ,0 ,0],[ -1, 0, 0, 0],[ 0, 0, 1, 0],[ 0 ,0, 0 ,1]])
    Tform = np.transpose( Tsom.dot(Tmic) )
    ctx.module("Phase_to_velocity").field("MatrixArithmetic2.matrixB").setValue(Tform)
    original_image = ctx.module("Phase_to_velocity").field("Bypass.output0").image()
    original_image = original_image.getTile((0, 0, 0, 0, 0,0), (original_image.UseImageExtent, original_image.UseImageExtent, original_image.UseImageExtent, original_image.UseImageExtent, original_image.UseImageExtent,original_image.UseImageExtent) )
    worldvoxelmat = np.array(ctx.module("Phase_to_velocity").field("MatrixArithmetic2.outputMatrixC").value)
    original_image_inworld = np.transpose(worldvoxelmat[0:3,0:3].dot(np.transpose(original_image,(5,4,3,2,0,1))),(0,5,4,3,2,1))
    interface = ctx.module("Phase_to_velocity").module("PythonImage").call("getInterface")
    interface.setImage(original_image_inworld, minMaxValues = (np.min(original_image_inworld), np.max(original_image_inworld))) 
  # load T1w
  if len(t1w_filenames) > 0:
    ctx.module("Load_split_dicoms").field("ImageLoad1.filename").setValue(os.path.join(ctx.field("FileInformation3.dirname").value, t1w_filenames[0]))
    magnitude = ctx.module("Load_split_dicoms").field("OrthoReformat3.output1").image()
    magnitude = magnitude.getTile( (0, 0, 0, 0, 0, 0), (magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent) )   
    interface = ctx.module("Load_split_dicoms").module("PythonImage1").call("getInterface")
    interface.setImage(magnitude, minMaxValues = (np.min(magnitude), np.max(magnitude)), voxelToWorldMatrix = ctx.module("Load_split_dicoms").field("Info1.worldMatrix").value)    
  ctx.field("ImageLoad5.close").touch()
  ctx.field("ImageLoad4.close").touch()   
  ctx.field("prbar_loading_images").setValue(0.9)
  ctx.module("Render_Centerline").field("Vesselness.update").touch()
  ctx.module("Render_Centerline").field("DtfSkeletonization.update").touch()
  ctx.field("prbar_loading_images").setValue(1)
  get_LUT_set_LUT()  
  
def delete_all_markers():
  #deletes all drawn markers in the second tab of the GUI.
  global backgrounds
  ctx.module("DTF_path_and_MPR").field("BaseBypass.bypass").setValue(True)
  ctx.module("DTF_path_and_MPR").field("XMarkerListContainer.deleteAll").touch()
  ctx.module("DTF_path_and_MPR").field("XMarkerShortestPath.updateButton").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer1.deleteAll").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer.deleteAll").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer2.deleteAll").touch()
  ctx.module("Render_Centerline").field("MaskToMarkers.update").touch()
  ctx.module("Measurements_3D_Viewer").field("WEMIsoSurface.clear").touch()
  ctx.module("Select_Vessels_3D_Viewer").field("WEMIsoSurface.clear").touch()
  ctx.field("prbar_pulsatility").setValue(0)

  if 'backgrounds' in globals():
    subplot.canvas.restore_region(backgrounds)    
    subplot.canvas.draw()
  

def get_LUT_set_LUT():
  #sets the LUT values of the 3D rendering of the vessels based on thresholded magnitude image, and the standard deviation of the thresholded image values.
  colorpoints = (ctx.field("SoLUTEditor.colorPoints").stringValue()).replace(',', '')
  colorpointslist = colorpoints[1:-1].split()
  floatarray_color = np.array([float(i) for i in colorpointslist])
  reshapedfloatarray_color = np.reshape(floatarray_color,(-1,4))
  reshapedfloatarray_color[1, 0] = ctx.field("Arithmetic02.resultX").intValue()
  reshapedfloatarray_color[3, 0] = ctx.field("Arithmetic02.arg1X").intValue() + 5 * ctx.field("Arithmetic02.arg2X").intValue()
  reshapedfloatarray_color[2, 0] = (reshapedfloatarray_color[3, 0] - reshapedfloatarray_color[1, 0])/2+reshapedfloatarray_color[1, 0]
  stringarray_color = np.array2string(reshapedfloatarray_color)
  stringarray_color = stringarray_color.replace(']\n [', ', ').replace('[[', '[').replace(']]', ']')
  ctx.field("SoLUTEditor.colorPoints").setStringValue(stringarray_color)
  
  alphapoints= (ctx.field("SoLUTEditor.alphaPoints").stringValue()).replace(',', '')
  alphapointslist = alphapoints[1:-1].split()
  floatarray_alpha = np.array([float(i) for i in alphapointslist])
  reshapedfloatarray_alpha = np.reshape(floatarray_alpha,(-1,2))
  reshapedfloatarray_alpha[1, 0] = ctx.field("Arithmetic02.resultX").intValue()
  reshapedfloatarray_alpha[3, 0] = ctx.field("Arithmetic02.arg1X").intValue() + 5 * ctx.field("Arithmetic02.arg2X").intValue()
  reshapedfloatarray_alpha[2, 0] = (reshapedfloatarray_alpha[3, 0] - reshapedfloatarray_alpha[1, 0])/2+reshapedfloatarray_alpha[1, 0]
  stringarray_alpha = np.array2string(reshapedfloatarray_alpha)
  stringarray_alpha = stringarray_alpha.replace(']\n [', ',').replace('[[', '[').replace(']]', ']')
  ctx.field("SoLUTEditor.alphaPoints").setStringValue(stringarray_alpha)

 
def contourcreation():
  # Creates iso contours at the FWHM value, and initializes the graphs on which the velocity curves are drawn.
  global Results_and_settings, Iso_vals  
  if ctx.field("FieldListener.sourceFieldValue").value=="Valid":
    ctx.module("MPRTubularTracking").field("XMarkerListContainer3.deleteAll").touch()
    ctx.module("MPRTubularTracking").field("XMarkerListContainer3.add").touch()
    ctx.module("MPRTubularTracking").field("WorldVoxelConvert1.intVoxelCoords").setValue(False)
    ctx.module("MPRTubularTracking").field("TubularTracking.update").touch()
    ctx.module("MPRTubularTracking").field("SoWEMConvertInventor.apply").touch()
    if ctx.module("MPRTubularTracking").field("TrackedPoints.numItems").value>0:
      T1_available=1
    else:
      T1_available=0  
  else:
    T1_available=0
  # Only possible if two 3D markers are clicked in the 3D viewer.    
  if ctx.module("Start_end_marker_selection").field("So3DMarkerEditor.numItems").value == 2:
    ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.removeAllCSOsAndGroups").touch()
    #Take maximum of magnitude image across all time points
    MPR_reconstruction_mag = ctx.field("DTF_path_and_MPR.output0").image()
    labellist=[]

    MPR_reconstruction_mag = MPR_reconstruction_mag.getTile( (0,0,0,0,0), (MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent) )   
    cso_areas = np.zeros((MPR_reconstruction_mag.shape[0], MPR_reconstruction_mag.shape[2]))
    tubetracking_areas = np.zeros((1,MPR_reconstruction_mag.shape[2]))
    Iso_vals = np.zeros((MPR_reconstruction_mag.shape[0], MPR_reconstruction_mag.shape[2]))
    No_of_CSOs = 0
    for timepoint in range(0, MPR_reconstruction_mag.shape[0]):
      cur_magnitude = MPR_reconstruction_mag[timepoint*ctx.field("Use_variable_FWHM").value, ...]
      interface = ctx.module("PythonImage2").call("getInterface")
      interface.setImage(cur_magnitude, minMaxValues =(np.min(cur_magnitude), np.max(cur_magnitude)), voxelToWorldMatrix = ctx.module("FWHM_iso_value_calculation").field("Info2.worldMatrix").value )    

      #
      ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.timePoint").setValue(timepoint)
      if ctx.module("FWHM_iso_value_calculation").field("Arithmetic010.resultX").value == 0: #if the window size is larger than the number of slices
        ctx.module("FWHM_iso_value_calculation").field("itkOtsuThresholdImageFilter.update").touch()
        ctx.module("FWHM_iso_value_calculation").field("RunPythonScript.execute").touch()
      
      for slice in range(0, MPR_reconstruction_mag.shape[2]):  
        if T1_available:
           
          ctx.module("MPRTubularTracking").field("WorldVoxelConvert4.voxelZ").setValue(int(slice)+0.5)
          ctx.module("MPRTubularTracking").field("WEMClipPlaneToCSO.apply").touch()
          ctx.module("MPRTubularTracking").field("CSOInfo.apply").touch()

          MLAB.processEvents()
          MLAB.processInventorQueue()
          tubetracking_areas[0,slice]= ctx.module("MPRTubularTracking").field("CSOInfo.csoArea").value
          csolist = ctx.module("MPRTubularTracking").field("WEMClipPlaneToCSO.outCSOList").object()
          
          try:
            cur_cso = csolist.getCSOAt(0)
            cur_cso.applyTransformationMatrix(np.array(ctx.module("MPRTubularTracking").field("MatrixArithmetic.outputMatrixC").value))
            seeds = cur_cso.getSeedPointsAsNumPyArray()
            y_val = seeds[0,0]
            x_val = seeds[0,1]
          except:
            y_val = np.round(MPR_reconstruction_mag.shape[4]/2).astype(int)+3
            x_val = np.round(MPR_reconstruction_mag.shape[4]/2).astype(int)+3
        else:
          y_val = np.round(MPR_reconstruction_mag.shape[4]/2).astype(int)+3
          x_val = np.round(MPR_reconstruction_mag.shape[4]/2).astype(int)+3
        max_iso_val = np.max(MPR_reconstruction_mag[timepoint, 0, slice,:,:])-10 if np.max(MPR_reconstruction_mag[timepoint, 0, slice,:,:])-10>0 else 0
        ctx.module("FWHM_iso_value_calculation").field("Arithmetic08.arg2X").setValue(max_iso_val) 
        
        if not MPR_reconstruction_mag[0, 0, slice,x_val.astype(int),y_val.astype(int)] <1:
          if ctx.module("FWHM_iso_value_calculation").field("Arithmetic010.resultX").value == 1: #if the window size is smaller than the number of slices
            ctx.module("FWHM_iso_value_calculation").field("MirrorSubimage.Center").setValue(slice)
            ctx.module("FWHM_iso_value_calculation").field("itkOtsuThresholdImageFilter.update").touch()
            ctx.module("FWHM_iso_value_calculation").field("RunPythonScript.execute").touch()
          
          if not ctx.field("Use_variable_FWHM").value:
            ctx.module("FWHM_iso_value_calculation").field("Arithmetic08.arg2X").setValue(max_iso_val) 
            
          ctx.field("WorldVoxelConvert2.voxelPos").setStringValue(np.array2string(np.array([y_val,x_val, slice]))[1:-1])      
          ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.startPosition").setStringValue(ctx.field("WorldVoxelConvert2.worldPos").stringValue())
          ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.addCSOToGroupWithLabel").setValue(ctx.field("WorldVoxelConvert2.voxelZ").intValue())
          
          try:
            ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.apply").touch()

            #check contour area
            csolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
            if csolist.getNumCSOs()> No_of_CSOs:
              No_of_CSOs = csolist.getNumCSOs()
              ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOInfo.csoShowByIndexOrId").setValue(csolist.getCSOIdList()[-1])
  
              cso_areas[timepoint, slice] = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOInfo.csoArea").value
              labellist.append(ctx.field("WorldVoxelConvert2.voxelZ").intValue())
              Iso_vals[timepoint, slice] = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").value
          except:
              print('not ok')
          
      ctx.field("prbar_pulsatility").setValue(0.8*(timepoint+1)/MPR_reconstruction_mag.shape[0])
    if T1_available:
      d = np.nan_to_num(cso_areas - np.tile(tubetracking_areas,[MPR_reconstruction_mag.shape[0],1]))/np.tile(tubetracking_areas,[MPR_reconstruction_mag.shape[0],1])
      outliers = np.abs(d) > 0.3
    else:
      #Check if there are outliers with low areas. These could be incorrect contours in flow voids.
      if MPR_reconstruction_mag.shape[2] >= MPR_reconstruction_mag.shape[0]:
        d = np.nan_to_num(cso_areas -  running_median(cso_areas,3,1))
        #d = np.abs(cso_areas - np.expand_dims(np.median(cso_areas[:, ~np.all(cso_areas == 0, axis = 0)], axis = 1), 1))
        axe = 0
      elif MPR_reconstruction_mag.shape[2] < MPR_reconstruction_mag.shape[0]:
        d = np.nan_to_num(cso_areas -  running_median(cso_areas,3,0))
        axe = 1
      mdev = np.nanmedian(np.abs(d)[np.nan_to_num(cso_areas)!=0])
      s = np.abs(d) / mdev if mdev else 0.
      outliers = s < 6
      np.all(outliers, axis = axe)
      np.argwhere(outliers == 0)
      
    flow_void_candidates = np.argwhere(outliers == 0)

    done = scipy.ndimage.uniform_filter(np.abs(d), size=(3,0), origin=-1, mode='wrap')
    
    if np.squeeze(np.argwhere(np.all(cso_areas == 0,axis = 0))).size:
      flow_void_candidates = flow_void_candidates[np.invert(np.isin(flow_void_candidates[:, 1], np.squeeze(np.argwhere(np.all(cso_areas == 0,axis = 0))))), :]
    
    unique_labels = np.sort(np.unique(np.array(labellist)))
    
    print('Flow void candidates to be fixed [timepoint,slice]:');print(flow_void_candidates)
    if flow_void_candidates.size!= 0 and ctx.module("Start_end_marker_selection").field("BoolInt.boolValue").value:
      outlierhaslarger = (d[flow_void_candidates[:,0], flow_void_candidates[:,1]])/np.abs(d[flow_void_candidates[:,0], flow_void_candidates[:,1]])
      #try to fix flow void candidates
      for candidate in range(0, flow_void_candidates.shape[0]):
        timepoint, slice = flow_void_candidates[candidate, 0],flow_void_candidates[candidate, 1]
        cur_magnitude = MPR_reconstruction_mag[timepoint*ctx.field("Use_variable_FWHM").value, ...]
        interface = ctx.module("PythonImage2").call("getInterface")
        interface.setImage(cur_magnitude, minMaxValues = (np.min(cur_magnitude), np.max(cur_magnitude)), voxelToWorldMatrix = ctx.module("FWHM_iso_value_calculation").field("Info2.worldMatrix").value)    
        trytime=1
        isoutlier=1
        #Generate contours 
        while isoutlier and trytime < 30:
          if trytime==1:
            ctx.module("FWHM_iso_value_calculation").field("MirrorSubimage.Center").setValue(slice)
            ctx.module("FWHM_iso_value_calculation").field("itkOtsuThresholdImageFilter.update").touch()
            ctx.module("FWHM_iso_value_calculation").field("RunPythonScript.execute").touch()
            ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.timePoint").setValue(timepoint)
          
            slice_index = np.where(unique_labels==slice)[0][0]
            #get seedpoint contour one slice back
            csolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
    
            if slice_index > 0:
              temp, previous_index = find_cso_at_time_index(csolist, unique_labels[slice_index-1], timepoint)
            elif timepoint>0:
              temp, previous_index = find_cso_at_time_index(csolist, unique_labels[slice_index], timepoint - 1)
            else:
              temp, previous_index = find_cso_at_time_index(csolist, unique_labels[slice_index+1], timepoint)
              
            first_contour_seedpoint = csolist.getCSOById(previous_index).getSeedPointsAsNumPyArray()[0]
            ctx.field("WorldVoxelConvert2.worldPos").setStringValue(np.array2string(first_contour_seedpoint)[1:-1])
            ctx.field("WorldVoxelConvert2.voxelPos").setStringValue(np.array2string(np.array([ctx.field("WorldVoxelConvert2.voxelX").value, ctx.field("WorldVoxelConvert2.voxelY").value, slice]))[1:-1])
            ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.startPosition").setStringValue(ctx.field("WorldVoxelConvert2.worldPos").stringValue())
            ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.addCSOToGroupWithLabel").setValue(ctx.field("WorldVoxelConvert2.voxelZ").intValue())
          
          update_val = (30-trytime)/6*np.abs(d)[timepoint,slice]
          update_val = update_val if np.abs(update_val)<10 else 10
          iso_val = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").value+outlierhaslarger[candidate]*update_val
          max_val = ctx.module("FWHM_iso_value_calculation").field("Arithmetic08.arg2X").value
          iso_val = iso_val if iso_val > 0 else ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").value+0.5*outlierhaslarger[candidate]*(30-trytime)/3*np.abs(d)[timepoint,slice]
          if iso_val > max_val:
            iso_val = max_val 
          ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").setValue(iso_val)
          ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.apply").touch()
          #delete old flow void CSO
          if ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.numCSOs").value > ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.numGroups").value:
            temp, deleteindex = find_cso_at_time_index(csolist, slice, timepoint)
            csolist.removeCSO(deleteindex)
          #check contour area
          csolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
          ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOInfo.csoShowByIndexOrId").setValue(csolist.getCSOIdList()[-1])
          cso_areas[timepoint, slice] = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOInfo.csoArea").value
          Iso_vals[timepoint, slice] = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").value
          if MPR_reconstruction_mag.shape[2] >= MPR_reconstruction_mag.shape[0]:
            d = np.nan_to_num(cso_areas - running_median(cso_areas,3,1))
            axe=0
          elif MPR_reconstruction_mag.shape[2] < MPR_reconstruction_mag.shape[0]:
            d = np.nan_to_num(cso_areas - running_median(cso_areas,3,0))
            axe=1
          mdev = np.median(np.abs(d)[cso_areas!=0])
          s = np.abs(d) / mdev if mdev else 0.
          outliers = s < 6
          isoutlier = not ((np.abs(d[timepoint,slice]) / mdev if mdev else 0.)<6 )
          outlierhaslarger[candidate] = (d[timepoint,slice])/np.abs(d[timepoint,slice])
          trytime = trytime+1
      np.all(outliers, axis = axe)
      np.argwhere(outliers == 0)
      flow_void_candidates = np.argwhere(outliers == 0) 
      print('Flow void candidates left [timepoint,slice]:');print(flow_void_candidates[flow_void_candidates[:,1].argsort()])

    Results_and_settings = {
     'Centerline Positions':np.zeros((MPR_reconstruction_mag.shape[2], 1)),
     'Centerline Vectors':np.zeros((MPR_reconstruction_mag.shape[2], 1)),
     'leftright':'left',
     'Anatomical Location':np.zeros((MPR_reconstruction_mag.shape[2], 1)),
     'Velocity Profiles':np.zeros((MPR_reconstruction_mag.shape[2], 1)),
     'Flow Profiles':np.zeros((MPR_reconstruction_mag.shape[2], 1)),
     'Arithmetic':'None',
     'Velocity pulsatility':np.zeros((MPR_reconstruction_mag.shape[2], 1)),
     'Path length':0,
     'Area':np.zeros((MPR_reconstruction_mag.shape[2], 1)),
     'Slice world matrix':[],
     'MPR world matrix':[],
     'Variable FWHM': ctx.field("Use_variable_FWHM").value,
     'FWHM value': Iso_vals
    }
    
    #store the world matrices of all the cso slices 
    for slice in range(0, MPR_reconstruction_mag.shape[2]):      
      ctx.module("DTF_path_and_MPR").field("MPRPath1.currentKeyFrame").setValue(slice)
      Results_and_settings['Slice world matrix'].append(np.array(ctx.module("DTF_path_and_MPR").field("Info.worldMatrix").value))
    Results_and_settings['MPR world matrix'] = np.array(ctx.module("DTF_path_and_MPR").field("Info1.worldMatrix").value)
    ctx.field("prbar_pulsatility").setValue(0.9)        
    phase_to_velocity()
    ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.apply").touch()
    ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage1.apply").touch()
    plot_clear()
    plot_all_slice_velocities()
    
    ctx.field("prbar_pulsatility").setValue(1)   
    cso_wem_computation()

def phase_to_velocity():
 
  #Takes the normalized vector components of the MPR slices, and calculates the velocity perpendicular to the surface of the MPR planes, using the individual phase components and their respective contributions.
  MPR_FH = ctx.module("Phase_to_velocity").field("MPR_FH.output0").image()
  MPR_FH = MPR_FH.getTile( (0, 0, 0, 0, 0), (MPR_FH.UseImageExtent, MPR_FH.UseImageExtent, MPR_FH.UseImageExtent, MPR_FH.UseImageExtent, MPR_FH.UseImageExtent) )
  MPR_RL = ctx.module("Phase_to_velocity").field("MPR_RL.output0").image()
  MPR_RL = MPR_RL.getTile( (0, 0, 0, 0, 0), (MPR_RL.UseImageExtent, MPR_RL.UseImageExtent, MPR_RL.UseImageExtent, MPR_RL.UseImageExtent, MPR_RL.UseImageExtent) )
  MPR_AP = ctx.module("Phase_to_velocity").field("MPR_AP.output0").image()
  MPR_AP = MPR_AP.getTile((0, 0, 0, 0, 0), (MPR_AP.UseImageExtent, MPR_AP.UseImageExtent, MPR_AP.UseImageExtent, MPR_AP.UseImageExtent, MPR_AP.UseImageExtent) )
  mpl.rcParams['toolbar'] = 'None'
  MPR_velocity = np.zeros(MPR_AP.shape)

  MPR = ctx.module("Phase_to_velocity").field("MPRPath3.output1").image()
  MPR = MPR.getTile((0, 0, 0, 0, 0, 0), (MPR.UseImageExtent, MPR.UseImageExtent, MPR.UseImageExtent, MPR.UseImageExtent, MPR.UseImageExtent,MPR.UseImageExtent) )

  
  for slice in range(0,ctx.module("Phase_to_velocity").field("MPRPath3.maxKeyFrame").value + 1):
    MPR_slice = MPR[:,:,:,slice,...]
    ctx.module("Phase_to_velocity").field("MPRPath3.currentKeyFrame").setValue(slice)
    plane_normal = -np.array(ctx.module("Phase_to_velocity").field("DecomposeVector41.v").value)
    norm_vect_non_unit = np.array([plane_normal[0],plane_normal[1],plane_normal[2]])
    norm_vect = norm_vect_non_unit/np.sqrt(np.sum(np.power(plane_normal[0:3],2)))
    direction_angle = np.arccos(np.transpose(np.transpose(MPR_slice).dot(norm_vect_non_unit))/(np.sqrt(np.sum(np.square(MPR_slice), axis=0))*np.sqrt(np.sum(np.power(plane_normal[0:3],2)))))
    velocity_direction = np.where((direction_angle)<(np.pi/2),1,-1)
    MPR_velocity[:,:,slice,...]= velocity_direction*np.abs(norm_vect[2] * MPR_FH[:, :, slice, ...] + norm_vect[0] * MPR_RL[:, :, slice, ...] + norm_vect[1] * MPR_AP[:, :, slice, ...])
  interface = ctx.module("Phase_to_velocity").module("MPR_velocity").call("getInterface")
  interface.setImage(MPR_velocity, minMaxValues = (np.min(MPR_velocity), np.max(MPR_velocity)))

def select_CSO():
  #Sets the color of the clicked CSO. Clicking it to red will omit it from further calculations
  if ctx.field("currenttab").value == 2:
    if ctx.field("CSOManager2.csoSinglePathPointColor").value[1] > 0 and ctx.field("SoView2DCSOExtensibleEditor.isEditingExistingCSO").value == True: #Case if CSO is green
      ctx.field("CSOManager2.csoSinglePathPointColor").setValue(0.66, 0, 0)
    elif ctx.field("SoView2DCSOExtensibleEditor.isEditingExistingCSO").value == True:# If CSO is red
      ctx.field("CSOManager2.csoSinglePathPointColor").setValue(0, 0.66, 0)

def plot_clear():
  #Clears the plot
  global subplot, ax, backgrounds, subplot_allslices, ax_allslices, backgrounds_allslices
  control = ctx.control("canvas").object()
  control.figure().clear()
  # clear the user interaction of plotC example:
  control.myUserInteraction = None
  subplot = ctx.control("canvas").object().figure()
  ax = subplot.add_subplot(111)
  ax.set_xlabel('Cardiac cycle') 
  ax.set_ylabel('Max velocity (cm/s)') 
  ax.set_title('Velocity') 
  
  # We need to draw the canvas before we start animating...
  subplot.canvas.draw()
  backgrounds = subplot.canvas.copy_from_bbox(subplot.bbox)
  plot_velocity_curve()
  
  control = ctx.control("canvas_allslices").object()
  control.figure().clear()
  # clear the user interaction of plotC example:
  control.myUserInteraction = None
  
  subplot_allslices = ctx.control("canvas_allslices").object().figure()
  ax_allslices = subplot_allslices.add_subplot(111)
  ax_allslices.set_xlabel('MPR slice') 
  ax_allslices.set_ylabel('Velocity (cm/s)') 
  ax_allslices.set_title('Min-mean-max velocities') 
  
  # We need to draw the canvas before we start animating...
  subplot_allslices.canvas.draw()
  backgrounds_allslices = subplot_allslices.canvas.copy_from_bbox(subplot_allslices.bbox)
    
def plot_velocity_curve():
  #Plots the velocity curve for the currently selected slice in the MPR.
  global subplot, ax, backgrounds
  if ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.output0").isValid() and ctx.field("currenttab").value == 2 and ctx.field("currenttabplots").value == 0:
    velocity_ROI = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.output0").image()
    velocity_ROI = velocity_ROI.getTile( (0,0,0,0,0), (velocity_ROI.UseImageExtent,velocity_ROI.UseImageExtent,velocity_ROI.UseImageExtent,velocity_ROI.UseImageExtent,velocity_ROI.UseImageExtent) )
    velocity_ROI = velocity_ROI[:, :, ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").intValue(), ...]
    flattened_v_ROI = velocity_ROI.reshape(velocity_ROI.shape[0], -1)
    if ctx.field("velocity_component").value == 2:
      areas_velocity = np.squeeze(np.sum(flattened_v_ROI, axis = -1) / np.count_nonzero(flattened_v_ROI, axis = -1))
    elif ctx.field("velocity_component").value == 1:
      areas_velocity = np.squeeze(np.max(flattened_v_ROI, axis = -1))

    control = ctx.control("canvas").object()
    control.myUserInteraction = None
    
    x = np.linspace(0,areas_velocity.size-1, areas_velocity.size)
    style='r-'
    def plot(ax, style):
      return ax.plot(x, areas_velocity, linestyle = '-', color = style, animated = True)[0]
    line = plot(ax, 'r')

    subplot.canvas.restore_region(backgrounds)    
    line.set_ydata(areas_velocity)
    ax.set_ylim([0, np.max(areas_velocity) + 10])
    ax.set_xlim([0, areas_velocity.size - 1])
    subplot.canvas.draw()
    ax.draw_artist(line)
    
    subplot.canvas.blit(ax.bbox)

def plot_all_slice_velocities():
  # plots mean, max and minimum velocity for all slices. 
  global subplot_allslices,ax_allslices, backgrounds_allslices
  if ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.output0").isValid() and ctx.field("currenttab").value == 2 and ctx.field("currenttabplots").value == 1:
    velocity_ROI = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.output0").image()
    velocity_ROI = velocity_ROI.getTile( (0, 0, 0, 0, 0), (velocity_ROI.UseImageExtent, velocity_ROI.UseImageExtent, velocity_ROI.UseImageExtent, velocity_ROI.UseImageExtent, velocity_ROI.UseImageExtent) )
    means = np.zeros((velocity_ROI.shape[2], 1))
    maxes = np.zeros((velocity_ROI.shape[2], 1))
    mins = np.zeros((velocity_ROI.shape[2], 1))
    for slice in range(0, velocity_ROI.shape[2]):
      velocity_ROI_selection = velocity_ROI[:, :, slice, ...]
      flattened_v_ROI = velocity_ROI_selection.reshape(velocity_ROI_selection.shape[0], -1)
      if ctx.field("velocity_component").value == 2:
        areas = np.squeeze(np.sum(flattened_v_ROI, axis = -1) / np.count_nonzero(flattened_v_ROI, axis = -1))
      elif ctx.field("velocity_component").value == 1:
        areas = np.squeeze(np.max(flattened_v_ROI,axis = -1))
      
      means[slice] = np.mean(areas)
      maxes[slice] = np.max(areas)
      mins[slice] = np.min(areas)
    control = ctx.control("canvas").object()
    control.myUserInteraction = None

    x = np.linspace(0, maxes.size-1, maxes.size)
    style = 'r'
    def plot(ax, style):
      return ax_allslices.plot(x, maxes, linestyle = '-', color = style, animated = True)[0]
    line_max = plot(ax, 'r')
    line_mean = plot(ax, 'b')
    line_min = plot(ax, 'g')
    
    subplot_allslices.canvas.restore_region(backgrounds_allslices)    
    line_max.set_ydata(maxes)
    line_mean.set_ydata(means)
    line_min.set_ydata(mins)

    subplot_allslices.canvas.draw()

    ax_allslices.set_ylim([0, np.max(maxes) + 10])
    ax_allslices.set_xlim([0, maxes.size - 1])
    ax_allslices.draw_artist(line_max)
    ax_allslices.draw_artist(line_mean)
    ax_allslices.draw_artist(line_min)
    subplot_allslices.canvas.blit(ax_allslices.bbox)
    
def get_all_slice_values():
  global Results_and_settings
  # Reads the CSO image, makes a velocity profile for all slices, and calculates the pulsatility value for all slices, given the use of mean/max velocity.
  if ctx.field("currenttab").value == 2:
    velocity_map = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.output0").image()
    velocity_map = velocity_map.getTile( (0, 0, 0, 0, 0), (velocity_map.UseImageExtent, velocity_map.UseImageExtent, velocity_map.UseImageExtent, velocity_map.UseImageExtent, velocity_map.UseImageExtent) )
    Results_and_settings['Area'] = np.array(np.count_nonzero(np.reshape(velocity_map, (velocity_map.shape[0], velocity_map.shape[2], -1)), axis = -1) * ctx.module("DTF_path_and_MPR").field("Info1.voxelSizeX").value * ctx.module("DTF_path_and_MPR").field("Info1.voxelSizeY").value)
    velocity_map = np.transpose(velocity_map, np.linspace(velocity_map.ndim - 1, 0, velocity_map.ndim).astype(int))
    
    current_velocity = np.zeros((velocity_map.shape[2], velocity_map.shape[4]))
    mean_flow = np.zeros((velocity_map.shape[2], velocity_map.shape[4]))
    cur_pos = np.zeros((velocity_map.shape[2], 6))
    cur_vec = np.zeros((velocity_map.shape[2], 3))
    for slice in range(0, velocity_map.shape[2]):
      cur_pos[slice, :] = np.array(ctx.module("DTF_path_and_MPR").field("XMarkerListContainer2.outXMarkerList").object().getMarkers()[slice].pos)
      cur_vec[slice, :] = np.array(ctx.module("DTF_path_and_MPR").field("XMarkerListContainer2.outXMarkerList").object().getMarkers()[slice].vec)
      velocity_ROI = velocity_map[:, :, slice, ...]
      flattened_v_ROI = velocity_ROI.reshape(-1, velocity_ROI.shape[-1])
      mean_flow[slice, :] = np.squeeze(np.sum(flattened_v_ROI, axis = 0))* ctx.module("DTF_path_and_MPR").field("Info1.voxelSizeX").value * ctx.module("DTF_path_and_MPR").field("Info1.voxelSizeY").value 
      if ctx.field("velocity_component").value == 2:
        current_velocity[slice, :] = np.squeeze(np.sum(flattened_v_ROI, axis = 0) / np.count_nonzero(flattened_v_ROI, axis = 0))
        Results_and_settings['Arithmetic'] = 'Mean'
      elif ctx.field("velocity_component").value == 1:
        Results_and_settings['Arithmetic'] = 'Max'
        current_velocity[slice, :] = np.squeeze(np.max(flattened_v_ROI, axis = 0))
      
    Results_and_settings['Velocity Profiles'] = current_velocity
    Results_and_settings['Centerline Positions'] = cur_pos
    Results_and_settings['Centerline Vectors'] = cur_vec
    Results_and_settings['Flow Profiles'] = mean_flow
    pulsatility = (np.max(current_velocity, axis = 1) - np.min(current_velocity, axis = 1)) / np.mean(current_velocity, axis = 1)
    Results_and_settings['Velocity pulsatility'] = pulsatility
    Results_and_settings['Path length'] = ctx.field("Arithmetic03.resultString").value
    for slice in range(0, velocity_map.shape[2]):
      ctx.module("Render_Pulsatility_spheres").field("Visualized_pulsatility_spheres.index").setValue(slice)
      if not np.isnan(pulsatility[slice]): 
        ctx.module("Render_Pulsatility_spheres").field("ComposeVector31.x").setValue(pulsatility[slice])
      else:
        ctx.module("Render_Pulsatility_spheres").field("ComposeVector31.x").setValue(0)
    ctx.module("Highlight_current_marker").field("XMarkerListContainer2.update").touch()
      
def setvelocityoverlay():
  # Toggles the velocity overlay in the magnitude MPR image
   oldslicevalue = ctx.module("Magnitude_velocity_viewers").field("MPR_velocity_ROI.startSlice").value
   ctx.module("Magnitude_velocity_viewers").module("Velocity_Overlay_Creation").field("SoView2DOverlay.drawingOn").setBoolValue(ctx.field("overlayswitch").value)
   ##ctx.module("Magnitude_velocity_viewers").field("FieldBypass.inputString0").setValue(oldslicevalue)
   
def csolistener(): 
  #Listens if a new CSO is drawn, or an existing one is edited. Will update the CSOtoimage module, so velocity values across the newly drawn/edited CSO contour are used.
  localcsolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
  slicelist = localcsolist.getGroupByLabel(ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").value)
  
  if slicelist.getNumCSOs() > (ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.maxTimePoint").value + 1):
    csolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
    temp, deleteindex = find_cso_at_time_index(csolist, ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").value, ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.timePoint").value)
    csolist.removeCSO(deleteindex)
    
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.apply").touch()
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage1.apply").touch()
  get_all_slice_values()
  plot_velocity_curve()
  plot_all_slice_velocities()
  plot_values_as_annotations()
    
def plot_values_as_annotations():
  #Sets the annotations in the magnitude MPR viewer to show the area index and pulsatility index of the current slice
  global Results_and_settings
  
  areas = Results_and_settings['Area']

  cur_area = areas[:,ctx.module("Magnitude_velocity_viewers").field('MPR_magnitude.startSlice').intValue()]
  area_index = (np.max(cur_area)-np.min(cur_area))/np.mean(cur_area)
  ctx.module("Magnitude_velocity_viewers").field("SoView2DAnnotation.numInput01").setValue(area_index)
  ctx.module("Magnitude_velocity_viewers").field("SoView2DAnnotation.numInput02").setValue(Results_and_settings['Velocity pulsatility'][ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").intValue()])

def set_anatomical_locations():
  #Stores selected Anatomical Location.
  global Results_and_settings
  Results_and_settings['leftright'] = ctx.field("leftright").value
  temp_locations = Results_and_settings['Anatomical Location']
  temp_locations[ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").intValue()]=ctx.field("Cname").intValue()
  Results_and_settings['Anatomical Location'] = temp_locations
  if np.any(temp_locations==0):
    ctx.field("Warning").setValue('At least one slice does not have a selected C.. location!')
  else:
    ctx.field("Warning").setValue("")
  
def get_anatomical_locations():
  #Sets displayed Anatomical Location
  global Results_and_settings
  temp_location = Results_and_settings['Anatomical Location']  
  ctx.field("Cname").setValue(temp_location[ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").intValue()])

def autofill():
#Autofills the anatomical location array. The first selected location is propagated until the slice for which the next location was selected, and so forth. 
  global Results_and_settings
  temp_location = Results_and_settings['Anatomical Location']  
  auto_filled = np.array([])
  for ind in temp_location:
    auto_filled = np.append(auto_filled, ind if ind else auto_filled[-1])
  Results_and_settings['Anatomical Location'] = auto_filled
  
def save_results():
#Wrapper to save all the results and settings into the user defined file.
  global Results_and_settings
  np.save(ctx.field("fileName").value + '.npy',Results_and_settings)
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("ImageSave.filename").setValue(ctx.field("fileName").value + '_cso_contour_velocity_image.dcm')
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("ImageSave.save").touch()
  
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOSave.fileName").setValue(ctx.field("fileName").value)
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOSave.startTaskSynchronous").touch()
  ctx.field("Warning").setValue("Saved")
  # Script to save the most important quantitative parameters to a csv file, usable in Excel.
  no_of_slices = Results_and_settings['Velocity Profiles'].shape[0]
  max_timepoint = Results_and_settings['Velocity Profiles'].shape[1]
  slice_list = ['Slice ' +str(slice) for slice in range(no_of_slices)]
  timepoint_list = np.linspace(0,max_timepoint-1,max_timepoint).astype(int).astype(str).tolist()
  with open(ctx.field("fileName").value+'.csv', 'w') as f:
    f.write('Original save location: '+ctx.field("fileName").value+'.csv')
    f.write('\n'*1)  
    if ctx.field("velocity_component").value==2:
      f.write('Arithmetic used: Mean')
    else:
      f.write('Arithmetic used: Max')
    f.write('\n'*1)  
    if Results_and_settings['leftright']=='Right':
      f.write('Measurements performed on right carotid')
    else:
      f.write('Measurements performed on left carotid')
    
    location_list = ['C' +str(location.astype(int)) for location in Results_and_settings['Anatomical Location']]
    csvappend = DataFrame(location_list,index=slice_list)
    f.write('\n'*1)  
    csvappend.to_csv(f, header=['Anatomical location per slice:'])    
    
    csvappend = DataFrame(np.transpose(Results_and_settings['Velocity pulsatility']),index=slice_list)
    f.write('\n'*1)
    csvappend.to_csv(f, header=['Velocity pulsatility'])
    
    csvappend = DataFrame(Results_and_settings['Velocity Profiles'],columns=timepoint_list,index=slice_list) 
    f.write('\n'*1)
    csvappend.to_csv(f,index_label=['Velocity (cm/s), Timepoints:'])
  
  
    csvappend = DataFrame(Results_and_settings['Flow Profiles'],columns=timepoint_list,index=slice_list)
    f.write('\n'*1)
    csvappend.to_csv(f, index_label=['Flow (mm3/s), Timepoints:'])
  
    csvappend = DataFrame(np.transpose(Results_and_settings['Area']),columns=timepoint_list,index=slice_list)
    f.write('\n'*1)
    csvappend.to_csv(f, index_label=['Area (mm2), Timepoints:'])
    
def find_cso_at_time_index(csolist, slice, timepoint):
  #Returns index in group, and Id of the cso.
  curslicelist = csolist.getGroupByLabel(slice)
  numcsos = curslicelist.getNumCSOs()
  for cso in range(0, numcsos):
    curslice = curslicelist.getCSOAt(cso)
    curtime = curslice.getTimePointIndex()
    if curtime == timepoint:
      return cso, curslice.getId()    
  return 0, curslicelist.getCSOAt(0).getId()
  
def propagate():
  # Redraws all the remaining timepoints of the current slice with a new starting location/isocontour isovalue. If a variable FWHM is chosen the isovalue is scaled as follows: isovalues_timepoints_new = current_timepoint_isoval/current_timepoint_old_isoval * isovalues_timepoints_old. 
  # If no variable FWHM was chosen the timepoints are redrawn with the isovalue selected for the current timepoint.
  global Results_and_settings, Iso_vals
  #Take maximum of magnitude image across all time points
  ctx.module("Magnitude_velocity_viewers").field("GetVoxelValue.update").touch()
  MPR_reconstruction_mag = ctx.field("DTF_path_and_MPR.output0").image()
  MPR_reconstruction_mag = MPR_reconstruction_mag.getTile( (0, 0, 0, 0, 0), (MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent, MPR_reconstruction_mag.UseImageExtent) )   
  slice = ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").intValue()
  cur_timepoint = ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.timePoint").intValue()
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").setValue(ctx.module("Magnitude_velocity_viewers").field("GetVoxelValue.outputValue").value)
  cur_magnitude = MPR_reconstruction_mag[0, ...]
  cso_areas = Results_and_settings['Area']
  cur_Iso_vals = Iso_vals[:,slice]*ctx.module("Magnitude_velocity_viewers").field("GetVoxelValue.outputValue").value/Iso_vals[cur_timepoint,slice]
  for timepoint in range(0, MPR_reconstruction_mag.shape[0]):
    ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").setValue(cur_Iso_vals[timepoint])
    #Generate contours
    ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.timePoint").setValue(timepoint)
    ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.timePoint").setValue(timepoint)
    ctx.module("Magnitude_velocity_viewers").field("ComposeVector3.z").setValue(slice)
    if not cur_magnitude[0, slice, np.round(MPR_reconstruction_mag.shape[4] / 2).astype(int), np.round(MPR_reconstruction_mag.shape[4] / 2).astype(int)] == 0:
      ctx.field("WorldVoxelConvert2.voxelPos").setStringValue(np.array2string(np.array([np.round(MPR_reconstruction_mag.shape[4] / 2), np.round(MPR_reconstruction_mag.shape[4] / 2), slice]))[1:-1])
      ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.startPosition").setStringValue(ctx.module("Magnitude_velocity_viewers").field("annoReadPix.worldPosition").stringValue())
      ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.addCSOToGroupWithLabel").setValue(slice)
      ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.apply").touch()
      
      csolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
      ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOInfo.csoShowByIndexOrId").setValue(csolist.getCSOIdList()[-1])
      cso_areas[timepoint, slice] = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOInfo.csoArea").value
      Iso_vals[timepoint, slice] = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOIsoGenerator2.isoValue").value
      
      localcsolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
      slicelist = localcsolist.getGroupByLabel(ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").value)
      if slicelist.getNumCSOs() > (ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.maxTimePoint").value + 1):
        csolist = ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOManager2.outCSOList").object()
        temp, deleteindex = find_cso_at_time_index(csolist, ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.startSlice").value, ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.timePoint").value)
        csolist.removeCSO(deleteindex)
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.apply").touch()
  ctx.module("Magnitude_velocity_viewers").field("MPR_magnitude.timePoint").setValue(0)
  Results_and_settings['Area'] = cso_areas 
  get_all_slice_values()
  plot_velocity_curve()
  plot_all_slice_velocities()
  plot_values_as_annotations()
  
  ctx.module("Magnitude_velocity_viewers").field("ImageClear.clear").touch()
  ctx.module("Magnitude_velocity_viewers").field("ImageClear.update").touch()
  ctx.module("Magnitude_velocity_viewers").field("annoReadPix.voxelPosition").setValue([0, 0, 0])

def apply_phase_inversion():
    try:
      phase_images = ctx.module("Load_split_dicoms").field("Phase_images.output0").image()
      phase_images = phase_images.getTile( (0, 0, 0, 0, 0, 0), (phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent) )   
      interface = ctx.module("Load_split_dicoms").module("PythonImage2").call("getInterface")
      interface.setImage(phase_images, minMaxValues = (np.min(phase_images), np.max(phase_images)), voxelToWorldMatrix = ctx.module("Load_split_dicoms").field("Info2.worldMatrix").value)
      # set view orientation/world matrix PCA
      ds = ctx.field("DirectDicomImport2.output0").getDicomTree()
      vieworientation = ds.getTag("(2001,105f)").getSequenceItem(0).getTag("(2005,1081)").value()
      Tmic = np.array([[ 0, -1, 0, 0],[ -1, 0, 0, 0],[ 0, 0 ,1, 0],[ 0, 0, 0, 1]])
      if vieworientation == 'AP':
        Tsom = np.array([[0, -1, 0, 0],[ 0, 0, 1, 0],[ 1, 0, 0, 0],[ 0, 0, 0, 1]])
      elif vieworientation == 'RL':
        Tsom = np.array([[0, 0, -1, 0],[ 0, -1, 0, 0],[ 1, 0, 0, 0],[ 0, 0, 0, 1]])
      elif vieworientation == 'FH':
        Tsom = np.array([[0, -1, 0, 0],[ -1, 0, 0, 0],[ 0, 0, 1, 0],[ 0, 0, 0, 1]])
      else:
        print('View angle not found in dicom header. Assuming transversal (FH) view angle')
        Tsom = np.array([[0 ,-1 ,0 ,0],[ -1, 0, 0, 0],[ 0, 0, 1, 0],[ 0 ,0, 0 ,1]])
      Tform = np.transpose( Tsom.dot(Tmic) )
      ctx.module("Phase_to_velocity").field("MatrixArithmetic2.matrixB").setValue(Tform)
      original_image = ctx.module("Phase_to_velocity").field("Bypass.output0").image()
      original_image = original_image.getTile((0, 0, 0, 0, 0,0), (original_image.UseImageExtent, original_image.UseImageExtent, original_image.UseImageExtent, original_image.UseImageExtent, original_image.UseImageExtent,original_image.UseImageExtent) )
      worldvoxelmat = np.array(ctx.module("Phase_to_velocity").field("MatrixArithmetic2.outputMatrixC").value)
      original_image_inworld = np.transpose(worldvoxelmat[0:3,0:3].dot(np.transpose(original_image,(5,4,3,2,0,1))),(0,5,4,3,2,1))
      interface = ctx.module("Phase_to_velocity").module("PythonImage").call("getInterface")
      interface.setImage(original_image_inworld, minMaxValues = (np.min(original_image_inworld), np.max(original_image_inworld))) 
    except:
      print('')
def apply_image_registration():
    # applies the image registration to all of the 4D PCA images.
    os.chdir(ctx.field("FileInformation3.dirname").value)
    ppca_filenames = glob.glob("*_PPCA*")
    mffe_filenames = glob.glob("*_MFFE*")
    mpca_filenames = glob.glob("*_MPCA*")
    all_files = ppca_filenames + mpca_filenames + mffe_filenames
    ctx.module("Load_split_dicoms").field("ImageLoad1.close").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad2.close").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad3.close").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad4.close").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad5.close").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad6.close").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad7.close").touch()
    
    for file in all_files:
      ctx.field("ImageLoad3.filename").setValue(ctx.field("FileInformation3.dirname").value + file )
      ctx.field("ImageSave3.save").touch()
      ctx.field("ImageLoad3.close").touch()
      shutil.move(ctx.field("FileInformation3.dirname").value + file + '.temp.dcm', ctx.field("FileInformation3.dirname").value + file)
      
    ctx.field("ImageLoad3.close").touch()
    ctx.field("ImageLoad3.filename").setValue("")
    ctx.module("Load_split_dicoms").field("ImageLoad1.load").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad2.load").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad3.load").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad4.load").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad5.load").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad6.load").touch()
    ctx.module("Load_split_dicoms").field("ImageLoad7.load").touch()
    phase_images = ctx.module("Load_split_dicoms").field("Phase_images.output0").image()
    phase_images = phase_images.getTile( (0, 0, 0, 0, 0, 0), (phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent, phase_images.UseImageExtent) )   
    interface = ctx.module("Load_split_dicoms").module("PythonImage2").call("getInterface")
    interface.setImage(phase_images, minMaxValues = (np.min(phase_images), np.max(phase_images)), voxelToWorldMatrix = ctx.module("Load_split_dicoms").field("Info2.worldMatrix").value)    
    
    magnitude = ctx.module("Load_split_dicoms").field("PythonArithmetic.output0").image()
    magnitude = magnitude.getTile( (0,0,0,0,0,0), (magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent, magnitude.UseImageExtent) )   
    
    interface = ctx.module("Load_split_dicoms").module("PythonImage").call("getInterface")
    interface.setImage(magnitude, minMaxValues = (np.min(magnitude), np.max(magnitude)), voxelToWorldMatrix = ctx.module("Load_split_dicoms").field("Info.worldMatrix").value)    

def cso_wem_computation():
# used for obtaining a visualization of the segmented vessel wall from the first timepoint isocontour per slice.
  ctx.module("SegmentationVisualization").field("CSOManager.copyInputCSOList").touch()
  csolist = ctx.module("SegmentationVisualization").field("CSOManager.outCSOList").object()
  numgroups = csolist.getNumGroups()
  for groupid in range(1,numgroups+1):
    group = csolist.getGroupById(groupid)
    cur_cso = group.getCSOAt(0)
    
    csolist.addCSOCopy(group.getCSOIdAt(0))
    cur_cso = csolist.getCSOAt(csolist.getNumCSOs()-1)
    cur_cso.setLabel(group.label)
    csolist.removeGroup(groupid)
  remainingIDs = csolist.getCSOIdList()
  
  for cur_id in remainingIDs:
    cur_cso = csolist.getCSOById(cur_id)
    ctx.module("SegmentationVisualization").field("MPRPath1.currentKeyFrame").setValue(int(cur_cso.getLabel()))
    ctx.module("SegmentationVisualization").field("CSOAffineTransformationModificator.csoIdList").setValue(str(cur_id))
    ctx.module("SegmentationVisualization").field("CSOAffineTransformationModificator.apply").touch()
    ctx.module("SegmentationVisualization").field("CSOAffineTransformationModificator1.apply").touch()
  ctx.module("SegmentationVisualization").field("CSOConvertTo3DMask.apply").touch()

def load_results():
  #Script linked to Load results button. Which will set the saved results into the GUI.
  global Results_and_settings
  ctx.field("Loaded_results_loc").value
  ctx.module("DTF_path_and_MPR").field("BaseBypass.bypass").setValue(False)
  Results_and_settings = np.load((ctx.field("Loaded_results_loc").value).replace("_cso_contour_velocity_image.dcm",".npy")).item()
  ctx.module("DTF_path_and_MPR").field("XMarkerListContainer1.deleteAll").touch()
  centers = Results_and_settings['Centerline Positions']
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer1.deleteAll").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer1.add").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer1.posXYZ").setValue(centers[0,:])
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer1.add").touch()
  ctx.module("Start_end_marker_selection").field("XMarkerListContainer1.posXYZ").setValue(centers[-1,:])
  for center in centers:
    ctx.module("DTF_path_and_MPR").field("XMarkerListContainer1.add").touch()
    ctx.module("DTF_path_and_MPR").field("XMarkerListContainer1.posXYZ").setValue(center)
  ctx.module("DTF_path_and_MPR").field("XMarkerShortestPath.startPoint").setValue(centers[0,:])
  ctx.module("DTF_path_and_MPR").field("XMarkerShortestPath.endPoint").setValue(centers[-1,:]) 
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOLoad.fileName").setValue((ctx.field("Loaded_results_loc").value).replace("_cso_contour_velocity_image.dcm",".cso"))
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOLoad.startTaskSynchronous").touch()
  
  #Set all radiobuttons
  
  if Results_and_settings['Arithmetic']=='Mean':
    ctx.field("velocity_component").value=2
  else:
    ctx.field("velocity_component").value=1 
  ctx.field("leftright").value = Results_and_settings['leftright']
  
  #ctx.module("DTF_path_and_MPR").field("PathToKeyFrame.numSmoothes").setValue(5)
  #set start-end markers  

  phase_to_velocity()
  ctx.module("Magnitude_velocity_viewers").module("CSO_creation").field("CSOConvertToImage.apply").touch()
  #Initialize plots
  plot_clear()
  plot_velocity_curve()
 
  
def running_median(matrix, width, axes):
  matrix[matrix==0]=np.nan
  vals=np.zeros((matrix.shape[axes],1))
  for dim in range(0,matrix.shape[axes]):
    orig_dim=int(dim)
    if dim==0:
      dim=(width-1)/2
    elif dim==matrix.shape[axes]-1:
      dim=matrix.shape[axes]-((width-1)/2)
    if axes==0:
      vals[orig_dim] = np.nanmedian(matrix[int(dim-(width-1)/2):int(dim+(width-1)/2),:])
    elif axes==1:
      vals[orig_dim] = np.nanmedian(matrix[:,int(dim-(width-1)/2):int(dim+(width-1)/2)])
  if axes:
    return np.tile(np.transpose(vals),[matrix.shape[0],1])
  else:
    return np.tile(vals,[1,matrix.shape[1]])
  
  
def shouldClose():  
  closeAllowed = True
  
  if MLAB.isStandaloneApplication():
    closeAllowed = MLAB.showQuestion("Are you sure you want to quit? Any unsaved data will be lost.", "Quit DampingUI", ['Yes', 'No'], 1) == 0
    
  ctx.window().setCloseAllowed(closeAllowed)
  
