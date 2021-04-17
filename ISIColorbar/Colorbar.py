# ----------------------------------------------------------------------------

# 
#  \file    Colorbar.py
#  \author  K.M. van Hespen
#  \date    2018-12-13
#
#  Color bar function

# ----------------------------------------------------------------------------

from mevis import *
import numpy as np
import os
import scipy.ndimage 
import sklearn.linear_model
import glob
import shutil
import matplotlib.pyplot as plt 
import matplotlib as mpl

def init():
  ctx.field("XMarkerListContainer6.deleteAll").touch()

def setlabels():
  ctx.field("XMarkerListContainer6.deleteAll").touch()
  ticklabels = np.linspace(ctx.field("LUTInfo.maxIndex").value,ctx.field("LUTInfo.minIndex").value,ctx.field("tickcount").value)
  ticklabels = np.around(ticklabels, decimals=2)
  ticklocations = np.linspace(20,ctx.field("LUTToMLImage.rescaleWidth").value-20,ctx.field("tickcount").value)
  
  for tick in range(0,ctx.field("tickcount").value):
    currentloc = ticklocations[tick]
    if tick==0:
      currentloc = ticklocations[tick]
    if tick==np.max(range(0,ctx.field("tickcount").value)):
      currentloc = ticklocations[tick]
    ctx.field("XMarkerListContainer6.add").touch()
    ctx.field("WorldVoxelConvert4.voxelY").setValue(currentloc)
    ctx.field("XMarkerListContainer6.name").setStringValue(ticklabels[tick])
  ctx.field("XMarkerListContainer6.add").touch()
  