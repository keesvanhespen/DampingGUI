#!/usr/bin/env python

"""
Get subimages where the border is mirrored.

Written by:
Kees van Hespen, UMC Utrecht, The Netherlands
"""
from mevis import *
import numpy as np


def autorun():
  if ctx.field("FieldListener1.sourceFieldValue").value=='Valid':
    if ctx.field("Arithmetic08.resultX").value==0:
      ctx.field("PythonArithmetic.update").touch()
    else:
      ctx.field("PythonArithmetic1.update").touch()