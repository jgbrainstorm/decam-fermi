#! /usr/bin/env python

"""
This code performe the bias and flat field correction. All the frames are overscan subtracted. Before you perform this, you need to have masterBias and masterFlat readay
J. Hao @ FNAL, 9/30/2012
"""
import numpy as np
import pyfits as pf
import sys,time
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
import glob as gl
from DECamCCD_def import *
from DECamCCD import *


if len(sys.argv) == 1:
    print 'syntax: '
    print '   desImgReduction masterBias masterFlat ImgFileHead epxid  '
    print 'example:' 
    print '   desImgReduction masterBias.fits masterFlat.fits DECam 001234 0023456'

else:
    startTime=time.time()
    bias = sys.argv[1]
    flat = sys.argv[2]
    filehead = sys.argv[3]
    nimg=len(sys.argv) - 4
    
    for i in range(nimg):
        hdu = pf.open(filehead+'_'+sys.argv[4+i]+'.fits.fz') # need to funpack first
        for ext in range(1,63):
            print ext
            hdu[ext].data = (oscanSub(hdu[ext].data) - pf.getdata(bias,ext))/pf.getdata(flat,ext)
        hdu.write(filehead+'_'+sys.argv[4+i]+'_corrected.fits.fz')
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)

    
