#! /usr/bin/env python

"""
This code perform only the overscan subtraction for each image extension, no masterbias, master flat are used. 
J. Hao @ FNAL, 11/28/2012
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
    print '   desImgQuickReduction.py ImgFileHead epxid  '
    print 'or:'
    print '   desImgQuickReduction.py all'
    print '   which process all the files in the current directory'
    print 'example:' 
    print '   desImgQuickReduction.py DECam 001234 0023456'
    print '   desImgQuickReduction.py all'

else:
    startTime=time.time()
    if sys.argv[1] == 'all':
        filename = gl.glob('*.fits')
        nimg = len(filename)
    else:
        filehead = sys.argv[1]
        nimg=len(sys.argv) - 2
    for i in range(nimg):
        if sys.argv[1] == 'all':
            hdu = pf.open(filename[i],mode='update')
            correctedFilename = filename[i][0:-5]+'_reduced.fits'
        else:
            hdu = pf.open(filehead+'_'+sys.argv[2+i]+'.fits',mode='update')
            correctedFilename = filehead+'_'+sys.argv[2+i]+'_reduced.fits'
        hdu[0].header.update('PROCTYPE','Reduced')
        for ext in range(1,71):
            print ext
            hdu[ext].data = oscanSub(hdu[ext].data
        hdu.writeto(correctedFilename)
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)

    
