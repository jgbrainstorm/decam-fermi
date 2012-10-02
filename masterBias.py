#! /usr/bin/env python

"""
This code generate the master bias image after subtracting the overscan of each bias frame only for the imaging ccd area, not the focus ccd.
J. Hao @ FNAL, 9/30/2012
"""
import numpy as np
import pyfits as pf
import sys
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
import glob as gl
from DECamCCD import *
from DECamCCD_def import *
import time

if len(sys.argv) == 1:
    print 'syntax: '
    print '    masterBias.py biasFileHead expids'
    print 'example: '
    print '    masterBias.py DECam 12345 12346 12347'
    print 'The resulting median image will be named as masterBias.fits'
else:
    startTime=time.time()
    if sys.argv[2] == 'all':
        filenamelist = gl.glob('*.fits')
        nimg = len(filenamelist)
        hdu = pf.open(filenamelist[0],mode='update')
    else:
        filehead = sys.argv[1]
        nimg=len(sys.argv) - 2
        hdu = pf.open(filehead+'_'+sys.argv[2]+'.fits',mode='update') # need to funpack first
    hdu[0].header.update('PROCTYPE','Master Bias')
    for ext in range(1,63):
        print ext
        for j in range(0,nimg):
            b=[]
            if sys.argv[2] == 'all':
                filename = filenamelist[j]
            else:
                filename = filehead+'_'+sys.argv[2+j]+'.fits'
            imgext = pf.getdata(filename,ext)
            imgosub = oscanSub(imgext)
            b.append(imgosub)
        hdu[ext].data=np.median(b,axis=0)
        hdu[ext].header.update('bzero',0)
    endTime=time.time()
    elapseTime=endTime-startTime
    hdu.writeto('masterBias.fits')
    print '---elapsed time: ' + str(elapseTime)

