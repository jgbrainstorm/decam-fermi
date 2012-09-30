#! /usr/bin/env python

"""
This code generates the master Flat image. each flat image is oversacn subtracted and then subtracted from the master bias. Then, the master flat is the median of these corrected flat images. 
J. Hao @ FNAL, 9/30/2012
"""
import numpy as np
import pyfits as pf
import sys,time
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
import glob as gl
from DECamCCD_def import *



if len(sys.argv) == 1:
    print 'syntax: '
    print 'For a given set of images: masterFlat masterBiasName FlatFileHead expids '
    print 'example: masterFlat masterBias.fits decam 12345 12346 12347'
    print 'The resulting median image will be named as masterFlat.fits'
else:
    startTime=time.time()
    bias = sys.argv[1]
    filehead = sys.argv[2]
    nimg=len(sys.argv) - 3
    hdu = pf.open(filehead+'_'+sys.argv[3]+'.fits') # need to funpack first
    for ext in range(1,63):
        for j in range(0,nimg):
            b=[]
            filename = filehead+'_'+sys.argv[2+j]+'.fits'
            imgext = pf.getdata(filename,ext)
            imgosub = oscanSub(imgext)
            imgosub = imgosub - pf.getdata(bias,ext)    
            b.append(imgosub)
        hdu[ext].data=np.median(b,axis=0)
        hdu[ext].data = hdu[ext].data / robust_mean(hdu[ext].data)
        hdu[ext].header.update('bzero',0)
    hdu.writeto('masterFlat.fits')
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)
    return '----done ---'
