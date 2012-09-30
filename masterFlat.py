#! /usr/bin/env python

"""
This code generates the master Flat image. each flat image is oversacn subtracted and then subtracted from the master bias. Then, the master flat is the median of these corrected flat images. 
J. Hao @ FNAL, 9/30/2012
"""
import numpy as np
import pyfits as pf
import sys
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
import glob as gl
from DECamCCD_def import *


if len(sys.argv) == 1:
    print 'syntax: '
    print 'For a given set of images: masterFlat FlatFileHead expids masterBiasName'
    print 'example: masterBias decam_09302012 (12345,12346,12347) masterBias.fits'
    print 'The resulting median image will be named as masterFlat.fits'
else:
    filehead = sys.argv[1]
    bias = sys.argv[2]
    expid = eval(sys.argv[2])
    nimg=len(expid)
    hdu = pf.open(filehead+str(expid[0])+'.fits.fz')
    for ext in range(1,63):
        for j in range(0,nimg):
            b=[]
            filename = filehead+'_'+str(expid[j])+'.fits.fz'
            imgext = pf.getdata(filename,ext)
            imgosub = oscanSub(imgext)
            imgosub = imgosub - pf.getdata(bias,ext)
            b.append(imgosub)
        hdu[ext].data=np.median(b,axis=0)
        hdu[ext].header.update('bzero',0)
    hdu.writeto('masterFlat.fits')
