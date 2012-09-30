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
from DECamCCD_def import *

if len(sys.argv) == 1:
    print 'syntax: '
    print 'masterBias biasFileHead, expids'
    print 'example: masterBias decam (12345,12346,12347)'
    print 'The resulting median image will be named as masterBias.fits'
else:
    filehead = sys.argv[1]
    expid = eval(sys.argv[2])
    nimg=len(expid)
    hdu = pf.open(filehead+'_'+str(expid[0])+'.fits.fz') # need to funpack first
    for ext in range(1,63):
        for j in range(0,nimg):
            b=[]
            filename = filehead+str(expid[j])+'.fits'
            imgext = pf.getdata(filename,ext)
            imgosub = oscanSub(imgext)
            b.append(imgosub)
        hdu[ext].data=np.median(b,axis=0)
        hdu[ext].header.update('bzero',0)
    hdu.writeto('masterBias.fits')
