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
from DECamCCD import *

if len(sys.argv) == 1:
    print 'syntax: '
    print '   masterFlat.py masterBiasName FlatFileHead expids '
    print 'or:'
    print '   masterFlat.py masterBiasName all'
    print 'example: '
    print 'masterFlat.py masterBias.fits decam 12345 12346 12347'
    print 'masterFlat.py masterBias.fits all'
    print '   The resulting median image will be named as masterFlat.fits'
else:
    startTime=time.time()
    bias = sys.argv[1]
    sumpix = 0.
    if sys.argv[2] == 'all':
        filenamelist = gl.glob('*.fits')
        nimg = len(filenamelist)
        hdu = pf.open(filenamelist[0],mode='update')
    else:   
        filehead = sys.argv[2]
        nimg=len(sys.argv) - 3
        hdu = pf.open(filehead+'_'+sys.argv[3]+'.fits',mode='update') # need to funpack first
    hdu[0].header.update('PROCTYPE','Master Flat')
    for ext in range(1,71):
        print ext
        b=[]
        for j in range(0,nimg):
            if sys.argv[2] == 'all':
                filename = filenamelist[j]
            else:
                filename=filehead+'_'+sys.argv[3+j]+'.fits'
            imgext = pf.getdata(filename,ext)
            imgosub = oscanSub(imgext)
            imgosub = imgosub - pf.getdata(bias,ext)    
            b.append(imgosub)
        col0=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[0].split(':')[0])-1
        col1=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[0].split(':')[1]) 
        row0=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[1].split(':')[0])-1
        row1=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[1].split(':')[1]) 
        hdu[ext].data=np.median(b,axis=0)
        hdu[ext].data[hdu[ext].data == 0.] = 0.0000001 #avoid the blow up
        sumpix = sumpix+np.sum(hdu[ext].data[row0:row1,col0:col1])
        hdu[ext].header.update('bzero',0)
    globalMean = sumpix / (4096.*2048*62+4096.*1024*8)
    for ext in range(1,71):
        col0=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[0].split(':')[0])-1
        col1=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[0].split(':')[1]) 
        row0=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[1].split(':')[0])-1
        row1=int(hdu[ext].header['datasec'].split('[')[1].split(']')[0].split(',')[1].split(':')[1]) 
        hdu[ext].data[row0:row1,col0:col1] = hdu[ext].data[row0:row1,col0:col1]/globalMean
    hdu.writeto('masterFlat.fits')
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)

