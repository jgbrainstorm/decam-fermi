#! /usr/bin/env python

import sys
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
from DECamCCD import *

file = gl.glob('*.fits')
file.sort()
bias = '/home3/data_local/images/fits/masterBias/masterBias.fits'

def getSec(hdr=None,ext=None,datasec=None):
    col0=int(hdr[datasec].split('[')[1].split(']')[0].split(',')[0].split(':')[0])-1
    col1=int(hdr[datasec].split('[')[1].split(']')[0].split(',')[0].split(':')[1]) 
    row0=int(hdr[datasec].split('[')[1].split(']')[0].split(',')[1].split(':')[0])-1
    row1=int(hdr[datasec].split('[')[1].split(']')[0].split(',')[1].split(':')[1]) 
    return row0,row1,col0,col1


nfile =len(file)
for ext in range(1,63):
    exptime = np.zeros(nfile)
    abratio = np.zeros(nfile)
    medcount = np.zeros(nfile)
    print ext
    for i in range(0,nfile):
        imgext = pf.getdata(file[i],ext)
        hdr = pf.getheader(file[i],ext)
        imgosub = oscanSub(imgext)
        imgosub = imgosub - pf.getdata(bias,ext)
        row0,row1,col0,col1=getSec(hdr,ext,'dataseca')
        dataA = imgosub[row0:row1,col0:col1]
        row0,row1,col0,col1=getSec(hdr,ext,'datasecb')
        dataB = imgosub[row0:row1,col0:col1]
        abratio[i] =np.median(dataA)/np.median(dataB) 
        medcount[i] = np.median(dataB)
        exptime[i] = pf.getheader(file[i],0)['exptime']
    idx1 = np.arange(0,nfile,2)
    idx2 = np.arange(1,nfile+1,2)
    pl.subplot(2,1,1)
    pl.plot(exptime[idx1],abratio[idx1],'bo')
    pl.plot(exptime[idx2],abratio[idx2],'r*')
    pl.xlabel('exptime')
    pl.ylabel('A/B')
    pl.title(hdr['detpos'])
    pl.ylim(0.9,1.1)
    pl.subplot(2,1,2)
    pl.plot(medcount[idx1],abratio[idx1],'bo')
    pl.plot(medcount[idx2],abratio[idx2],'r*')
    pl.xlabel('Median count of amp A')
    pl.ylabel('A/B')
    pl.ylim(0.9,1.1)
    pl.savefig(hdr['detpos']+'.png')
    pl.close()
    
