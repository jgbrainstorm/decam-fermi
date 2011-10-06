#! /usr/bin/env python
import numpy as np
import pyfits as pf
import scipy.ndimage as nd
import pylab as pl
import sys

def getxy(data,title=None):
    #datas=nd.filters.gaussian_filter(data,5)
    datas=data
    ok=datas>=datas.mean()+1.5*datas.std()
    good=nd.binary_opening(ok,structure=np.ones((60,60)))    
    datagood=datas*good
    structuring_element=np.ones((3,3))
    segmentation,segments=nd.label(good,structure=structuring_element)    
    coords=nd.center_of_mass(datagood,segmentation,range(1,segments+1))
    xcoords=np.array([x[1] for x in coords])
    ycoords=np.array([x[0] for x in coords])
    pl.imshow(datagood,origin='lower')
    pl.plot(xcoords,ycoords,'g.')
    pl.xlim(0,data.shape[1])
    pl.ylim(0,data.shape[0])
    pl.title(title)
    return xcoords,ycoords

if len(sys.argv) == 1:
    print 'syntax: centroid_pos imageName1 imageName2 extension'
else:
    data1,hdr1=pf.getdata(sys.argv[1],int(sys.argv[3]),header=True)
    pl.figure(figsize=(14,9))
    pl.subplot(2,2,1)
    x1,y1=getxy(data1,sys.argv[1]+'['+sys.argv[3]+']')
    pl.subplot(2,2,2)
    data2=pf.getdata(sys.argv[2],int(sys.argv[3]))
    x2,y2=getxy(data2,sys.argv[2]+'['+sys.argv[3]+']')
    pl.subplot(2,2,3)
    pl.plot(x1,x2-x1,'bo')
    pl.xlabel('x1')
    pl.ylabel('x2-x1')
    pl.subplot(2,2,4)
    pl.plot(y1,y2-y1,'bo')
    pl.xlabel('y1')
    pl.ylabel('y2-y1')
    sep=np.sqrt((x1-x2)**2+(y1-y2)**2)
    print hdr1['detpos']+': '+str(np.round(sep,3))
    pl.show()
