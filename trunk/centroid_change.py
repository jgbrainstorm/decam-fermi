#! /usr/bin/env python
import numpy as np
import pyfits as pf
import scipy.ndimage as nd
import pylab as pl
import sys


def getxy(datas):
    ok=datas>=datas.mean()+1.5*datas.std()
    good=nd.binary_opening(ok,structure=np.ones((60,60)))    
    datagood=datas*good
    structuring_element=np.ones((3,3))
    segmentation,segments=nd.label(good,structure=structuring_element)    
    coords=nd.center_of_mass(datagood,segmentation,range(1,segments+1))
    xcoords=np.array([x[1] for x in coords])
    ycoords=np.array([x[0] for x in coords])
    return xcoords,ycoords

if len(sys.argv) == 1:
    print 'syntax: centroid_pos imageName1 imageName2'
else:
    ext=[2,3,4,5,19,20,21,37,38,39,41,42,43,47,48,49,52,45]
    #ext=[2,37,41]
    print '------'+sys.argv[1]+' vs. '+sys.argv[2]+'----'
    print 'CCD Position------dot change in pixels'
    for i in ext:
        data1,hdr1=pf.getdata(sys.argv[1],i,header=True)
        x1,y1=getxy(data1)
        data2,hdr2=pf.getdata(sys.argv[2],i,header=True)
        x2,y2=getxy(data2)
        sep=np.sqrt((x1-x2)**2+(y1-y2)**2)
        print hdr1['detpos']+' x2-x1 :'+str(np.round(x2-x1,3))
        print hdr1['detpos']+' y2-y1 :'+str(np.round(y2-y1,3))
        print hdr1['detpos']+' separation: '+str(np.round(sep,3))
