#! /usr/bin/env python
# this one is for the scrumbled exposure sequence


import sys,os
sys.path.append('/home/jghao/fmccd')
from DECamCCD import *
import cPickle as p

ini_expid = 199076
end_expid = 199237
nfile = end_expid - ini_expid +1

expid = np.arange(ini_expid,end_expid+1)

dir = '/blue-orch/DTS/src/20130419/src/'
exptime =np.zeros(nfile) 
NameFits = []
for i in range(nfile):
    NameFits.append(dir+'DECam_00'+str(expid[i])+'.fits.fz')
    exptime[i] = pf.getheader(dir+'DECam_00'+str(expid[i])+'.fits.fz',0)['exptime']
            
NameFits = np.array(NameFits)

if not os.path.isfile('median.fits'):
    os.system('medianImg '+NameFits[0]+' '+NameFits[1]+' '+NameFits[2]+' '+NameFits[3]+' '+NameFits[4])

NameFits = NameFits[5:]
exptime = exptime[5:]

idxSeq = (exptime != 0)*(exptime != 2)*(exptime != 0.05)
exptime = exptime[idxSeq]
NameFits = NameFits[idxSeq]

idxsort = np.argsort(exptime)
NameFits = NameFits[idxsort]

NameBias = 'median.fits'
hdu=pf.open(NameBias)
Channel=range(1,len(hdu))

ccd=[]
fullwell=[]
gain=[]

fullwellData = []
for i in Channel:
    print '----- channel: ',i,'-----'
    data,hdr=pf.getdata(NameBias,i,header=True)
    oscanL=np.mean(data[500:1500,10:50])
    #detname=hdr['detser']
    detname=hdr['detpos']
    if oscanL != 0:     
        resL=linearity(NameFits,NameBias,i,shift=0,left=1)
        pl.savefig('fig/ptc_'+detname+'_channel_'+str(i)+'_left.png')
        ccd.append(detname+'L')
        fullwell.append(resL[1])
        gain.append(resL[0])
        fullwellData.append(resL)
    oscanR=np.mean(data[500:1500,2110:2150])
    if oscanR != 0:
        resR=linearity(NameFits,NameBias,i,shift=0,left=0)
        pl.savefig('fig/ptc_'+detname+'_channel_'+str(i)+'_right.png')
        ccd.append(detname+'R')
        fullwell.append(resR[1])
        gain.append(resR[0])
        fullwellData.append(resR)

fullwell=np.array(fullwell)
gain=np.array(gain)
f = open('fullwell.txt', 'w')
for j in range(len(ccd)):
    f.write(ccd[j]+'   '+str(round(fullwell[j]))+'   '+str(round(gain[j],4))+'\n')
f.close()
    
p.dump(fullwellData,open('fullwellData.p','w'))
    

print '---- linearity analysis complete ----'
