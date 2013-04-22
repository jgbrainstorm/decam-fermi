#! /usr/bin/env python
import sys
sys.path.append('/home/jghao/fmccd')
from DECamCCD import *
import cPickle as p

dir=os.getcwd()+'/'

NameFits=gl.glob(dir+'/*.fits*')
NameBias=dir+'/bias/median.fits'
NameFits.sort()
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
        pl.savefig(dir+'fig/ptc_'+detname+'_channel_'+str(i)+'_left.png')
        ccd.append(detname+'L')
        fullwell.append(resL[1])
        gain.append(resL[0])
        fullwellData.append(resL)
    oscanR=np.mean(data[500:1500,2110:2150])
    if oscanR != 0:
        resR=linearity(NameFits,NameBias,i,shift=0,left=0)
        pl.savefig(dir+'fig/ptc_'+detname+'_channel_'+str(i)+'_right.png')
        ccd.append(detname+'R')
        fullwell.append(resR[1])
        gain.append(resR[0])
        fullwellData.append(resR)
fullwell=np.array(fullwell)
gain=np.array(gain)
f = open(dir+'fullwell.txt', 'w')
for j in range(len(ccd)):
    f.write(ccd[j]+'   '+str(round(fullwell[j]))+'   '+str(round(gain[j],4))+'\n')
f.close()
    
p.dump(fullwellData,open('fullwellData.p','w'))
    

print '---- linearity analysis complete ----'
