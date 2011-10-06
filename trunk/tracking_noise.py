#! /usr/bin/env python
from fermiMCCD import *

dir=
NameFits=gl.glob(dir+'*.fits')
NameFits.sort()

Nfile=len(NameFits)
noise=np.zeros(Nfile)
Nchannel=len(NameFits[0])

for j in range(1,Nchannel):
    for i in range(0:Nfile):
        noise[i]=noise_channel(NameFits[i],Channel[j])
    pl.plot(noise,'bo')
    pl.xlabel('Sequence Number')
    pl.ylabel('noise (ADU)')
    pl.title('extension: '+str(Channel[j]))
    pl.savefig(dir+'fig/noise_channel_'+str(Channel[j])+'.png')
