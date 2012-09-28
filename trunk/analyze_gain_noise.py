import pylab as pl
import numpy as np

#----noise ----

noisedata = np.genfromtxt('noise.txt',dtype=None,names=['ccdname','noise'])
ccdname = noisedata['ccdname']
noise = noisedata['noise']

n = len(noise)
pl.figure(figsize=(16,10))
pl.subplot(2,1,1)
pl.bar(np.arange(n/2),noise[0:n/2],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccdname[0:n/2],rotation=90)
pl.ylim(0,3)
pl.ylabel('Readout Noise (ADU)')

pl.subplot(2,1,2)
pl.bar(np.arange(n/2),noise[n/2:],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccdname[n/2:],rotation=90)
pl.ylim(0,3)
pl.ylabel('Readout Noise (ADU)')
pl.savefig('readout_noise.png')
pl.close()

#---gain and fullwell

data = np.genfromtxt('fullwell.txt',dtype=None,names=['ccd','fw','gain'])
fw = data['fw']
ccd = data['ccd']
gain = data['gain']

n = len(gain)

pl.figure(figsize=(16,10))
pl.subplot(2,1,1)
pl.bar(np.arange(n/2),fw[0:n/2],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccd[0:n/2],rotation=90)
pl.ylim(0,60000)
pl.ylabel('Full Well in ADU')

pl.subplot(2,1,2)
pl.bar(np.arange(n/2),fw[n/2:],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccd[n/2:],rotation=90)
pl.ylim(0,60000)
pl.ylabel('Full Well in ADU')
pl.savefig('fullwell_adu.png')
pl.close()

pl.figure(figsize=(16,10))
pl.subplot(2,1,1)
pl.bar(np.arange(n/2),fw[0:n/2]/gain[0:n/2],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccd[0:n/2],rotation=90)
pl.ylim(0,250000)
pl.hlines(130000,0,n/2,color='blue',label='spec')
pl.ylabel('Full Well in electron')
pl.legend(loc='best')

pl.subplot(2,1,2)
pl.bar(np.arange(n/2),fw[n/2:]/gain[0:n/2],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccd[n/2:],rotation=90)
pl.ylim(0,250000)
pl.ylabel('Full Well in electron')
pl.hlines(130000,0,n/2,color='blue',label='spec')
pl.legend(loc='best')
pl.savefig('fullwell_electron.png')
pl.close()


pl.figure(figsize=(16,10))
pl.subplot(2,1,1)
pl.bar(np.arange(n/2),1./gain[0:n/2],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccd[0:n/2],rotation=90)
pl.ylim(0,5)
pl.ylabel('Gain (e-/ADU)')

pl.subplot(2,1,2)
pl.bar(np.arange(n/2),1./gain[n/2:],color=['pink','green'])
pl.xticks(np.arange(n/2)+0.5,ccd[n/2:],rotation=90)
pl.ylim(0,5)
pl.ylabel('Gain (e-/ADU)')
pl.savefig('gain.png')
pl.close()
