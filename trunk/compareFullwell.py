import pylab as pl
import numpy as np

data1 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_12_7_2011.txt',dtype=None,names=['ccd','fw','gain'])
data2 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_7_26_2012.txt',dtype=None,names=['ccd','fw','gain'])
data3 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_9_27_2012.txt',dtype=None,names=['ccd','fw','gain'])


fw1 = data1['fw']
ccd1 = data1['ccd']
gain1 = data1['gain']

idx1 = np.argsort(ccd1)
fw1 = fw1[idx1]
ccd1 = ccd1[idx1]
gain1 = gain1[idx1]
fwe1 = fw1/gain1


fw2 = data2['fw']
ccd2 = data2['ccd']
gain2 = data2['gain']

idx2 = np.argsort(ccd2)
fw2 = fw2[idx2]
ccd2 = ccd2[idx2]
gain2 = gain2[idx2]
fwe2 = fw2/gain2


fw3 = data3['fw']
ccd3 = data3['ccd']
gain3 = data3['gain']

idx3 = np.argsort(ccd3)
fw3 = fw3[idx3]
ccd3 = ccd3[idx3]
gain3 = gain3[idx3]
fwe3 = fw3/gain3

n = len(gain1)

pl.figure(figsize=(20,10))
pl.subplot(2,1,1)
pl.bar(np.arange(n/2),fwe1[0:n/2],width=0.2,color='pink',label='12/7/2012')
pl.bar(np.arange(n/2)+0.2,fwe2[0:n/2],width=0.2,color='green',label='7/26/2012')
pl.bar(np.arange(n/2)+0.4,fwe3[0:n/2],width=0.2,color='blue',label='9/27/2012')
pl.xticks(np.arange(n/2)+0.5,ccd1[0:n/2],rotation=90)
pl.hlines(130000,0,n/2,color='red',label='spec')
pl.legend(loc = 'best')
pl.ylim(-3000,250000)
pl.ylabel('Fullwell (e-)')

pl.subplot(2,1,2)
pl.bar(np.arange(n/2),fwe1[0:n/2],width=0.2,color='pink',label='12/7/2012')
pl.bar(np.arange(n/2)+0.2,fwe2[0:n/2],width=0.2,color='green',label='7/26/2012')
pl.bar(np.arange(n/2)+0.4,fwe3[0:n/2],width=0.2,color='blue',label='9/27/2012')
pl.xticks(np.arange(n/2)+0.5,ccd1[0:n/2],rotation=90)
pl.hlines(130000,0,n/2,color='red',label='spec')
pl.legend(loc = 'best')
pl.ylim(-3000,250000)
pl.ylabel('Fullwell (e-)')
pl.savefig('fullwell_comp.png')

#---in adu---

pl.figure(figsize=(20,10))
pl.subplot(2,1,1)
pl.bar(np.arange(n/2),fw1[0:n/2],width=0.2,color='pink',label='12/7/2011')
pl.bar(np.arange(n/2)+0.2,fw2[0:n/2],width=0.2,color='green',label='7/26/2012')
pl.bar(np.arange(n/2)+0.4,fw3[0:n/2],width=0.2,color='blue',label='9/27/2012')
pl.xticks(np.arange(n/2)+0.5,ccd1[0:n/2],rotation=90)
pl.legend(loc = 'best')
pl.ylim(0,50000)
pl.ylabel('Fullwell (ADU)')

pl.subplot(2,1,2)
pl.bar(np.arange(n/2),fw1[0:n/2],width=0.2,color='pink',label='12/7/2011')
pl.bar(np.arange(n/2)+0.2,fw2[0:n/2],width=0.2,color='green',label='7/26/2012')
pl.bar(np.arange(n/2)+0.4,fw3[0:n/2],width=0.2,color='blue',label='9/27/2012')
pl.xticks(np.arange(n/2)+0.5,ccd1[0:n/2],rotation=90)
pl.legend(loc = 'best')
pl.ylim(0,50000)
pl.ylabel('Fullwell (ADU)')
pl.savefig('fullwell_comp_adu.png')
