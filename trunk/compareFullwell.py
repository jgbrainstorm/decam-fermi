#! /usr/bin/env python
import pylab as pl
import numpy as np
import cPickle as p

data1 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_12_7_2011.txt',dtype=None,names=['ccd','fw','gain'])
data2 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_7_26_2012.txt',dtype=None,names=['ccd','fw','gain'])
data3 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_9_27_2012.txt',dtype=None,names=['ccd','fw','gain'])
data4 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_10_02_2012.txt',dtype=None,names=['ccd','fw','gain'])
data5 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_10_05_2012_old_seqr.txt',dtype=None,names=['ccd','fw','gain'])
data6 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_10_6_2012_veryshort.txt',dtype=None,names=['ccd','fw','gain'])
data7 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_11_09_2012.txt',dtype=None,names=['ccd','fw','gain'])
data8 = np.genfromtxt('/home/jghao/research/ccd/CTI_new/new_in_ctio/fullwell_11_19_2012.txt',dtype=None,names=['ccd','fw','gain'])


fw1 = data1['fw']
ccd1 = data1['ccd']
gain1 = data1['gain']
fwe1 = fw1/gain1


fw2 = data2['fw']
ccd2 = data2['ccd']
gain2 = data2['gain']
fwe2 = fw2/gain2


fw3 = data3['fw']
ccd3 = data3['ccd']
gain3 = data3['gain']
fwe3 = fw3/gain3

fw4 = data4['fw']
ccd4 = data4['ccd']
gain4 = data4['gain']
fwe4 = fw4/gain4

fw5 = data5['fw']
ccd5 = data5['ccd']
gain5 = data5['gain']
fwe5 = fw5/gain5


fw6 = data6['fw']
ccd6 = data6['ccd']
gain6 = data6['gain']
fwe6 = fw6/gain6

fw7 = data7['fw']
ccd7 = data7['ccd']
gain7 = data7['gain']
fwe7 = fw7/gain7

fw8 = data8['fw']
ccd8 = data8['ccd']
gain8 = data8['gain']
fwe8 = fw8/gain8


n = len(gain1)

pl.figure(figsize=(20,20))
pl.subplot(4,1,1)
pl.bar(np.arange(n/4),fwe1[0:n/4],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fwe2[0:n/4],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fwe3[0:n/4],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fwe4[0:n/4],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fwe5[0:n/4],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fwe7[0:n/4],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fwe8[0:n/4],width=0.1,color='c',label='11/19/2012')

pl.xticks(np.arange(n/4)+0.5,ccd1[0:n/4],rotation=90)
pl.hlines(130000,0,n/4,color='red',label='spec')
pl.legend(loc = 'upper right',bbox_to_anchor=(1.,1.5))
pl.ylim(-3000,300000)
pl.grid()
pl.ylabel('Fullwell (e-)')

pl.subplot(4,1,2)
pl.bar(np.arange(n/4),fwe1[n/4:n/2],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fwe2[n/4:n/2],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fwe3[n/4:n/2],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fwe4[n/4:n/2],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fwe5[n/4:n/2],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fwe7[n/4:n/2],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fwe8[n/4:n/2],width=0.1,color='c',label='11/19/2012')

pl.xticks(np.arange(n/4)+0.5,ccd1[n/4:n/2],rotation=90)
pl.hlines(130000,0,n/4,color='red',label='spec')
#pl.legend(loc = 'best')
pl.ylim(-3000,300000)
pl.ylabel('Fullwell (e-)')
pl.grid()

pl.subplot(4,1,3)
pl.bar(np.arange(n/4),fwe1[n/4+n/4:n/2+n/4],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fwe2[n/4+n/4:n/2+n/4],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fwe3[n/4+n/4:n/2+n/4],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fwe4[n/4+n/4:n/2+n/4],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fwe5[n/4+n/4:n/2+n/4],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fwe7[n/4+n/4:n/2+n/4],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fwe8[n/4+n/4:n/2+n/4],width=0.1,color='c',label='11/19/2012')
pl.xticks(np.arange(n/4)+0.5,ccd1[n/4+n/4:n/2+n/4],rotation=90)
pl.hlines(130000,0,n/4,color='red',label='spec')
#pl.legend(loc = 'best')
pl.ylim(-3000,300000)
pl.ylabel('Fullwell (e-)')
pl.grid()

pl.subplot(4,1,4)
pl.bar(np.arange(n/4),fwe1[n/2+n/4:n/2+n/2],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fwe2[n/2+n/4:n/2+n/2],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fwe3[n/2+n/4:n/2+n/2],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fwe4[n/2+n/4:n/2+n/2],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fwe5[n/2+n/4:n/2+n/2],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fwe7[n/2+n/4:n/2+n/2],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fwe8[n/2+n/4:n/2+n/2],width=0.1,color='c',label='11/19/2012')
pl.xticks(np.arange(n/4)+0.5,ccd1[n/2+n/4:n/2+n/2],rotation=90)
pl.hlines(130000,0,n/4,color='red',label='spec')
#pl.legend(loc = 'best')
pl.ylim(-3000,300000)
pl.ylabel('Fullwell (e-)')
pl.grid()

pl.savefig('fullwell_comp.png')
pl.close()






#---in adu---

pl.figure(figsize=(20,20))

pl.subplot(4,1,1)
pl.bar(np.arange(n/4),fw1[0:n/4],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fw2[0:n/4],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fw3[0:n/4],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fw4[0:n/4],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fw5[0:n/4],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fw7[0:n/4],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fw8[0:n/4],width=0.1,color='c',label='11/19/2012')

pl.xticks(np.arange(n/4)+0.5,ccd1[0:n/4],rotation=90)
pl.legend(loc = 'upper right',bbox_to_anchor=(1.,1.5))
pl.ylim(0,60000)
pl.grid()
pl.ylabel('Fullwell (ADU)')

pl.subplot(4,1,2)
pl.bar(np.arange(n/4),fw1[n/4:n/2],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fw2[n/4:n/2],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fw3[n/4:n/2],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fw4[n/4:n/2],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fw5[n/4:n/2],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fw7[n/4:n/2],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fw8[n/4:n/2],width=0.1,color='c',label='11/19/2012')

pl.xticks(np.arange(n/4)+0.5,ccd1[n/4:n/2],rotation=90)
#pl.legend(loc = 'best')
pl.ylim(0,60000)
pl.ylabel('Fullwell (ADU)')
pl.grid()

pl.subplot(4,1,3)
pl.bar(np.arange(n/4),fw1[n/4+n/4:n/2+n/4],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fw2[n/4+n/4:n/2+n/4],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fw3[n/4+n/4:n/2+n/4],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fw4[n/4+n/4:n/2+n/4],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fw5[n/4+n/4:n/2+n/4],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fw7[n/4+n/4:n/2+n/4],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fw8[n/4+n/4:n/2+n/4],width=0.1,color='c',label='11/19/2012')
pl.xticks(np.arange(n/4)+0.5,ccd1[n/4+n/4:n/2+n/4],rotation=90)

#pl.legend(loc = 'best')
pl.ylim(0,60000)
pl.ylabel('Fullwell (ADU)')
pl.grid()

pl.subplot(4,1,4)
pl.bar(np.arange(n/4),fw1[n/2+n/4:n/2+n/2],width=0.1,color='pink',label='12/7/2012')
pl.bar(np.arange(n/4)+0.1,fw2[n/2+n/4:n/2+n/2],width=0.1,color='green',label='7/26/2012')
pl.bar(np.arange(n/4)+0.2,fw3[n/2+n/4:n/2+n/2],width=0.1,color='blue',label='9/27/2012')
pl.bar(np.arange(n/4)+0.3,fw4[n/2+n/4:n/2+n/2],width=0.1,color='red',label='10/2/2012')
pl.bar(np.arange(n/4)+0.4,fw5[n/2+n/4:n/2+n/2],width=0.1,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/4)+0.5,fw7[n/2+n/4:n/2+n/2],width=0.1,color='black',label='11/9/2012')
pl.bar(np.arange(n/4)+0.6,fw8[n/2+n/4:n/2+n/2],width=0.1,color='c',label='11/19/2012')
pl.xticks(np.arange(n/4)+0.5,ccd1[n/2+n/4:n/2+n/2],rotation=90)

#pl.legend(loc = 'best')
pl.ylim(0,60000)
pl.ylabel('Fullwell (ADU)')
pl.grid()


pl.savefig('fullwell_comp_adu.png')

pl.close()

#-----compare the gain ---
"""
pl.figure(figsize=(20,10))
pl.subplot(2,1,1)
pl.bar(np.arange(n/2),1./gain1[0:n/2],width=0.15,color='pink',label='12/7/2012')
pl.bar(np.arange(n/2)+0.15,1./gain2[0:n/2],width=0.15,color='green',label='7/26/2012')
pl.bar(np.arange(n/2)+0.3,1./gain3[0:n/2],width=0.15,color='blue',label='9/27/2012')
pl.bar(np.arange(n/2)+0.45,1./gain4[0:n/2],width=0.15,color='red',label='10/2/2012')
pl.bar(np.arange(n/2)+0.6,1./gain5[0:n/2],width=0.15,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/2)+0.75,1./gain6[0:n/2],width=0.15,color='k',label='10/6/2012(short exp)')

pl.xticks(np.arange(n/2)+0.5,ccd1[0:n/2],rotation=90)
pl.legend(loc = 'best')
pl.ylim(0,10)
pl.grid()
pl.ylabel('Gain (e-/ADU)')

pl.subplot(2,1,2)
pl.bar(np.arange(n/2),1./gain1[n/2:],width=0.15,color='pink',label='12/7/2012')
pl.bar(np.arange(n/2)+0.15,1./gain2[n/2:],width=0.15,color='green',label='7/26/2012')
pl.bar(np.arange(n/2)+0.3,1./gain3[n/2:],width=0.15,color='blue',label='9/27/2012')
pl.bar(np.arange(n/2)+0.45,1./gain4[n/2:],width=0.15,color='red',label='10/2/2012')
pl.bar(np.arange(n/2)+0.6,1./gain5[n/2:],width=0.15,color='yellow',label='10/5/2012')
pl.bar(np.arange(n/2)+0.75,1./gain6[n/2:],width=0.15,color='k',label='10/6/2012(short exp)')

pl.xticks(np.arange(n/2)+0.5,ccd1[n/2:],rotation=90)
pl.legend(loc = 'best')
pl.ylim(0,10)
pl.ylabel('Gain (e-/ADU)')
pl.grid()
pl.savefig('gain_comp.png')
pl.close()
#---in adu---


#-----get average gain-----

gainall = np.array([1./gain1,1./gain2,1./gain3,1./gain4,1./gain5,1./gain6]).T
gain_mean = gainall.mean(axis=1)
gainstd = gainall.std(axis=1)/np.sqrt(6.)
np.savetxt('averageGain.txt',np.array([ccd1,gain_mean,gainstd]).T,fmt='%s',delimiter = ',')
p.dump(np.array([ccd1,gain_mean,gainstd]).T,open('averageGain.p','w'),2)
"""
