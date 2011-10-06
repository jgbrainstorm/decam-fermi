from fermiMCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *


ccd=['n4ul','n4ur','n4ll','n4lr']
fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))

cat=[]

b=gl.glob('/home/jghao/research/ccd/imager/flatness_7_8_11/upperleft/*i4n*.fits')
cat.append(b[0])

b=gl.glob('/home/jghao/research/ccd/imager/flatness_7_8_11/upperright/*i4n*.fits')
cat.append(b[0])

b=gl.glob('/home/jghao/research/ccd/imager/flatness_7_8_11/lowerleft/*i4n*.fits')
cat.append(b[0])

b=gl.glob('/home/jghao/research/ccd/imager/flatness_7_8_11/lowerright/*i4n*.fits')
cat.append(b[0])


offset_LR(CCD1='n4ul',CCD2='n4ur',cat=cat,xadd=0,yadd=0,sep=70,rot=0,crit_f=0,ccd2fileNo=ccd2fileNo)
 


#-------------UL vs UR -----------

l0,lerr0=offset_LR(CCD1='n4ul',CCD2='n4ur',cat=cat,xadd=70,yadd=-600,sep=70,rot=0,crit_f=0,ccd2fileNo=ccd2fileNo)

l1,lerr1=offset_LR(CCD1='n4ur',CCD2='n4ul',cat=cat,xadd=-60,yadd=600,sep=70,rot=0,crit_f=0,ccd2fileNo=ccd2fileNo)

zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.


#----------UL VS LL ---------

l0,lerr0=offset_LR(CCD1='n4ul',CCD2='n4ll',cat=cat,xadd=70,yadd=-620,sep=70,rot=180,crit_f=0,ccd2fileNo=ccd2fileNo)

l1,lerr1=offset_LR(CCD1='n4ll',CCD2='n4ul',cat=cat,xadd=-40,yadd=560,sep=70,rot=-180,crit_f=0,ccd2fileNo=ccd2fileNo)

zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.

#----------UL VS LR -------Not good--

l0,lerr0=offset_LR(CCD1='n4ul',CCD2='n4lr',cat=cat,xadd=-60,yadd=-690,sep=70,rot=180,crit_f=0,ccd2fileNo=ccd2fileNo)

l1,lerr1=offset_LR(CCD1='n4lr',CCD2='n4ul',cat=cat,xadd=70,yadd=600,sep=70,rot=-180,crit_f=0,ccd2fileNo=ccd2fileNo)

zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.



#----------LL VS LR ---------

l0,lerr0=offset_LR(CCD1='n4ll',CCD2='n4lr',cat=cat,xadd=-130,yadd=-50,sep=70,rot=0,crit_f=0,ccd2fileNo=ccd2fileNo)

l1,lerr1=offset_LR(CCD1='n4lr',CCD2='n4ll',cat=cat,xadd=130,yadd=50,sep=70,rot=0,crit_f=0,ccd2fileNo=ccd2fileNo)

zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.


#-------results------

UR -> UL: -11.903
LL -> UL: -136.914
LR -> LL: -79.436
LR -> UL: -136.914 - 79.436 = 216.350
