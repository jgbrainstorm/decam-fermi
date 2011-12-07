#! /usr/bin/env python
from DECamCCD import *
#from scipy.optimize import leastsq
#from enthought.mayavi.mlab import *

#import rpy2.robjects as robjects
#import rpy2.robjects.numpy2ri
#r=robjects.r

ccd=['s7','s6','s5','s4','n4','n5','n6','n7','n13','n12','n11','n17','n18','n19','n24','n23','n22','n27','n28','n31','n30']


fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))


baseDir='/usr/remote/user/sispi/jiangang/test_flat/lowerleft/'
#baseDir = '/home/jghao/research/ccd/imager/flatness_7_29_11/lowerleft/'
cat=gl.glob(baseDir+'Image*catalog.fits')
cat.sort()
x=[]
y=[]
z=[]
zerr=[]
name=[]




#-------row 1---------
l0,lerr0=offset(CCD1='n4',CCD2='n30',cat=cat,xadd=-200,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n30',CCD2='n4',cat=cat,xadd=200,yadd=30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n30','n4')[0])
y.append(ccd_offset('n30','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n30')

l0,lerr0=offset(CCD1='n4',CCD2='n31',cat=cat,xadd=-200,yadd=-80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n31',CCD2='n4',cat=cat,xadd=200,yadd=80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n31','n4')[0])
y.append(ccd_offset('n31','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n31')


#----------row 2 ----------------

l0,lerr0=offset(CCD1='n4',CCD2='n27',cat=cat,xadd=-150,yadd=-40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n27',CCD2='n4',cat=cat,xadd=150,yadd=40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n27','n4')[0])
y.append(ccd_offset('n27','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n27')


l0,lerr0=offset(CCD1='n4',CCD2='n28',cat=cat,xadd=-165,yadd=-100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n28',CCD2='n4',cat=cat,xadd=165,yadd=100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n28','n4')[0])
y.append(ccd_offset('n28','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n28')

#----------row 3 ----------------

l0,lerr0=offset(CCD1='n4',CCD2='n22',cat=cat,xadd=-80,yadd=-10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n22',CCD2='n4',cat=cat,xadd=80,yadd=10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n22','n4')[0])
y.append(ccd_offset('n22','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n22')


l0,lerr0=offset(CCD1='n4',CCD2='n23',cat=cat,xadd=-120,yadd=-70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n23',CCD2='n4',cat=cat,xadd=120,yadd=70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n23','n4')[0])
y.append(ccd_offset('n23','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n23')

l0,lerr0=offset(CCD1='n4',CCD2='n24',cat=cat,xadd=-140,yadd=-120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n24',CCD2='n4',cat=cat,xadd=140,yadd=120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n24','n4')[0])
y.append(ccd_offset('n24','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n24')


#------row 4 --------------

l0,lerr0=offset(CCD1='n4',CCD2='n17',cat=cat,xadd=-70,yadd=-40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n17',CCD2='n4',cat=cat,xadd=70,yadd=40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n17','n4')[0])
y.append(ccd_offset('n17','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n17')

l0,lerr0=offset(CCD1='n4',CCD2='n18',cat=cat,xadd=-90,yadd=-80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n18',CCD2='n4',cat=cat,xadd=90,yadd=80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n18','n4')[0])
y.append(ccd_offset('n18','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n18')

l0,lerr0=offset(CCD1='n4',CCD2='n19',cat=cat,xadd=-100,yadd=-140,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n19',CCD2='n4',cat=cat,xadd=100,yadd=140,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n19','n4')[0])
y.append(ccd_offset('n19','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n19')


#-------row 5 ------------------

l0,lerr0=offset(CCD1='n4',CCD2='n11',cat=cat,xadd=-50,yadd=-20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n11',CCD2='n4',cat=cat,xadd=50,yadd=20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n11','n4')[0])
y.append(ccd_offset('n11','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n11')

l0,lerr0=offset(CCD1='n4',CCD2='n12',cat=cat,xadd=-50,yadd=-80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n12',CCD2='n4',cat=cat,xadd=50,yadd=80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n12','n4')[0])
y.append(ccd_offset('n12','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n12')

l0,lerr0=offset(CCD1='n4',CCD2='n13',cat=cat,xadd=-70,yadd=-140,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n13',CCD2='n4',cat=cat,xadd=70,yadd=140,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n13','n4')[0])
y.append(ccd_offset('n13','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n13')




#-----------row 6 ---------------

x.append(0)
y.append(0)
z.append(0)
zerr.append(0)
name.append('n4')



l0,lerr0=offset(CCD1='n4',CCD2='n5',cat=cat,xadd=0,yadd=-60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n5',CCD2='n4',cat=cat,xadd=0,yadd=60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n5','n4')[0])
y.append(ccd_offset('n5','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n5')


l0,lerr0=offset(CCD1='n4',CCD2='n6',cat=cat,xadd=-20,yadd=-110,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n6',CCD2='n4',cat=cat,xadd=20,yadd=110,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n6','n4')[0])
y.append(ccd_offset('n6','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n6')

l0,lerr0=offset(CCD1='n4',CCD2='n7',cat=cat,xadd=-30,yadd=-160,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n7',CCD2='n4',cat=cat,xadd=30,yadd=160,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n7','n4')[0])
y.append(ccd_offset('n7','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n7')


#--------row 7 --------------------
"""
l0,lerr0=offset(CCD1='n4',CCD2='s4',cat=cat,xadd=50,yadd=30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s4',CCD2='n4',cat=cat,xadd=-50,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s4','n4')[0])
y.append(ccd_offset('s4','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s4')

l0,lerr0=offset(CCD1='n4',CCD2='s5',cat=cat,xadd=50,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s5',CCD2='n4',cat=cat,xadd=-50,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s5','n4')[0])
y.append(ccd_offset('s5','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s5')

l0,lerr0=offset(CCD1='n4',CCD2='s6',cat=cat,xadd=70,yadd=-20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s6',CCD2='n4',cat=cat,xadd=-70,yadd=20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s6','n4')[0])
y.append(ccd_offset('s6','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s6')

l0,lerr0=offset(CCD1='n4',CCD2='s7',cat=cat,xadd=70,yadd=-50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s7',CCD2='n4',cat=cat,xadd=-70,yadd=50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s7','n4')[0])
y.append(ccd_offset('s7','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s7')

"""




#---------------





x=np.array(x)
y=np.array(y)
z=np.array(z)
zerr=np.array(zerr)

def bfplane(p,x,y,z,zerr):
    A=p[0]
    B=p[1]
    C=p[2]
    res=((A*x+B*y+C)-z)/zerr
    return res

def plane(p,x,y):
    A=p[0]
    B=p[1]
    C=p[2]
    res=(A*x+B*y+C)
    return res

p=np.array([0.0007606,-0.0034225,-0.6859953])
zoffset=z-plane(p,x,y)

c1=pf.Column(name='x',format='D',array=x)
c2=pf.Column(name='y',format='D',array=y)
c3=pf.Column(name='z',format='D',array=z)
c4=pf.Column(name='zerr',format='D',array=zerr)
c5=pf.Column(name='zoffset',format='D',array=zoffset)
c6=pf.Column(name='CCD',format='10A',array=name)

hdu=pf.new_table([c1,c2,c3,c4,c5,c6])
os.system('rm '+baseDir+'flatness_lowerleft.fit')
hdu.writeto(baseDir+'flatness_lowerleft.fit')

data=zip(x,y,z)
data=np.array(data)
np.savetxt(baseDir+'flatness_lowerleft.txt',data)



