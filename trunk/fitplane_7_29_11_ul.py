#! /usr/bin/env python
from DECamCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *

#import rpy2.robjects as robjects
#import rpy2.robjects.numpy2ri
#r=robjects.r

ccd=['s30','s31','s28','s27','s22','s23','s24','s19','s18','s17','s11','s12','s13','s7','s6','s5','s4','n4','n5','n6','n7']


fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))



baseDir = '/home/jghao/research/ccd/imager/flatness_7_29_11/upperleft/'
cat=gl.glob(baseDir+'Image*catalog.fits')
cat.sort()
x=[]
y=[]
z=[]
zerr=[]
name=[]





#-------row 1---------
l0,lerr0=offset(CCD1='s4',CCD2='s30',cat=cat,xadd=200,yadd=30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s30',CCD2='s4',cat=cat,xadd=-200,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s30','s4')[0])
y.append(ccd_offset('s30','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s30')

l0,lerr0=offset(CCD1='s4',CCD2='s31',cat=cat,xadd=170,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s31',CCD2='s4',cat=cat,xadd=-170,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s31','s4')[0])
y.append(ccd_offset('s31','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s31')


#----------row 2 ----------------

l0,lerr0=offset(CCD1='s4',CCD2='s27',cat=cat,xadd=150,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s27',CCD2='s4',cat=cat,xadd=-150,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s27','s4')[0])
y.append(ccd_offset('s27','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s27')


l0,lerr0=offset(CCD1='s4',CCD2='s28',cat=cat,xadd=130,yadd=-50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s28',CCD2='s4',cat=cat,xadd=-130,yadd=50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s28','s4')[0])
y.append(ccd_offset('s28','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s28')

#----------row 3 ----------------

l0,lerr0=offset(CCD1='s4',CCD2='s22',cat=cat,xadd=120,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s22',CCD2='s4',cat=cat,xadd=-120,yadd=0,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s22','s4')[0])
y.append(ccd_offset('s22','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s22')


l0,lerr0=offset(CCD1='s4',CCD2='s23',cat=cat,xadd=100,yadd=-40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s23',CCD2='s4',cat=cat,xadd=-100,yadd=40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s23','s4')[0])
y.append(ccd_offset('s23','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s23')

l0,lerr0=offset(CCD1='s4',CCD2='s24',cat=cat,xadd=90,yadd=-80,sep=70,crit_f=15,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s24',CCD2='s4',cat=cat,xadd=-90,yadd=80,sep=70,crit_f=15,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s24','s4')[0])
y.append(ccd_offset('s24','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s24')


#------row 4 --------------

l0,lerr0=offset(CCD1='s4',CCD2='s17',cat=cat,xadd=70,yadd=-10,sep=70,crit_f=5,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s17',CCD2='s4',cat=cat,xadd=-70,yadd=10,sep=70,crit_f=5,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s17','s4')[0])
y.append(ccd_offset('s17','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s17')

l0,lerr0=offset(CCD1='s4',CCD2='s18',cat=cat,xadd=70,yadd=-60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s18',CCD2='s4',cat=cat,xadd=-70,yadd=60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s18','s4')[0])
y.append(ccd_offset('s18','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s18')

l0,lerr0=offset(CCD1='s4',CCD2='s19',cat=cat,xadd=50,yadd=-110,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s19',CCD2='s4',cat=cat,xadd=-50,yadd=110,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s19','s4')[0])
y.append(ccd_offset('s19','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s19')


#-------row 5 ------------------

l0,lerr0=offset(CCD1='s4',CCD2='s11',cat=cat,xadd=30,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s11',CCD2='s4',cat=cat,xadd=-30,yadd=30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s11','s4')[0])
y.append(ccd_offset('s11','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s11')

l0,lerr0=offset(CCD1='s4',CCD2='s12',cat=cat,xadd=30,yadd=-70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s12',CCD2='s4',cat=cat,xadd=-30,yadd=70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s12','s4')[0])
y.append(ccd_offset('s12','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s12')

l0,lerr0=offset(CCD1='s4',CCD2='s13',cat=cat,xadd=10,yadd=-120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s13',CCD2='s4',cat=cat,xadd=-10,yadd=120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s13','s4')[0])
y.append(ccd_offset('s13','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s13')

#--------row 6 --------------------

x.append(0)
y.append(0)
z.append(0)
zerr.append(0)
name.append('s4')


l0,lerr0=offset(CCD1='s4',CCD2='s5',cat=cat,xadd=-10,yadd=-40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s5',CCD2='s4',cat=cat,xadd=10,yadd=40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s5','s4')[0])
y.append(ccd_offset('s5','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s5')

l0,lerr0=offset(CCD1='s4',CCD2='s6',cat=cat,xadd=-20,yadd=-100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s6',CCD2='s4',cat=cat,xadd=20,yadd=100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s6','s4')[0])
y.append(ccd_offset('s6','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s6')

l0,lerr0=offset(CCD1='s4',CCD2='s7',cat=cat,xadd=-30,yadd=-150,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s7',CCD2='s4',cat=cat,xadd=30,yadd=150,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s7','s4')[0])
y.append(ccd_offset('s7','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s7')



#-----------row 7 ---------------
"""
x.append(0)
y.append(0)
z.append(0)
zerr.append(0)
name.append('s4')



l0,lerr0=offset(CCD1='s4',CCD2='n5',cat=cat,xadd=-60,yadd=-20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n5',CCD2='s4',cat=cat,xadd=60,yadd=20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n5','s4')[0])
y.append(ccd_offset('n5','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n5')


l0,lerr0=offset(CCD1='s4',CCD2='n6',cat=cat,xadd=-60,yadd=-60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n6',CCD2='s4',cat=cat,xadd=60,yadd=60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n6','s4')[0])
y.append(ccd_offset('n6','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n6')

l0,lerr0=offset(CCD1='s4',CCD2='n7',cat=cat,xadd=-60,yadd=-120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n7',CCD2='s4',cat=cat,xadd=60,yadd=120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n7','s4')[0])
y.append(ccd_offset('n7','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n7')


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


p=np.array([0.0004993,-0.0033122,0.1510956])

zoffset=z-plane(p,x,y)

c1=pf.Column(name='x',format='D',array=x)
c2=pf.Column(name='y',format='D',array=y)
c3=pf.Column(name='z',format='D',array=z)
c4=pf.Column(name='zerr',format='D',array=zerr)
c5=pf.Column(name='zoffset',format='D',array=zoffset)
c6=pf.Column(name='CCD',format='10A',array=name)

hdu=pf.new_table([c1,c2,c3,c4,c5,c6])
os.system('rm '+baseDir+'flatness_upperleft.fit')
hdu.writeto(baseDir+'flatness_upperleft.fit')

data=zip(x,y,z)
data=np.array(data)
np.savetxt(baseDir+'flatness_upperleft.txt',data)



"""
----analysis----
from fermiMCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *

baseDir = '/home/jghao/research/ccd/imager/flatness_7_29_11/upperleft/'
data=pf.getdata(baseDir+'flatness_upperleft.fit')
x=data.field('x')
y=data.field('y')
z=data.field('z')
zerr=data.field('zerr')
zoffset=data.field('zoffset')
name=data.field('CCD')

pl.figure(figsize=(11,6))
pl.plot(np.arange(len(name)),zoffset,'b.')
pl.errorbar(np.arange(len(name)),zoffset,yerr=zerr,ecolor='b',fmt='.')
pl.xticks(np.arange(len(name)),name)
pl.hlines(np.median(zoffset),-1,len(x)+1,'r')
pl.hlines(np.median(zoffset)+2.1,-1,len(x)+1,'g')
pl.hlines(np.median(zoffset)-2.1,-1,len(x)+1,'g')
pl.ylabel('offset (pixels)')
pl.xlabel('CCDs')
pl.savefig('/home/jghao/research/ccd/imager/flatness_7_29_11/results/upperleft.png')

scale=2048.
xx=x/scale
yy=y/scale
zz=z/20.
zzoffset=zoffset/20.

barchart(xx[0:-4],yy[0:-4],zz[0:-4])
barchart(xx,yy,zzoffset)


"""
