#! /usr/bin/env python
from DECamCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *

#import rpy2.robjects as robjects
#import rpy2.robjects.numpy2ri
#r=robjects.r

ccd=['s30','s29','s25','s26','s22','s21','s20','s14','s15','s16','s10','s9','s8','s1','s2','s3','s4','n4','n3','n2','n1']


fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))



baseDir = '/home/jghao/research/ccd/imager/flatness_7_29_11/upperright/'
cat=gl.glob(baseDir+'Image*catalog.fits')
cat.sort()
x=[]
y=[]
z=[]
zerr=[]
name=[]




#-------row 1---------
l0,lerr0=offset(CCD1='s4',CCD2='s30',cat=cat,xadd=260,yadd=30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s30',CCD2='s4',cat=cat,xadd=-260,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s30','s4')[0])
y.append(ccd_offset('s30','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s30')

l0,lerr0=offset(CCD1='s4',CCD2='s29',cat=cat,xadd=280,yadd=70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s29',CCD2='s4',cat=cat,xadd=-280,yadd=-70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s29','s4')[0])
y.append(ccd_offset('s29','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s29')


#----------row 2 ----------------

l0,lerr0=offset(CCD1='s4',CCD2='s26',cat=cat,xadd=160,yadd=50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s26',CCD2='s4',cat=cat,xadd=-160,yadd=-50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s26','s4')[0])
y.append(ccd_offset('s26','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s26')


l0,lerr0=offset(CCD1='s4',CCD2='s25',cat=cat,xadd=160,yadd=90,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s25',CCD2='s4',cat=cat,xadd=-160,yadd=-90,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s25','s4')[0])
y.append(ccd_offset('s25','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s25')

#----------row 3 ----------------

l0,lerr0=offset(CCD1='s4',CCD2='s22',cat=cat,xadd=120,yadd=20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s22',CCD2='s4',cat=cat,xadd=-120,yadd=-20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s22','s4')[0])
y.append(ccd_offset('s22','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s22')


l0,lerr0=offset(CCD1='s4',CCD2='s21',cat=cat,xadd=130,yadd=60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s21',CCD2='s4',cat=cat,xadd=-130,yadd=-60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s21','s4')[0])
y.append(ccd_offset('s21','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s21')

l0,lerr0=offset(CCD1='s4',CCD2='s20',cat=cat,xadd=130,yadd=110,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s20',CCD2='s4',cat=cat,xadd=-130,yadd=-110,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s20','s4')[0])
y.append(ccd_offset('s20','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s20')


#------row 4 --------------

l0,lerr0=offset(CCD1='s4',CCD2='s16',cat=cat,xadd=80,yadd=40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s16',CCD2='s4',cat=cat,xadd=-80,yadd=-40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s16','s4')[0])
y.append(ccd_offset('s16','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s16')

l0,lerr0=offset(CCD1='s4',CCD2='s15',cat=cat,xadd=90,yadd=80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s15',CCD2='s4',cat=cat,xadd=-90,yadd=-80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s15','s4')[0])
y.append(ccd_offset('s15','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s15')

l0,lerr0=offset(CCD1='s4',CCD2='s14',cat=cat,xadd=100,yadd=140,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s14',CCD2='s4',cat=cat,xadd=-100,yadd=-140,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s14','s4')[0])
y.append(ccd_offset('s14','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s14')


#-------row 5 ------------------

l0,lerr0=offset(CCD1='s4',CCD2='s10',cat=cat,xadd=50,yadd=30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s10',CCD2='s4',cat=cat,xadd=-50,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s10','s4')[0])
y.append(ccd_offset('s10','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s10')

l0,lerr0=offset(CCD1='s4',CCD2='s9',cat=cat,xadd=50,yadd=70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s9',CCD2='s4',cat=cat,xadd=-50,yadd=-70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s9','s4')[0])
y.append(ccd_offset('s9','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s9')

l0,lerr0=offset(CCD1='s4',CCD2='s8',cat=cat,xadd=70,yadd=130,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s8',CCD2='s4',cat=cat,xadd=-70,yadd=-130,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s8','s4')[0])
y.append(ccd_offset('s8','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s8')

#--------row 6 --------------------

x.append(0)
y.append(0)
z.append(0)
zerr.append(0)
name.append('s4')


l0,lerr0=offset(CCD1='s4',CCD2='s3',cat=cat,xadd=0,yadd=40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s3',CCD2='s4',cat=cat,xadd=0,yadd=-40,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s3','s4')[0])
y.append(ccd_offset('s3','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s3')

l0,lerr0=offset(CCD1='s4',CCD2='s2',cat=cat,xadd=20,yadd=100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s2',CCD2='s4',cat=cat,xadd=-20,yadd=-100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s2','s4')[0])
y.append(ccd_offset('s2','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s2')

l0,lerr0=offset(CCD1='s4',CCD2='s1',cat=cat,xadd=30,yadd=150,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='s1',CCD2='s4',cat=cat,xadd=-30,yadd=-150,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('s1','s4')[0])
y.append(ccd_offset('s1','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('s1')



#-----------row 7 ---------------
"""
x.append(0)
y.append(0)
z.append(0)
zerr.append(0)
name.append('s4')



l0,lerr0=offset(CCD1='s4',CCD2='n3',cat=cat,xadd=0,yadd= 30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n3',CCD2='s4',cat=cat,xadd=0,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n3','s4')[0])
y.append(ccd_offset('n3','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n3')


l0,lerr0=offset(CCD1='s4',CCD2='n2',cat=cat,xadd=-20,yadd=60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n2',CCD2='s4',cat=cat,xadd=20,yadd=-60,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n2','s4')[0])
y.append(ccd_offset('n2','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n2')

l0,lerr0=offset(CCD1='s4',CCD2='n1',cat=cat,xadd=-30,yadd=105,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n1',CCD2='s4',cat=cat,xadd=30,yadd=-105,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n1','s4')[0])
y.append(ccd_offset('n1','s4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n1')

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

p=np.array([0.0005774,-0.0030075,0.4109771])
zoffset=z-plane(p,x,y)

c1=pf.Column(name='x',format='D',array=x)
c2=pf.Column(name='y',format='D',array=y)
c3=pf.Column(name='z',format='D',array=z)
c4=pf.Column(name='zerr',format='D',array=zerr)
c5=pf.Column(name='zoffset',format='D',array=zoffset)
c6=pf.Column(name='CCD',format='10A',array=name)

hdu=pf.new_table([c1,c2,c3,c4,c5,c6])
os.system('rm '+baseDir+'flatness_upperright.fit')
hdu.writeto(baseDir+'flatness_upperright.fit')

data=zip(x,y,z)
data=np.array(data)
np.savetxt(baseDir+'flatness_upperright.txt',data)



"""
----analysis----
from fermiMCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *

baseDir = '/home/jghao/research/ccd/imager/flatness_7_29_11/upperright/'
data=pf.getdata(baseDir+'flatness_upperright.fit')
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
pl.savefig('/home/jghao/research/ccd/imager/flatness_7_29_11/results/upperright.png')

scale=2048.
xx=x/scale
yy=y/scale
zz=z/20.
zzoffset=zoffset/20.

barchart(xx,yy,zz)
barchart(xx,yy,zzoffset)


"""
