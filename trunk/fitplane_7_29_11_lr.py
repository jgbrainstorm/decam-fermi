#! /usr/bin/env python
from fermiMCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
r=robjects.r

ccd=['s1','s2','s3','s4','n4','n3','n2','n1','n8','n9','n10','n16','n15','n14','n20','n21','n22','n26','n25','n29','n30']



fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))



baseDir = '/home/jghao/research/ccd/imager/flatness_7_29_11/lowerright/'
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


l0,lerr0=offset(CCD1='n4',CCD2='n29',cat=cat,xadd=-180,yadd=20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n29',CCD2='n4',cat=cat,xadd=180,yadd=-20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n29','n4')[0])
y.append(ccd_offset('n29','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n29')


#----------row 2 ----------------

l0,lerr0=offset(CCD1='n4',CCD2='n26',cat=cat,xadd=-150,yadd=10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n26',CCD2='n4',cat=cat,xadd=150,yadd=-10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n26','n4')[0])
y.append(ccd_offset('n26','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n26')


l0,lerr0=offset(CCD1='n4',CCD2='n25',cat=cat,xadd=-130,yadd=50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n25',CCD2='n4',cat=cat,xadd=130,yadd=-50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n25','n4')[0])
y.append(ccd_offset('n25','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n25')

#----------row 3 ----------------

l0,lerr0=offset(CCD1='n4',CCD2='n22',cat=cat,xadd=-110,yadd=-10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n22',CCD2='n4',cat=cat,xadd=110,yadd=10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n22','n4')[0])
y.append(ccd_offset('n22','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n22')


l0,lerr0=offset(CCD1='n4',CCD2='n21',cat=cat,xadd=-100,yadd=30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n21',CCD2='n4',cat=cat,xadd=100,yadd=-30,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n21','n4')[0])
y.append(ccd_offset('n21','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n21')

l0,lerr0=offset(CCD1='n4',CCD2='n20',cat=cat,xadd=-100,yadd=70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n20',CCD2='n4',cat=cat,xadd=100,yadd=-70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n20','n4')[0])
y.append(ccd_offset('n20','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n20')


#------row 4 --------------

l0,lerr0=offset(CCD1='n4',CCD2='n16',cat=cat,xadd=-70,yadd=10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n16',CCD2='n4',cat=cat,xadd=70,yadd=-10,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n16','n4')[0])
y.append(ccd_offset('n16','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n16')

l0,lerr0=offset(CCD1='n4',CCD2='n15',cat=cat,xadd=-70,yadd=70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n15',CCD2='n4',cat=cat,xadd=70,yadd=-70,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n15','n4')[0])
y.append(ccd_offset('n15','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n15')

l0,lerr0=offset(CCD1='n4',CCD2='n14',cat=cat,xadd=-50,yadd=120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n14',CCD2='n4',cat=cat,xadd=50,yadd=-120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n14','n4')[0])
y.append(ccd_offset('n14','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n14')


#-------row 5 ------------------

l0,lerr0=offset(CCD1='n4',CCD2='n10',cat=cat,xadd=-30,yadd=20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n10',CCD2='n4',cat=cat,xadd=30,yadd=-20,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n10','n4')[0])
y.append(ccd_offset('n10','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n10')

l0,lerr0=offset(CCD1='n4',CCD2='n9',cat=cat,xadd=-20,yadd=80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n9',CCD2='n4',cat=cat,xadd=20,yadd=-80,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n9','n4')[0])
y.append(ccd_offset('n9','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n9')

l0,lerr0=offset(CCD1='n4',CCD2='n8',cat=cat,xadd=-20,yadd=120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n8',CCD2='n4',cat=cat,xadd=20,yadd=-120,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n8','n4')[0])
y.append(ccd_offset('n8','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n8')




#-----------row 6 ---------------

x.append(0)
y.append(0)
z.append(0)
zerr.append(0)
name.append('n4')



l0,lerr0=offset(CCD1='n4',CCD2='n3',cat=cat,xadd=20,yadd=50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n3',CCD2='n4',cat=cat,xadd=-20,yadd=-50,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n3','n4')[0])
y.append(ccd_offset('n3','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n3')


l0,lerr0=offset(CCD1='n4',CCD2='n2',cat=cat,xadd=10,yadd=100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n2',CCD2='n4',cat=cat,xadd=-10,yadd=-100,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n2','n4')[0])
y.append(ccd_offset('n2','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n2')

l0,lerr0=offset(CCD1='n4',CCD2='n1',cat=cat,xadd=20,yadd=150,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
l1,lerr1=offset(CCD1='n1',CCD2='n4',cat=cat,xadd=-20,yadd=-150,sep=70,crit_f=0,ccd2fileNo=ccd2fileNo,rot=0)
zm=(l0-l1)/2.
zmerr=abs(lerr0+lerr1)/2.
x.append(ccd_offset('n1','n4')[0])
y.append(ccd_offset('n1','n4')[1])
z.append(zm)
zerr.append(zmerr)
name.append('n1')


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


p=np.array([0.0007819,-0.0031909,-0.3362989])
zoffset=z-plane(p,x,y)

c1=pf.Column(name='x',format='D',array=x)
c2=pf.Column(name='y',format='D',array=y)
c3=pf.Column(name='z',format='D',array=z)
c4=pf.Column(name='zerr',format='D',array=zerr)
c5=pf.Column(name='zoffset',format='D',array=zoffset)
c6=pf.Column(name='CCD',format='10A',array=name)

hdu=pf.new_table([c1,c2,c3,c4,c5,c6])
os.system('rm '+baseDir+'flatness_lowerright.fit')
hdu.writeto(baseDir+'flatness_lowerright.fit')

data=zip(x,y,z)
data=np.array(data)
np.savetxt(baseDir+'flatness_lowerright.txt',data)



"""
----analysis----
from fermiMCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *

baseDir = '/home/jghao/research/ccd/imager/flatness_7_29_11/lowerright/'
data=pf.getdata(baseDir+'flatness_lowerright.fit')
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
pl.savefig('/home/jghao/research/ccd/imager/flatness_7_29_11/results/lowerright.png')


scale=2048.
xx=x/scale
yy=y/scale
zz=z/20.
zzoffset=zoffset/20.

barchart(xx,yy,zz)
barchart(xx,yy,zzoffset)


"""
