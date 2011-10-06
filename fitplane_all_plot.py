from fermiMCCD import *
from scipy.optimize import leastsq
from enthought.mayavi.mlab import *

#--------plot all ---------------


ul=pf.getdata('/home/jghao/research/ccd/imager/flatness_7_29_11/upperleft/flatness_upperleft.fit')
ur=pf.getdata('/home/jghao/research/ccd/imager/flatness_7_29_11/upperright/flatness_upperright.fit')

ll=pf.getdata('/home/jghao/research/ccd/imager/flatness_7_29_11/lowerleft/flatness_lowerleft.fit')
lr=pf.getdata('/home/jghao/research/ccd/imager/flatness_7_29_11/lowerright/flatness_lowerright.fit')

x=np.append(ul.field('x'),ur.field('x'))
x=np.append(x,ll.field('x'))
x=np.append(x,lr.field('x'))

y=np.append(ul.field('y'),ur.field('y'))
y=np.append(y,ll.field('y'))
y=np.append(y,lr.field('y'))

zoffset=np.append(ul.field('zoffset'),ur.field('zoffset'))
zoffset=np.append(zoffset,ll.field('zoffset'))
zoffset=np.append(zoffset,lr.field('zoffset'))


CCD=np.append(ul.field('CCD'),ur.field('CCD'))
CCD=np.append(CCD,ll.field('CCD'))
CCD=np.append(CCD,lr.field('CCD'))


scale=2048.
xx=x/scale
yy=y/scale
zzoffset=zoffset/20.
barchart(xx,yy,zzoffset)



pl.figure(figsize=(18,6))
pl.plot(np.arange(len(CCD)),zoffset,'b.')
pl.errorbar(np.arange(len(CCD)),zoffset,yerr=zerr,ecolor='b',fmt='.')
pl.xticks(np.arange(len(CCD)),CCD,rotation=90)
pl.hlines(np.median(zoffset),-1,len(x)+1,'r')
pl.hlines(np.median(zoffset)+2.1,-1,len(x)+1,'g')
pl.hlines(np.median(zoffset)-2.1,-1,len(x)+1,'g')
pl.ylabel('offset (pixels)')
pl.xlabel('CCDs')
pl.title('Offset w.r.t. the best fit plane of each quandrant')
pl.savefig('/home/jghao/research/ccd/imager/flatness_7_29_11/results/allccd.png')
