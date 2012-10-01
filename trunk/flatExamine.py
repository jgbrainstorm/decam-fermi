import sys
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi/pyRaytrace')
from decamRay import *

hdu = pf.open('masterFlat.fits')
data = []
ccd = []
gain2ap = np.genfromtxt('../fullwell.txt',dtype=None,names=['ccd','fw','gain'])['gain'][0:124]
gain=[]

for i in range(0,len(gain2ap),2):
    gain.append((gain2ap[i]+gain2ap[i+1])*0.5)

gain = np.array(gain)

for ext in range(1,63):
    x = eval(hdu[ext].header['detpos'])[1]
    y = eval(hdu[ext].header['detpos'])[2]
    col0=500
    col1=1500
    row0=500
    row1=3500
    mn = hdu[ext].data[row0:row1,col0:col1].mean()
    sd = hdu[ext].data[row0:row1,col0:col1].std()
    data.append([x,y,mn,sd])
    ccd.append(eval(hdu[ext].header['detpos'])[0])
data = np.array(data)
ccd = np.array(ccd)
srtidx = np.argsort(ccd)

pl.figure(figsize=(16,8))
pl.errorbar(np.arange(1,63),data[srtidx,2],yerr=data[srtidx,3],fmt='g.',label='flat ADU count')
pl.xticks(np.arange(1,63)+0.5,ccd[srtidx],rotation=-90)

pl.errorbar(np.arange(1,63),data[srtidx,2]/gain[srtidx],yerr=data[srtidx,3],fmt='b.',label='flat electron count')
pl.xticks(np.arange(1,63)+0.5,ccd[srtidx],rotation=-90)

pl.errorbar(np.arange(1,63),1./gain[srtidx],yerr=data[srtidx,3],fmt='r.',label='gain')
pl.xticks(np.arange(1,63)+0.5,ccd[srtidx],rotation=-90)
pl.legend(loc='best')
pl.savefig('flatfield.png')
pl.close()


pl.scatter(data[:,0],data[:,1],s=abs(data[:,2]-data[:,2].mean())*4000)
pl.xlabel('West')
pl.ylabel('North')
pl.grid()
pl.title('Variations of ADUs on')
pl.savefig('adu_variation.png')
pl.close()

pl.figure(figsize=(10,10))
pl.scatter(data[:,0],data[:,1],s=abs(1./gain-np.mean(1./gain))*1000,color='red')
pl.xlabel('West')
pl.ylabel('North')
pl.grid()
pl.title('Variations of Gains (e/adu)on masterFlat')
pl.savefig('gain_variation.png')
pl.close()

pl.figure(figsize=(10,10))
pl.scatter(data[:,0],data[:,1],s=abs(data[:,2]/gain-np.mean(data[:,2]/gain))*1000,color='green')
pl.xlabel('West')
pl.ylabel('North')
pl.grid()
pl.title('Variations of electrons on masterFlat')
pl.savefig('electron_variation.png')
pl.close()


betaMN,betaErr,R2_adj = zernikeFit(data[:,0],data[:,1],data[:,2]-data[:,2].mean(),max_order=10)

betaMN,betaErr,R2_adj = zernikeFit(data[:,0],data[:,1],data[:,2],max_order=40)


betaSD,betaErr,R2_adj = zernikeFit(data[:,0],data[:,1],data[:,3],max_order=20)

showZernike(betaMN)
