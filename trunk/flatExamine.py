import sys
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi/pyRaytrace')
from decamRay import *
import pickle as p

hdu = pf.open('masterFlat.fits')
data = []
ccd = np.genfromtxt('../fullwell.txt',dtype=None,names=['ccd','fw','gain'])['ccd'][0:124]
gain = np.genfromtxt('../fullwell.txt',dtype=None,names=['ccd','fw','gain'])['gain'][0:124]

gain = gain[0:124]

for ext in range(1,63):
    xhigh = eval(hdu[ext].header['detpos'])[1]+15.
    xlow = eval(hdu[ext].header['detpos'])[1]-15.
    y = eval(hdu[ext].header['detpos'])[2]
    detector = hdu[ext].header['DETPOS']   
    if detector[0]=='S':
        # --- left---:
        colmin=1360
        colmax=1860
        rowmin=500  
        rowmax=3500
        mn = hdu[ext].data[rowmin:rowmax,colmin:colmax].mean()
        sd = hdu[ext].data[rowmin:rowmax,colmin:colmax].std()
        mne = mn/gain[2*(ext-1)]
        sde = sd/gain[2*(ext-1)]
        x = xlow
        data.append([x,y,mn,mne,sd,sde,gain[2*(ext-1)]])
        # ---right ---:
        colmin=300
        colmax=800
        rowmin=500
        rowmax=3500
        mn = hdu[ext].data[rowmin:rowmax,colmin:colmax].mean()
        sd = hdu[ext].data[rowmin:rowmax,colmin:colmax].std()
        mne = mn/gain[2*(ext-1)+1]
        sde = sd/gain[2*(ext-1)+1]
        x = xhigh
        data.append([x,y,mn,mne,sd,sde,gain[2*(ext-1)+1]])
 
    if detector[0]=='N':
        # ---left ---:
        colmin=300
        colmax=800
        rowmin=500
        rowmax=3500
        mn = hdu[ext].data[rowmin:rowmax,colmin:colmax].mean()
        sd = hdu[ext].data[rowmin:rowmax,colmin:colmax].std()
        mne = mn/gain[2*(ext-1)]
        sde = sd/gain[2*(ext-1)]
        x = xhigh
        data.append([x,y,mn,mne,sd,sde,gain[2*(ext-1)]])

        # --- right ---:
        colmin=1360
        colmax=1860
        rowmin=500
        rowmax=3500
        mn = hdu[ext].data[rowmin:rowmax,colmin:colmax].mean()
        sd = hdu[ext].data[rowmin:rowmax,colmin:colmax].std()
        mne = mn/gain[2*(ext-1)+1]
        sde = sd/gain[2*(ext-1)+1]
        x = xlow
        data.append([x,y,mn,mne,sd,sde,gain[2*(ext-1)+1]])


data = np.array(data)
srtidx = np.argsort(ccd)

data = data[srtidx,:]
np.savetxt('flatField.txt',data,fmt='%10.5f',delimiter=',')

pl.figure(figsize=(16,8))
pl.subplot(2,1,1)
pl.bar(np.arange(0,62),data[0:62,3],yerr=data[0:62,5],color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[0:62],rotation=-90)
pl.title('Photon Count of the Flat Field')

pl.subplot(2,1,2)
pl.bar(np.arange(0,62),data[62:125,3],yerr=data[62:125,5],color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[62:125],rotation=-90)
pl.savefig('photon_count_flat.png')

pl.figure(figsize=(16,8))
pl.subplot(2,1,1)
pl.bar(np.arange(0,62),(data[0:62,3] - data[0:62,3].mean())/data[0:62,3].mean(),color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[0:62],rotation=-90)
pl.title('Photon Count fractional variation of the Flat Field')

pl.subplot(2,1,2)
pl.bar(np.arange(0,62),(data[62:125,3] - data[62:125,3].mean())/data[62:125,3].mean(),color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[62:125],rotation=-90)
pl.savefig('fractional_photon_count_flat.png')

