#! /usr/bin/env python

import sys
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi/pyRaytrace')
from decamRay import *
import pickle as p

hdu = pf.open('/data_local/images/fits/masterFlat/masterFlat.fits')
data = []
ccd = np.genfromtxt('../fullwell.txt',dtype=None,names=['ccd','fw','gain'])['ccd'][0:124]
gain = np.genfromtxt('../fullwell.txt',dtype=None,names=['ccd','fw','gain'])['gain'][0:124]

gain = gain[0:124]

for ext in range(1,63):
    xhigh = eval(hdu[ext].header['detpos'])[1]+7.5
    xlow = eval(hdu[ext].header['detpos'])[1]-7.5
    y = eval(hdu[ext].header['detpos'])[2]
    detector = hdu[ext].header['DETPOS']   
    if detector[0]=='S':
        # --- left---:
        colmin=1360
        colmax=1860
        rowmin=500  
        rowmax=3500
        mn = np.median(hdu[ext].data[rowmin:rowmax,colmin:colmax])
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
        mn = np.median(hdu[ext].data[rowmin:rowmax,colmin:colmax])
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
        mn = np.median(hdu[ext].data[rowmin:rowmax,colmin:colmax])
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
        mn = np.median(hdu[ext].data[rowmin:rowmax,colmin:colmax])
        sd = hdu[ext].data[rowmin:rowmax,colmin:colmax].std()
        mne = mn/gain[2*(ext-1)+1]
        sde = sd/gain[2*(ext-1)+1]
        x = xlow
        data.append([x,y,mn,mne,sd,sde,gain[2*(ext-1)+1]])


data = np.array(data)

np.savetxt('flatField.txt',data,fmt='%10.5f',delimiter=',')

# ---photon count variation ----

pl.figure(figsize=(16,8))
pl.subplot(2,1,1)
pl.bar(np.arange(0,62),data[0:62,3],yerr=data[0:62,5],color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[0:62],rotation=-90)
pl.title('Photon Count of the Flat Field')

pl.subplot(2,1,2)
pl.bar(np.arange(0,62),data[62:125,3],yerr=data[62:125,5],color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[62:125],rotation=-90)
pl.savefig('photon_count_flat.png')

#-----fractional change ---
pl.figure(figsize=(16,8))
pl.subplot(2,1,1)
pl.bar(np.arange(0,62),(data[0:62,3] - data[0:62,3].mean())/data[0:62,3].mean(),color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[0:62],rotation=-90)
pl.title('Photon Count fractional variation of the Flat Field')

pl.subplot(2,1,2)
pl.bar(np.arange(0,62),(data[62:125,3] - data[62:125,3].mean())/data[62:125,3].mean(),color=['red','green'])
pl.xticks(np.arange(0,62)+0.5,ccd[62:125],rotation=-90)
pl.savefig('fractional_photon_count_flat.png')


#---make new plot based on the new flat -------------
pl.figure(figsize=(10,10))
b=np.genfromtxt('flatField_averageGain.txt',delimiter=',')
pl.scatter(b[:,0],b[:,1],s=200*(b[:,3]-b[:,3].mean()),color='red',label='positive')
pl.scatter(b[:,0],b[:,1],s=-200*(b[:,3]-b[:,3].mean()),color='blue',label='negative')
pl.legend(loc='best')
pl.title('photon count - mean photon count')
pl.grid()
pl.savefig('photon_distribution.png')
pl.close()

pl.figure(figsize=(10,10))
beta,betaErr,R2 = zernikeFit(b[:,0],b[:,1],b[:,3]-b[:,3].mean(),max_order=15)
showZernike(beta)
pl.grid()
pl.savefig('photon_distribution_zernike.png')
pl.close()
