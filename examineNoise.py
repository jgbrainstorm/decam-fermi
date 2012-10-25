# this code examine the readout noise on the ccds close to the ion pump controller.
# J. Hao, 10/5/2012 @ CTIO
import sys
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
import numpy as np
import pyfits as pf
import pylab as pl
from DECamCCD import *
import scipy.fftpack as fft


def fftnoise(filename=None,ext='N4',ionpump='on'):
    if filename != None:
        b = pf.getdata(filename,ext)
        hdr = pf.getheader(filename,ext)
    elif ionpump == 'on':
        b=pf.getdata('/home3/data_local/images/fits/ptc_10_2_hao/bias/DECam_00137808.fits',ext)
        hdr=pf.getheader('/home3/data_local/images/fits/ptc_10_2_hao/bias/DECam_00137808.fits',ext)
    else:
        b = pf.getdata('/home3/data_local/images/fits/ptc_9_27_hao/bias/DECam_00136714.fits',ext)
        hdr = pf.getheader('/home3/data_local/images/fits/ptc_9_27_hao/bias/DECam_00136714.fits',ext)
    timestep = 4e-6# each pixel read out is 4 ms-> 250kpix/sec

    col0=int(hdr['BIASSECA'].split('[')[1].split(']')[0].split(',')[0].split(':')[0])-1
    col1=int(hdr['BIASSECA'].split('[')[1].split(']')[0].split(',')[0].split(':')[1])
    row0=int(hdr['BIASSECA'].split('[')[1].split(']')[0].split(',')[1].split(':')[0])-1
    row1=int(hdr['BIASSECA'].split('[')[1].split(']')[0].split(',')[1].split(':')[1])
    oscanA = b[row0:row1,col0:col1]
    oscanA1 = oscanA.reshape(oscanA.shape[0]*oscanA.shape[1])
    oscanA1FFT = np.fft.fft(oscanA1)
    col0=int(hdr['BIASSECB'].split('[')[1].split(']')[0].split(',')[0].split(':')[0])-1
    col1=int(hdr['BIASSECB'].split('[')[1].split(']')[0].split(',')[0].split(':')[1])
    row0=int(hdr['BIASSECB'].split('[')[1].split(']')[0].split(',')[1].split(':')[0])-1
    row1=int(hdr['BIASSECB'].split('[')[1].split(']')[0].split(',')[1].split(':')[1])
    oscanB = b[row0:row1,col0:col1]
    oscanB1 = oscanB.reshape(oscanB.shape[0]*oscanB.shape[1])
    oscanB1FFT = np.fft.fft(oscanB1) 
    freq = np.fft.fftfreq(len(oscanB1), d=timestep)
    ok = freq > 0
    pl.figure(figsize=(14,9))
    pl.subplot(2,1,1)
    pl.plot(np.abs(oscanA1FFT[ok]),'b-')
    pl.semilogy()
    xtickidx=np.arange(0,len(oscanA1),10000)[ok]
    pl.xticks(xtickidx,np.repeat('',len(xtickidx)))
    pl.grid()
    pl.title('Noise Spectra, ccd: '+ext)
    pl.ylabel('Amp: A')
    pl.subplot(2,1,2)
    pl.plot(np.abs(oscanB1FFT[ok]),'b-')
    pl.semilogy()
    pl.xticks(xtickidx,np.round(freq[xtickidx]),rotation=-60)
    pl.grid()
    pl.ylabel('Amp: B')
    pl.figtext(0.7,0.8,'Ion Pump: '+ionpump,color='red',fontsize=17)
    pl.savefig('noiseSpectra_'+ext+'ionpump_'+ionpump+'.png')
    return '---done!---'

if len(sys.argv) == 1:
    print 'syntax: examineNoise filename ext ionpumpStatus'
    print 'example: examineNoise bias.fits S2 off'
    
else:
    filename = sys.argv[1]
    ext = sys.argv[2]
    ionpump = sys.argv[3]
    t = fftnoise(filename=filename,ext=ext,ionpump=ionpump)
