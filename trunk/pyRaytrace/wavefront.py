from math import *
import numpy as np
import numpy.fft
import sys
from scipy.special import gamma
from numpy.random import normal
import scipy.fftpack as fft

def factorial(n):
  return gamma(n+1)

def zernike(j,npix=256,phase=0.0):
    if (j > 820):
      print "For n < 40, pick j < 820"
      sys.exit()
    x = np.arange(-npix/2,npix/2,dtype='d')
    y = np.arange(-npix/2,npix/2,dtype='d')

    xarr = np.outer(np.ones(npix,dtype='d'),x)
    yarr = np.outer(y,np.ones(npix,dtype='d'))

    rarr = np.sqrt(np.power(xarr,2) + np.power(yarr,2))/(npix/2)
    thetarr = np.arctan2(yarr,xarr) + phase

    outside = np.where(rarr > 1.0)

    narr = np.arange(40)
    jmax = (narr+1)*(narr+2)/2
    wh = np.where(j <= jmax)
    n = wh[0][0]
    mprime = j - n*(n+1)/2
    if ((n % 2) == 0):
      m = 2*int(floor(mprime/2))
    else:
      m = 1 + 2*int(floor((mprime-1)/2))

    radial = np.zeros((npix,npix),dtype='d')
    zarr = np.zeros((npix,npix),dtype='d')

    for s in range((n-m)/2 + 1):
      tmp = pow(-1,s) * factorial(n-s)
      tmp /= factorial(s)*factorial((n+m)/2 - s)*factorial((n-m)/2 - s)
      radial += tmp*np.power(rarr,n-2*s)

    if (m == 0):
      zarr = radial
    else:
      if ((j % 2) == 0):
        zarr = sqrt(2.0)*radial*np.cos(m*thetarr)
      else:
        zarr = sqrt(2.0)*radial*np.sin(m*thetarr)

    zarr *= sqrt(n+1)
    zarr[outside] = 0.0

    return zarr

def aperture(npix=256, cent_obs=0.0, spider=0):

  illum = np.ones((npix,npix),dtype='d')
  x = np.arange(-npix/2,npix/2,dtype='d')
  y = np.arange(-npix/2,npix/2,dtype='d')

  xarr = np.outer(np.ones(npix,dtype='d'),x)
  yarr = np.outer(y,np.ones(npix,dtype='d'))

  rarr = np.sqrt(np.power(xarr,2) + np.power(yarr,2))/(npix/2)
  outside = np.where(rarr > 1.0)
  inside = np.where(rarr < cent_obs)

  illum[outside] = 0.0
  if np.any(inside[0]):
    illum[inside] = 0.0

  if (spider > 0):
   start = npix/2 - int(spider)/2
   illum[start:start+int(spider),:] = 0.0
   illum[:,start:start+int(spider)] = 0.0

  return illum

def plane_wave(npix=256):
  wf = np.zeros((npix,npix),dtype='d')
  return wf

def seeing(d_over_r0, npix=256, nterms=15, level=None, quiet=False):
  scale = pow(d_over_r0,5.0/3.0)
  if (nterms < 10):
    print "C'mon, at least use ten terms..."
    sys.exit()
  if level:
    narr = np.arange(400,dtype='d') + 2
    coef = np.sqrt(0.2944*scale*(np.power((narr-1),-0.866) - np.power(narr,-0.866)))
    wh = np.where(coef < level)
    n = wh[0][0]
    norder = int(ceil(sqrt(2*n)-0.5))
    nterms = norder*(norder+1)/2
    if (nterms < 15):
      nterms = 15
  wf = np.zeros((npix,npix),dtype='d')
  resid = np.zeros(nterms,dtype='d')
  coeff = np.zeros(nterms,dtype='d')
  resid[0:10] = [1.030,0.582,0.134,0.111,0.088,0.065,0.059,0.053,0.046,0.040]
  if (nterms > 10):
    for i in range(10,nterms):
      resid[i] = 0.2944*pow(i+1,-0.866)
  for j in range(2,nterms+1):
    coeff[j-1] = sqrt((resid[j-2]-resid[j-1])*scale)
    wf += coeff[j-1]*normal()*zernike(j,npix=npix)
  if not quiet:
    print "Computed Zernikes to term %d and RMS %f" % (nterms,coeff[nterms-1])
  return wf

def psf(aperture, wavefront, overfill=1):
  npix = len(wavefront)
  nbig = npix*overfill
  wfbig = np.zeros((nbig,nbig),dtype='d')
  half = (nbig - npix)/2
  wfbig[half:half+npix,half:half+npix] = wavefront
  illum = np.zeros((nbig,nbig),dtype='d')
  illum[half:half+npix,half:half+npix] = aperture
  phase = np.exp(wfbig*(0.+1.j))
  input = illum*phase
  ft = fft.fft2(input)
  ft = fft.fftshift(ft)
  powft = abs(ft)**2
  crop =  powft[half:half+npix,half:half+npix]
  #fluxrat = np.sum(crop)/np.sum(sorted)
  #print "Cropped PSF has %.2f%% of the flux" % (100*fluxrat)
  return crop

def flux_in(img,ctrx,ctry,rad):
  flux = 0.0
  xp = np.outer(np.ones(10,dtype='d'),np.linspace(-0.45,0.45,10))
  yp = np.outer(np.linspace(-0.45,0.45,10),np.ones(10,dtype='d'))
  for x in range(np.size(img,0)):
    for y in range(np.size(img,1)):
      r = sqrt(pow(x-ctrx,2) + pow(y-ctry,2))
      if (r - rad < 1.0):
        xgrid = x + xp - ctrx
        ygrid = y + yp - ctry
        rgrid = np.sqrt(xgrid*xgrid + ygrid*ygrid)
        whin = np.where(rgrid < rad)
        count = len(whin[0])
        flux += img[x][y]*count/100.0
  return flux


def overlapIndices(a1, a2, 
                   shiftx, shifty):
    if shiftx >=0:
        a1xbeg=shiftx
        a2xbeg=0
        a1xend=a1.shape[0]
        a2xend=a1.shape[0]-shiftx
    else:
        a1xbeg=0
        a2xbeg=-shiftx
        a1xend=a1.shape[0]+shiftx
        a2xend=a1.shape[0]

    if shifty >=0:
        a1ybeg=shifty
        a2ybeg=0
        a1yend=a1.shape[1]
        a2yend=a1.shape[1]-shifty
    else:
        a1ybeg=0
        a2ybeg=-shifty
        a1yend=a1.shape[1]+shifty
        a2yend=a1.shape[1]

    return (a1xbeg, a1xend, a1ybeg, a1yend), (a2xbeg, a2xend, a2ybeg, a2yend)



def hogbom(dirty,
           psf,
           window,
           gain,
           thresh,
           niter):
    """
    Hogbom clean

    :param dirty: The dirty image, i.e., the image to be deconvolved

    :param psf: The point spread-function

    :param window: Regions where clean components are allowed. If
    True, thank all of the dirty image is assumed to be allowed for
    clean components

    :param gain: The "loop gain", i.e., the fraction of the brightest
    pixel that is removed in each iteration

    :param thresh: Cleaning stops when the maximum of the absolute
    deviation of the residual is less than this value

    :param niter: Maximum number of components to make if the
    threshold "thresh" is not hit
    """
    comps=np.zeros(dirty.shape)
    res=np.array(dirty)
    if window is True:
        window=np.ones(dirty.shape,
                          np.bool)
    for i in range(niter):
        mx, my=np.unravel_index(np.fabs(res[window]).argmax(), res.shape)
        mval=res[mx, my]*gain
        comps[mx, my]+=mval
        a1o, a2o=overlapIndices(dirty, psf,
                                mx-dirty.shape[0]/2,
                                my-dirty.shape[1]/2)
        res[a1o[0]:a1o[1],a1o[2]:a1o[3]]-=psf[a2o[0]:a2o[1],a2o[2]:a2o[3]]*mval
        if np.fabs(res).max() < thresh:
            break
    return comps, res
