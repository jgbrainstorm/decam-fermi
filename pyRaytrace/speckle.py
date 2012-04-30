#!/usr/bin/env python

from math import *
import numpy
import wavefront
import Image
import pyfits
import os

N = 256						# number of pixels across
D_over_r0 = 20.0				# seeing is  lambda/r0 radians

wf = wavefront.seeing(D_over_r0,npix=N,level=0.05)	# get wavefront

wfmax = numpy.max(wf)
wfmin = numpy.min(wf)
wfrng = wfmax - wfmin

print wfmin,wfmax				# reported in radians

# create GIF image of wavefront
im = Image.new('L',(N,N))
for i in range(N):
  for j in range(N):
    im.putpixel((i,j),127.0 + wf[i][j]*256/wfrng)
im.save('wavefront.gif')

# establish aperture illumination pattern
illum = wavefront.aperture(npix=N, cent_obs = 0.3,spider=2)
#illum = wavefront.aperture(npix=N, cent_obs = 0.0,spider=0)

# create GIF image of aperture illumination
im = Image.new('L',(N,N))
for i in range(N):
  for j in range(N):
    im.putpixel((i,j),illum[i][j]*255)
im.save('aperture.gif')

psf_scale = 3	# pads aperture by this amount; resultant pix scale is
		# lambda/D/psf_scale, so for instance full frame 256 pix
		# for 3.5 m at 532 nm is 256*5.32e-7/3.5/3 = 2.67 arcsec
		# for psf_scale = 3

# generate speckle pattern given my wavefront and aperture illumination
psf = wavefront.psf(illum,wf,overfill=psf_scale)

psfmax = numpy.max(psf)
psfmin = numpy.min(psf)
psfrng = psfmax - psfmin

# create GIF image of  speckle pattern
im = Image.new('L',(N,N))
for i in range(N):
  for j in range(N):
    im.putpixel((i,j),63.0 + psf[i][j]*192/psfrng)
im.save('psf.gif')

# create FITS image of speckle pattern
if 'psf.fits' in os.listdir('.'):
  os.unlink('psf.fits')
hdu = pyfits.PrimaryHDU(psf)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('psf.fits')

# make a plane wavefront and create telescope diffraction-limited pattern
pwf = wavefront.plane_wave(npix=N)
diffrac = wavefront.psf(illum,pwf,overfill=psf_scale)

diffracmax = numpy.max(diffrac)
diffracmin = numpy.min(diffrac)
diffracrng = diffracmax - diffracmin

# create GIF image of diffraction pattern
im = Image.new('L',(N,N))
for i in range(N):
  for j in range(N):
    im.putpixel((i,j),63.0 + diffrac[i][j]*1920/diffracrng)
im.save('diffrac.gif')

# create FITS image of diffraction pattern
if 'diffrac.fits' in os.listdir('.'):
  os.unlink('diffrac.fits')
hdu = pyfits.PrimaryHDU(diffrac)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('diffrac.fits')

# print some flux stats
influx = numpy.sum(illum)
outflux = influx*pow(N*psf_scale,2)
print "Cropped PSF has %.2f%% of the flux; cropped airy has %.2f%%" % \
       (numpy.sum(psf)/outflux*100,numpy.sum(diffrac)/outflux*100)
flux = wavefront.flux_in(diffrac,N/2,N/2,1.22*psf_scale)
print "Inside first ring: ",flux/outflux
#print numpy.sum(illum)
