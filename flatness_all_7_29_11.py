#!/usr/bin/env python
# examine the flatness
import sys
sys.path.append('/home/jghao/fmccd')
from DECamCCD import *


#---------------------upper left-----------

#step 1: set the directory

baseDir = '/data/jiangang/flat_allccd_7_29_11/upperleft/'

#step2: separate extensions to files

NameFits=gl.glob(baseDir+'*.fits')
NameFits.sort()
ccd=['s30','s31','s28','s27','s22','s23','s24','s19','s18','s17','s11','s12','s13','s7','s6','s5','s4','n4','n5','n6','n7']
ext=['s30','s31','s28','s27','s22','s23','s24','s19','s18','s17','s11','s12','s13','s7','s6','s5','s4','n4','n5','n6','n7']

fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))

for i in range(len(ext)):
    t=extract_extension(NameFits[i],ext[i])

NameBias=baseDir+'bias/bias.fits'
for i in range(len(ext)):
    t=extract_extension(NameBias,ext[i])

#step 2: subtract the bias by extension

NameExt=gl.glob(baseDir+'*exten*.fits')
NameExt.sort()
for i in range(len(ext)):
    hdu=pf.open(NameExt[i])
    bias=pf.getdata(baseDir+'bias/bias_exten_'+ext[i]+'.fits')
    hdu[0].data = hdu[0].data - bias
    hdu[0].header.update('CTYPE2', 'DEC--TAN')
    hdu[0].header.update('CTYPE1', 'RA--TAN')    
    output=NameExt[i][0:-5]+'_sub.fits'
    hdu.writeto(output)

#step 3: extract catalog

name = gl.glob(baseDir+'*sub.fits')
name.sort()
Nfile=len(name)
sexdir='/home/jghao/software/sextractor-2.5.0/config/'

for i in range(Nfile):
        img_name=name[i]
        output=img_name[0:-5]+'_catalog.fits'
        ckimg=img_name[0:-5]+'check.fits'
	t=sex(img_name,output,thresh=2.0,fwhm=4,gain=None,zmag=None,sexdir=sexdir,scale=0.125,check_img=ckimg,sigma_map=None,config="default.sex")
	
 

#---------------------upper right-----------

#step 1: set the directory

baseDir = '/data/jiangang/flat_allccd_7_29_11/upperright/'

#step2: separate extensions to files

NameFits=gl.glob(baseDir+'*.fits')
NameFits.sort()


ccd=['s30','s29','s25','s26','s22','s21','s20','s14','s15','s16','s10','s9','s8','s1','s2','s3','s4','n4','n3','n2','n1']

ext=['s30','s29','s25','s26','s22','s21','s20','s14','s15','s16','s10','s9','s8','s1','s2','s3','s4','n4','n3','n2','n1']

fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))

for i in range(len(ext)):
    t=extract_extension(NameFits[i],ext[i])

NameBias=baseDir+'bias/bias.fits'
for i in range(len(ext)):
    t=extract_extension(NameBias,ext[i])

#step 2: subtract the bias by extension

NameExt=gl.glob(baseDir+'*exten*.fits')
NameExt.sort()
for i in range(len(ext)):
    hdu=pf.open(NameExt[i])
    bias=pf.getdata(baseDir+'bias/bias_exten_'+ext[i]+'.fits')
    hdu[0].data = hdu[0].data - bias
    hdu[0].header.update('CTYPE2', 'DEC--TAN')
    hdu[0].header.update('CTYPE1', 'RA--TAN')    
    output=NameExt[i][0:-5]+'_sub.fits'
    hdu.writeto(output)

#step 3: extract catalog

name = gl.glob(baseDir+'*sub.fits')
name.sort()
Nfile=len(name)
sexdir='/home/jghao/software/sextractor-2.5.0/config/'

for i in range(Nfile):
        img_name=name[i]
        output=img_name[0:-5]+'_catalog.fits'
        ckimg=img_name[0:-5]+'check.fits'
	t=sex(img_name,output,thresh=2.0,fwhm=4,gain=None,zmag=None,sexdir=sexdir,scale=0.125,check_img=ckimg,sigma_map=None,config="default.sex")
	
 
#---------------------lower left-----------

#step 1: set the directory

baseDir = '/data/jiangang/flat_allccd_7_29_11/lowerleft/'

#step2: separate extensions to files

NameFits=gl.glob(baseDir+'*.fits')
NameFits.sort()

ccd=['s7','s6','s5','s4','n4','n5','n6','n7','n13','n12','n11','n17','n18','n19','n24','n23','n22','n27','n28','n31','n30']

ext=['s7','s6','s5','s4','n4','n5','n6','n7','n13','n12','n11','n17','n18','n19','n24','n23','n22','n27','n28','n31','n30']


fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))

for i in range(len(ext)):
    t=extract_extension(NameFits[i],ext[i])

NameBias=baseDir+'bias/bias.fits'
for i in range(len(ext)):
    t=extract_extension(NameBias,ext[i])

#step 2: subtract the bias by extension

NameExt=gl.glob(baseDir+'*exten*.fits')
NameExt.sort()
for i in range(len(ext)):
    hdu=pf.open(NameExt[i])
    bias=pf.getdata(baseDir+'bias/bias_exten_'+ext[i]+'.fits')
    hdu[0].data = hdu[0].data - bias
    hdu[0].header.update('CTYPE2', 'DEC--TAN')
    hdu[0].header.update('CTYPE1', 'RA--TAN')    
    output=NameExt[i][0:-5]+'_sub.fits'
    hdu.writeto(output)

#step 3: extract catalog

name = gl.glob(baseDir+'*sub.fits')
name.sort()
Nfile=len(name)
sexdir='/home/jghao/software/sextractor-2.5.0/config/'

for i in range(Nfile):
        img_name=name[i]
        output=img_name[0:-5]+'_catalog.fits'
        ckimg=img_name[0:-5]+'check.fits'
	t=sex(img_name,output,thresh=2.0,fwhm=4,gain=None,zmag=None,sexdir=sexdir,scale=0.125,check_img=ckimg,sigma_map=None,config="default.sex")
	
 

#---------------------lower right-----------

#step 1: set the directory

baseDir = '/data/jiangang/flat_allccd_7_29_11/lowerright/'

#step2: separate extensions to files

NameFits=gl.glob(baseDir+'*.fits')
NameFits.sort()

ccd=['s1','s2','s3','s4','n4','n3','n2','n1','n8','n9','n10','n16','n15','n14','n20','n21','n22','n26','n25','n29','n30']

ext=['s1','s2','s3','s4','n4','n3','n2','n1','n8','n9','n10','n16','n15','n14','n20','n21','n22','n26','n25','n29','n30']

fileNo=range(len(ccd))
ccd2fileNo=dict(zip(ccd,fileNo))

for i in range(len(ext)):
    t=extract_extension(NameFits[i],ext[i])

NameBias=baseDir+'bias/bias.fits'
for i in range(len(ext)):
    t=extract_extension(NameBias,ext[i])

#step 2: subtract the bias by extension

NameExt=gl.glob(baseDir+'*exten*.fits')
NameExt.sort()
for i in range(len(ext)):
    hdu=pf.open(NameExt[i])
    bias=pf.getdata(baseDir+'bias/bias_exten_'+ext[i]+'.fits')
    hdu[0].data = hdu[0].data - bias
    hdu[0].header.update('CTYPE2', 'DEC--TAN')
    hdu[0].header.update('CTYPE1', 'RA--TAN')    
    output=NameExt[i][0:-5]+'_sub.fits'
    hdu.writeto(output)

#step 3: extract catalog

name = gl.glob(baseDir+'*sub.fits')
name.sort()
Nfile=len(name)
sexdir='/home/jghao/software/sextractor-2.5.0/config/'

for i in range(Nfile):
        img_name=name[i]
        output=img_name[0:-5]+'_catalog.fits'
        ckimg=img_name[0:-5]+'check.fits'
	t=sex(img_name,output,thresh=2.0,fwhm=4,gain=None,zmag=None,sexdir=sexdir,scale=0.125,check_img=ckimg,sigma_map=None,config="default.sex")
	
 


























