#! /usr/bin/env python
#---this is a python scrit that encapsulate the raytracing code from steve to make it more interactive/scriptable. 
# J Hao 3/28/2012 @ FNAL

try:
    import numpy as np
    import pyfits as pf
    import pylab as pl
    import os
    from DECamCCD_def import *
    from sklearn.decomposition import PCA
    from scipy.optimize import leastsq
    import mahotas as mh
    import scipy.ndimage as nd
except ImportError:
    print 'the required packages are: numpy, pyfits,pylab,scikit,scipy,mahotas'
    raise

#---------------calcuate moments ---------------
def moments(data):
    """
    Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments
    """
    total = data.sum()
    if total != 0.:
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        col = data[:, int(y)]
        width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
        row = data[int(x), :]
        width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
        height = data.max()
    else:
        height=0
        x=0
        y=0
        width_x=0
        width_y=0
    return height, x, y, width_x, width_y

def gaussianMoments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = leastsq(errorfunction, params)
    return p

#-----define the parameters --------
# this parameter will be written into the header

raypattern = 18
npix = 40
scale = 0.27
fwhm = 0.5
zenith = 0.
filter = 'g'
theta = 0.
corrector = 'corrector'
x = None
y = None
z = 0.1
output='temp.fit'
#------------------------------

def decamspot(xmm=None,ymm=None,seeing=0,npix=40,zenith=0,filter='g', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None):
    #---generating the .par file------
    file = open('temp.par','w')
    file.write('RAYPATTERN '+str(raypattern) +'\n')
    file.write('NPIX '+str(npix) +'\n')
    file.write('SCALE '+str(scale)+'\n')
    file.write('FWHM '+str(fwhm)+'\n')
    file.write('ZENITH '+str(zenith)+'\n')
    file.write('FILTER '+filter +'\n')
    file.write('XMM '+str(xmm)+'\n')
    file.write('YMM '+str(ymm)+'\n')
    if suband is None:
        file.write('WEIGHTS 1 0.9 0.8 0.7 0.6 \n')
    elif suband == 1:
        file.write('WEIGHTS 1 0.00001 0.00001 0.00001 0.00001\n') #approximately monochromatic
    elif suband == 2:
        file.write('WEIGHTS 0.00001 1 0.00001 0.00001 0.00001 \n')
    elif suband == 3:
        file.write('WEIGHTS 0.00001 0.00001 1 0.00001 0.00001 \n')
    elif suband == 4:
        file.write('WEIGHTS 0.00001 0.00001 0.00001 1 0.00001 \n')
    elif suband == 5:
        file.write('WEIGHTS 0.00001 0.00001 0.00001 0.00001 1 \n')
    file.write('THETA '+corrector+' '+str(theta)+'\n')
    if x is not None:
        file.write('X '+corrector+' '+str(x)+'\n')
    if y is not None:
        file.write('Y '+corrector+' '+str(y)+'\n')
    if z is not None:
        file.write('Z '+corrector+' '+str(z)+'\n')
    file.write('OUTPUT '+output+'\n')
    file.close()
    #---execute the raytrace code ------
    os.system('raytrace-3.13/decamspot '+'temp.par')
    #---output the result as an image vector
    b=pf.getdata('temp.fit')
    seeing = seeing/scale # convert arcseconds to pixel
    b=nd.filters.gaussian_filter(b,(seeing,seeing))
    hdr = pf.getheader('temp.fit')
    bb = b.reshape(npix*npix)
    pos = np.array([xmm,ymm])
    return np.concatenate((pos,bb)),hdr


def genImgV(filename=None,Nstar=None,ccd=None,seeing=0,npix=40,zenith=0,filter='g', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None):
    """
    seeing is the rms in arcseconds
    """
    datalist = []
    hdrlist = []
    randfactor=np.array([-1,1])
    if ccd is None:
        for i in range(Nstar):
            xmm = np.random.rand()*225*randfactor[np.random.randint(0,2)]
            ymm = np.random.rand()*225*randfactor[np.random.randint(0,2)]
            res = decamspot(xmm=xmm,ymm=ymm,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z,suband=suband)
            datalist.append(res[0])
            hdrlist.append(res[1])
        data = np.array(datalist)
    else:
        for i in range(Nstar):
            if i == 0:
                xmm = ccd[1]
                ymm = ccd[2]
            else:
                xmm = np.random.rand()*15*randfactor[np.random.randint(0,2)] + ccd[1]
                ymm = np.random.rand()*30*randfactor[np.random.randint(0,2)] + ccd[2]
            res = decamspot(xmm=xmm,ymm=ymm,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z,suband=suband)
            datalist.append(res[0])
            hdrlist.append(res[1])
        data = np.array(datalist)
    if filename is not None:
        hdu = pf.PrimaryHDU(data)
        hdu.header.add_comment('RAYPATTERN: '+str(raypattern))
        hdu.header.add_comment('NPIX: '+str(npix))
        hdu.header.add_comment('SCALE: '+str(scale))
        hdu.header.add_comment('FWHM: '+str(fwhm))
        hdu.header.add_comment('ZENITH: '+str(zenith))
        hdu.header.add_comment('FILTER: '+filter)
        hdu.header.add_comment('THETA: '+str(theta))
        hdu.header.add_comment('CORRECTOR: '+str(corrector))
        hdu.header.add_comment('X: '+str(x))
        hdu.header.add_comment('Y: '+str(y))
        hdu.header.add_comment('Z: '+str(z))
        hdu.header.add_comment('Seeing: '+str(seeing))
        if os.path.exists(filename):
            os.system('rm '+filename)
            hdu.writeto(filename)
        else:
            hdu.writeto(filename)        
    return data,hdrlist


def genImgVfixedPos(filename=None,seeing=0,npix=40,zenith=0,filter='g', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None):
    """
    seeing is the rms in arcseconds
    """
    datalist = []
    hdrlist = []
    xmm,ymm = np.meshgrid([-120,-80,-40,-20,0,20,40,80,120],[-120,-80,-40,-20,0,20,40,80,120])
    xmm = xmm.flatten()
    ymm = ymm.flatten()
    Nstar = len(xmm)
    for i in range(Nstar):
        res = decamspot(xmm=xmm[i],ymm=ymm[i],seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z,suband=suband)
        datalist.append(res[0])
        hdrlist.append(res[1])
    data = np.array(datalist)
    if filename is not None:
        hdu = pf.PrimaryHDU(data)
        hdu.header.add_comment('RAYPATTERN: '+str(raypattern))
        hdu.header.add_comment('NPIX: '+str(npix))
        hdu.header.add_comment('SCALE: '+str(scale))
        hdu.header.add_comment('FWHM: '+str(fwhm))
        hdu.header.add_comment('ZENITH: '+str(zenith))
        hdu.header.add_comment('FILTER: '+filter)
        hdu.header.add_comment('THETA: '+str(theta))
        hdu.header.add_comment('CORRECTOR: '+str(corrector))
        hdu.header.add_comment('X: '+str(x))
        hdu.header.add_comment('Y: '+str(y))
        hdu.header.add_comment('Z: '+str(z))
        hdu.header.add_comment('Seeing: '+str(seeing))
        if os.path.exists(filename):
            os.system('rm '+filename)
            hdu.writeto(filename)
        else:
            hdu.writeto(filename)        
    return data,hdrlist


   
def disImg(data=None,colorbar=False):
    """
    data is a vector with first, second as the center position, then is an image vector.
    """
    #pl.figure()
    size = np.sqrt(len(data[2:]))
    xmm = data[0]
    ymm = data[1]
    pl.matshow(data[2:].reshape(size,size),fignum=0)
    if colorbar == True:
        pl.colorbar()
    pl.xlim(0,size-1)
    pl.ylim(0,size-1)
    pl.xlabel('Pixels')
    pl.ylabel('Pixels')
    #pl.title('x ='+str(round(xmm,4)) +', y = '+str(round(ymm,4)) + ' (mm)')
    pl.grid(color='yellow')


def disImgAll(imgV=None):
    nrow,ncol = imgV.shape
    for i in range(nrow):
        x=imgV[i,0]*1000./15.
        y=imgV[i,1]*1000./15.
        size = np.sqrt(len(imgV[i,2:]))
        img = imgV[i,2:].reshape(size,size)
        xm = np.round(np.arange(size) - x,3)
        ym = np.round(np.arange(size) - y,3)
        pl.matshow(img)

def disImgCCD(imgV=None,ccd=None):
    """
    the CCD names are S1,...S31, N1...N31 
    """
    x = imgV[:,0]
    y = imgV[:,1]
    ok = (x >= ccd[1]-15)*(x<=ccd[1]+15)*(y<=ccd[2]+30)*(y>=ccd[2]-30)
    imgVV = imgV[ok,:]
    xx = imgVV[:,0]
    yy = imgVV[:,1]
    dd = xx**2+yy**2
    idx = np.argsort(dd)
    imgVV = imgVV[idx,:]
    nstar = imgVV.shape[0]
    ncol = int(round(np.sqrt(nstar)))
    if nstar/float(ncol)>nstar/ncol:
        nrow = nstar/ncol + 1
    else:
        nrow = nstar/ncol
    pl.figure(figsize=(nrow*5,ncol*5))
    for i in range(nstar):
        pl.subplot(min(nrow,ncol),max(nrow,ncol),i+1)
        disImg(imgVV[i,:],colorbar=False)
    return 'The image is done!'


def imgCCDctr(ccd=None,filter='g',seeing=0,z=None,theta=None):
    xmm = ccd[1]
    ymm = ccd[2]
    res = genImgV(Nstar=1, ccd = ccd,seeing=seeing,theta=theta,z=z)
    data = res[0]
    xcen = res[1][0]['xcen']
    ycen = res[1][0]['ycen']
    disImgCCD(data,ccd)
    pl.figtext(0.2,0.85,'CCD: '+ccd[0], color='r')
    pl.figtext(0.2,0.8,'Filter: '+filter, color='r')
    pl.figtext(0.2,0.75, 'xtarget: '+str(xmm) + ',  ytarget: '+str(ymm), color='r')
    pl.figtext(0.2,0.7, 'xcen: '+str(xcen) + ',  ycen: '+str(ycen), color='r')
    height, x, y, width_x, width_y = moments(data[0][2:].reshape(npix,npix))
    return height, x, y, width_x, width_y

def decompPCA(data=None,comp=None):
    img = data[:,2:].T
    size = np.sqrt(img.shape[0])
    pca = PCA(n_components=20)
    imgNew = pca.fit_transform(img)
    pl.matshow(imgNew[:,comp-1].reshape(size,size))
    pl.title('The '+str(comp)+' PC')
    pl.colorbar()
    return pca.explained_variance_ratio_

def fitPCAimg(coef=None, data = None, maxcomp = None):
    """
    get the coefficients of an image in terms of the PCA eigen images
    """
    img = data[:,2:].T
    size = np.sqrt(img.shape[0])
    pca = PCA(n_components=20)
    imgBasis = pca.fit_transform(img)
    nimg = img.shape[1]
    #def residule(basis=None,maxcomp=None):

def centroidChange(ccd=None, filter=None, suband=None):
    res = genImgV(Nstar=1,ccd=ccd,filter=filter,suband=suband)
    data = res[0]
    xcen = res[1][0]['xcen']
    ycen = res[1][0]['ycen']
    size = np.sqrt(len(data[0][2:]))
    img = data[0][2:].reshape(size,size)
    xmm = data[0][0]
    ymm = data[0][1]
    xcentroid,ycentroid = nd.center_of_mass(img)
    return xcentroid, ycentroid, xcen, ycen

def addseeing(imgV=None,seeing=None):
    img = imgV[0]

    
def decompZernike(ccd=None,seeing=0.,theta=0.,x=None,y=None,z=None,filter='g'):
    imgV=genImgV(Nstar=1, ccd = ccd,seeing=seeing,theta = theta,filter=filter,x=x,y=y,z=z)[0]
    size = np.sqrt(len(imgV[0][2:]))
    img = imgV[0][2:].reshape(size,size)        
    coeff = mh.zernike.zernike_moments(img,30,degree = 3, cm=mh.center_of_mass(img))
    coeff = coeff
    names = ['piston','tilt','defocus','astigmatism','coma','trefoil']
    res = dict(zip(names,coeff))
    disImgCCD(imgV,ccd)
    return res

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


def psfSizeCCD(ccd=None,filter='g',seeing=0,theta=0., zenith = 0., x=None, y=None,z=None):
    res = genImgV(Nstar=1, ccd = ccd,seeing=seeing,theta=theta,zenith=zenith,x=x,y=y,z=z)
    imgV = res[0]
    size = np.sqrt(len(imgV[0][2:]))
    img = imgV[0][2:].reshape(size,size)        
    xcen = res[1][0]['xcen']
    ycen = res[1][0]['ycen']
    height, x, y, width_x, width_y = fitgaussian(img)
    disImgCCD(imgV,ccd)
    pl.figtext(0.2,0.85,'CCD: '+ccd[0], color='r')
    pl.figtext(0.2,0.8,'Filter: '+filter, color='r')
    return width_x, width_y


def bfplane(x, y, z):
    """
    fit a best fit plane z = a*x +b*y +c
    return: a, b, c
    """
    n = float(len(x))
    A = np.array([[sum(x*x),sum(x*y),sum(x)],[sum(x*y),sum(y*y),sum(y)],[sum(x),sum(y),n]])
    B = np.array([sum(x*z),sum(y*z),sum(z)])
    res = np.linalg.solve(A,B)
    return res

def psfSizeAll(Nstar=None,filter='g',npix=40,seeing=0,theta=0., zenith = 0.,corrector='corrector', x=None, y=None,z=None):
    res = genImgV(Nstar=Nstar,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z)
    imgV = res[0]
    x = imgV[:,0]
    y = imgV[:,1]
    x_width = np.zeros(Nstar)
    y_width = np.zeros(Nstar)
    for i in range(Nstar):
        img = imgV[i,2:].reshape(npix,npix)
        hight,xtmp,tmp, x_width[i],y_width[i] = moments(img)
    return x, y, x_width, y_width, np.sqrt(x_width**2+y_width**2)

def psfSizeAllZernike(Nstar=None,filter='g',npix=40,seeing=0,theta=0., zenith = 0.,corrector='corrector', x=None, y=None,z=None,rand=False):
    if rand is False:
        res = genImgVfixedPos(seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z)
    else:
        res=genImgV(Nstar=Nstar,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z)
    imgV = res[0]
    x = imgV[:,0]
    y = imgV[:,1]
    coeff=[]
    Nstar = len(x)
    for i in range(Nstar):
        img = imgV[i,2:].reshape(npix,npix)
        coeff.append(mh.zernike.zernike_moments(img,30,degree = 3, cm=mh.center_of_mass(img)))
    return x, y, np.array(coeff)

def centroidChangeband(side=None):
    """
    This code calculate the centroid change of stars at different positions of the FP as the band filter changes
    """
    Nccd = len(side)
    xmmg = np.zeros(Nccd)
    ymmg = np.zeros(Nccd)
    xceng = np.zeros(Nccd)
    yceng = np.zeros(Nccd)
    xmmr = np.zeros(Nccd)
    ymmr = np.zeros(Nccd)
    xcenr = np.zeros(Nccd)
    ycenr = np.zeros(Nccd)
    xmmi = np.zeros(Nccd)
    ymmi = np.zeros(Nccd)
    xceni = np.zeros(Nccd)
    yceni = np.zeros(Nccd)
    xmmz = np.zeros(Nccd)
    ymmz = np.zeros(Nccd)
    xcenz = np.zeros(Nccd)
    ycenz = np.zeros(Nccd)
    ccdname = []
    for i in range(Nccd):
        print i
        xmmg[i],ymmg[i],xceng[i],yceng[i] = centroidChange(ccd=side[i], filter='g')
        xmmr[i],ymmr[i],xcenr[i],ycenr[i] = centroidChange(ccd=side[i], filter='r')
        xmmi[i],ymmi[i],xceni[i],yceni[i] = centroidChange(ccd=side[i], filter='i')
        xmmz[i],ymmz[i],xcenz[i],ycenz[i] = centroidChange(ccd=side[i], filter='z')
        ccdname.append(side[i][0])
    xrg = xcenr - xceng
    xig = xceni - xceng
    xzg = xcenz - xceng
    yrg = ycenr - yceng
    yig = yceni - yceng
    yzg = ycenz - yceng
    pl.subplot(2,1,1)
    pl.plot(xrg*1000./15.,'b.',label='r-band vs. g-band')
    pl.plot(xig*1000./15.,'r.',label='i-band vs. g-band')
    pl.plot(xzg*1000./15.,'g.',label='z-band vs. g-band')
    pl.xticks(np.arange(Nccd),ccdname)
    pl.xlabel('CCD position')
    pl.ylabel('x centroid difference (Pixels)')
    pl.legend(loc='best')
    pl.hlines(0,-1,31,linestyle='dashed',colors='k')
    pl.xlim(-1,31)
    pl.ylim(-1.5,1.5)
    pl.subplot(2,1,2)
    pl.plot(yrg*1000./15.,'b.',label='r-band vs. g-band')
    pl.plot(yig*1000./15.,'r.',label='i-band vs. g-band')
    pl.plot(yzg*1000./15.,'g.',label='z-band vs. g-band')
    pl.xticks(np.arange(Nccd),ccdname)
    pl.xlabel('CCD position')
    pl.ylabel('y centroid difference (Pixels)')
    pl.legend(loc='best')
    pl.hlines(0,-1,31,linestyle='dashed',colors='k')
    pl.xlim(-1,31)
    pl.ylim(-1.5,1.5)
    return '--- done!---'

def centroidChangeFP(filter='g',suband=None):
    xmm,ymm = np.meshgrid([-230,-180,-120,-40,0,40,120,180,230],[-230,-180,-120,-40,0,40,120,180,230])
    xmm = xmm.flatten()
    ymm = ymm.flatten()
    Nstar = len(xmm)
    xcen = np.zeros(Nstar)
    ycen = np.zeros(Nstar)
    for i in range(Nstar):
        res=decamspot(xmm=xmm[i],ymm=ymm[i],seeing=0,npix=40,zenith=0,filter=filter, theta=0., corrector='corrector',x=None,y=None,z=None,suband=suband)
        xcen[i] = res[1]['xcen']
        ycen[i] = res[1]['ycen']
    return xcen,ycen


def centroidSuband(filter=None):
    """
    This code calcualte the centroid change in each subband. 
    """
    xceng1,yceng1 = centroidChangeFP(filter=filter,suband = 1)
    xceng2,yceng2 = centroidChangeFP(filter=filter,suband = 2)
    xceng3,yceng3 = centroidChangeFP(filter=filter,suband = 3)
    xceng4,yceng4 = centroidChangeFP(filter=filter,suband = 4)
    xceng5,yceng5 = centroidChangeFP(filter=filter,suband = 5)
    r = np.sqrt(xceng1**2 + yceng1**2)
    pl.figure(figsize=(10,7))
    pl.subplot(2,1,1)
    pl.plot(r,(xceng2-xceng1)*1000./15.,'b.',label='sub2 - sub1')
    pl.plot(r,(xceng3-xceng1)*1000./15.,'r.',label='sub3 - sub1')
    pl.plot(r,(xceng4-xceng1)*1000./15.,'g.',label='sub4 - sub1')
    pl.plot(r,(xceng5-xceng1)*1000./15.,'c.',label='sub5 - sub1')
    pl.legend(loc='lower left')
    pl.hlines(0,0,300,color='k',linestyle='dashed')
    pl.xlim(0,300)
    pl.xlabel('distance to the FP center (mm)')
    pl.ylabel('x centroid difference (pixel)')
    pl.title('Centroid Change for subands of filter: '+filter)
    pl.subplot(2,1,2)
    pl.plot(r,(yceng2-yceng1)*1000./15.,'b.',label='sub2 - sub1')
    pl.plot(r,(yceng3-yceng1)*1000./15.,'r.',label='sub3 - sub1')
    pl.plot(r,(yceng4-yceng1)*1000./15.,'g.',label='sub4 - sub1')
    pl.plot(r,(yceng5-yceng1)*1000./15.,'c.',label='sub5 - sub1')
    pl.legend(loc='lower left')
    pl.hlines(0,0,300,color='k',linestyle='dashed')
    pl.xlim(0,300)
    pl.xlabel('distance to the FP center (mm)')
    pl.ylabel('y centroid difference (pixel)')
    return xceng1, yceng1, xceng2, yceng2,xceng3,yceng3,xceng4,yceng4,xceng5,yceng5
    



if __name__ == '__main__':
    import healpy as hp
    # ----the centroid change for different filters -----
    """
    S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23,S24,S25,S26,S27,S28,S29,S30,S31]

    N=[N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31]
    
    pl.figure(figsize=(15,10))
    centroidChangeband(S)
    pl.savefig('/home/jghao/research/decamFocus/figures/Scentroid_band.png')
    pl.figure(figsize=(15,10))
    centroidChangeband(N)
    pl.savefig('/home/jghao/research/decamFocus/figures/Ncentroid_band.png')


    #-----as function of radial -----
    # ----different suband ----
    colnames = ['xcen_sub1','ycen_sub1','xcen_sub2','ycen_sub2','xcen_sub3','ycen_sub3','xcen_sub4','ycen_sub4','xcen_sub5','ycen_sub5']
    
    res = centroidSuband(filter='g')
    pl.savefig('/home/jghao/research/decamFocus/figures/centroid_change_suband_g.png')
    hp.mwrfits('/home/jghao/research/decamFocus/data_g.fit',res,colnames = colnames)
    res = 0
    pl.close()

    res=centroidSuband(filter='r')
    pl.savefig('/home/jghao/research/decamFocus/figures/centroid_change_suband_r.png')
    hp.mwrfits('/home/jghao/research/decamFocus/data_r.fit',res,colnames = colnames)
    res = 0
    pl.close()

   
    res=centroidSuband(filter='i')
    pl.savefig('/home/jghao/research/decamFocus/figures/centroid_change_suband_i.png')
    hp.mwrfits('/home/jghao/research/decamFocus/data_i.fit',res,colnames = colnames)
    res = 0
    pl.close()
    
  
    res=centroidSuband(filter='z')
    pl.savefig('/home/jghao/research/decamFocus/figures/centroid_change_suband_z.png')
    hp.mwrfits('/home/jghao/research/decamFocus/data_z.fit',res,colnames = colnames)
    res = 0
    pl.close()
    
    
    #----for different band----
    xceng,yceng = centroidChangeFP(filter='g')
    xcenr,ycenr = centroidChangeFP(filter='r')
    xceni,yceni = centroidChangeFP(filter='i')
    xcenz,ycenz = centroidChangeFP(filter='z')
   
    r = np.sqrt(xceng**2 + yceng**2)
    pl.figure(figsize=(10,7))
    pl.subplot(2,1,1)
    pl.plot(r,(xcenr - xceng)*1000./15.,'b.',label = 'r - g')
    pl.plot(r,(xceni - xceng)*1000./15.,'g.',label = 'i - g')
    pl.plot(r,(xcenz - xceng)*1000./15.,'r.',label = 'z - g')
    pl.xlabel('distance to the FP center (mm)')
    pl.ylabel('x centroid difference (pixel)')
    pl.legend(loc='best')
    pl.hlines(0,0,300,color='k',linestyle='dashed')
    pl.xlim(0,300)
    pl.subplot(2,1,2)
    pl.plot(r,(ycenr - yceng)*1000./15.,'b.',label = 'r - g')
    pl.plot(r,(yceni - yceng)*1000./15.,'g.',label = 'i - g')
    pl.plot(r,(ycenz - yceng)*1000./15.,'r.',label = 'z - g')
    pl.xlabel('distance to the FP center (mm)')
    pl.ylabel('y centroid difference (pixel)')
    pl.legend(loc='best')
    pl.hlines(0,0,300,color='k',linestyle='dashed')
    pl.xlim(0,300)
    pl.savefig('/home/jghao/research/decamFocus/figures/centroid_change_fp.png')
    """

#----generate the data ------
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar1000_notilt_nodefocus.fit',Nstar=1000,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_notilt_z0.01mm.fit',Nstar=1000,z=0.01,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_notilt_z0.1mm.fit',Nstar=1000,z=0.1,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_tilt_10arcsec_nodefocus.fit',Nstar=1000,theta=10.,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_tilt_100arcsec_nodefocus.fit',Nstar=1000,theta=100.,seeing=0.9)
