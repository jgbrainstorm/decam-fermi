#! /usr/bin/env python
#---this is a python scrit that encapsulate the raytracing code from steve to make it more interactive/scriptable.
# the zernike polynomial definition codes is adopted and modified from:
#http://www.staff.science.uu.nl/~werkh108/docs/teach/2011b_python/python102/examples/py102-example2-zernike.py
# J Hao 3/28/2012 @ FNAL
# This is copied from decamspotPY on 6/18/2012 to delete those unnecessary routines. 

try:
    import numpy as np
    import pyfits as pf
    import pylab as pl
    import os,sys
    from DECamCCD_def import *
    import scipy.ndimage as nd
    import healpy as hp
    import glob as gl
    from scipy.misc import factorial as fac
    from scipy.signal import convolve2d
    from scipy.signal import fftconvolve
    import scipy.signal as sg
    import binplot as bp
    from scipy.optimize import leastsq
    import cPickle as p
    import sklearn.neighbors as nb
    from sklearn.svm import SVR
except ImportError:
    print 'the required packages are: numpy, pyfits,pylab,scikit,scipy,mahotas'
    raise

#-----define the parameters ----------------------------------
# this parameter will be written into the header
install_dir ='/home/jghao/research/ggsvn/decam-fermi/pyRaytrace/'
raypattern = 18
npix = 160
scale = 0.27/4.  # arcseconds/pix
#scale = 0.27
diffusionfwhm = 0.4
zenith = 0.
filter = 'r'
theta = 0.
corrector = 'corrector'
x = None
y = None
z = None
output='temp.fit'
#----------------------------------------------------------------


#---------------calcuate moments ---------------

def rebin(a, new_shape):
    """
    example:
    c = rebin(b, (2, 3))
    """
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)

def sech(x,width,height):
    """
    hyperbolic secant function
    """
    z = x/width
    res = height*2./(np.exp(z)+np.exp(-z))
    return res

def s2profile(r,r0,A,B):
    """
    hyperbolic secant square function
    """
    x = r/r0
    res = A*4./(np.exp(x)+np.exp(-x))**2 + B
    return res

def gprofile(r,sig,A,B):
    """
    Fit the binned distribution to a gaussian profile with a constant
    """
    res = A*np.exp(-0.5*(r/sig)**2)+B
    return res

def mprofile(r, alpha, beta,A,B):
    """
    Fit the light distribution to a Moffat profile
    """
    res = A*(1+(r/alpha)**2)**(-beta)+B
    return res

def gfwhm(img):
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    A0,sig0 = moments(img)
    radius = np.sqrt(row**2+col**2)
    img = img.flatten()
    ok = img >0
    img = img[ok]
    radius = radius.flatten()
    radius = radius[ok]
    def residualg(p,r,I):
        sig,A,B = p
        Ierr = np.sqrt(abs(I))
        res = (gprofile(radius,sig,A,B) - I)/Ierr
        return res 
    B0 = 0.
    p0=np.array([sig0,A0,B0])
    p = leastsq(residualg,p0,args=(radius,img))[0]
    sig,A,B = p
    fwhm_gauss= 2. * sig * np.sqrt(2. * np.log(2.))
    return sig,A,B,fwhm_gauss

def s2fwhm(img):
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    A0,r0_0 = moments(img)
    radius = np.sqrt(row**2+col**2)
    img = img.flatten()
    ok = img >0
    img = img[ok]
    radius = radius.flatten()
    radius = radius[ok]
    def residuals2(p,r,I):
        r0,A,B = p
        Ierr = np.sqrt(abs(I))
        res = (s2profile(radius,r0,A,B) - I)/Ierr
        return res 
    B0 = 0.
    p0=np.array([r0_0,A0,B0])
    p = leastsq(residuals2,p0,args=(radius,img))[0]
    r0,A,B = p
    fwhm_sech2= 1.7627471*r0 # obtained by solving the equation
    return r0,A,B,fwhm_sech2

def gaussian2d(x,y,xc,yc,sigmax,sigmay,rho,A,B):
    res = A*np.exp(-0.5/(1-rho**2)*(x**2/sigmax**2+y**2/sigmay**2-2.*rho*x*y/(sigmax*sigmay)))+B
    return res

def g2dfwhm(img):
    """
    x is col, y is row
    sigmax - sigmac, sigmay - sigmar
    """
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    A0,sigmac0 = moments(img)
    sigmar0 = sigmac0
    rho0 = 0.
    B0 = 0.
    p0=np.array([sigmac0,sigmar0,rho0,A0, B0])
    def residualg2d(p,x,y,xc,yc,I):
        sigmax,sigmay,rho,A,B = p
        Ierr = np.sqrt(abs(I))+0.00001 # to avoid those = 0, add a small number 
        res = (gaussian2d(x,y,xc,yc,sigmax,sigmay,rho,A,B) - I)/Ierr
        return res.flatten()
    p = leastsq(residualg2d,p0,args=(col,row,colCen,rowCen,img))[0]
    sigmac,sigmar,rho,A,B = p
    Mcc = sigmac**2
    Mrr = sigmar**2
    Mrc = rho**2*Mcc*Mrr
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    whiskerLength = np.sqrt(np.abs(M22))
    lambdap = 0.5*(M20 + abs(M22))
    lambdam = 0.5*(M20 - abs(M22))
    fwhm_g2d = np.sqrt(2.*np.log(2.))*(np.sqrt(lambdap)+np.sqrt(lambdam))
    #fwhm = np.sqrt(M20/2.)*2.35482*scale
    return A, B, whiskerLength, fwhm_g2d

def mfwhm(img):
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,1.1/scale)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    radius = np.sqrt(row**2+col**2)
    A0,alpha0 = moments(img)
    beta0=1.5
    B0 = 0.
    p0=np.array([alpha0,beta0,A0, B0])
    img = img.flatten()
    ok = img >0
    img = img[ok]
    radius = radius.flatten()
    radius = radius[ok]
    def residualm(p,r,I):
        alpha,beta,A,B = p
        Ierr = np.sqrt(abs(I))
        res = (mprofile(radius,alpha,beta,A,B) - I)/Ierr
        return res 
    p = leastsq(residualm,p0,args=(radius,img))[0]
    #rad,im,imerr=bp.bin_scatter(radius,img,binsize=1,fmt='b.',plot=False)
    alpha,beta,A,B = p
    fwhm_moffat= 2. * abs(alpha) * np.sqrt(2.**(1./beta)-1)
    return alpha,beta,A,B,fwhm_moffat

def wfwhm(img,sigma):
    """
    calculate the fwhm using the weighted moments after subtracting the weights.
    This is using the corrected moments calcualtion, i.e., adaptiveMomentsNew
    """
    nrow,ncol=img.shape
    Isum = img.sum()
    Icol = img.sum(axis=0) # sum over all rows
    Irow = img.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = img*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        IWsum = IWmat.sum()
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/IWsum
        dcolmean = np.sum((colgrid-colmean)*IWcol)/IWsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mrr = np.sum(rowgrid**2*IWrow)/IWsum
    Mcc = np.sum(colgrid**2*IWcol)/IWsum
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/IWsum
    Cm = np.matrix([[Mcc,Mrc],[Mrc,Mrr]])
    Cw = np.matrix([[sigma**2,0.],[0.,sigma**2]])
    Cimg = (Cm.I - Cw.I).I
    Mcc = Cimg[0,0]
    Mrr = Cimg[1,1]
    Mrc = Cimg[0,1]
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    e1 = M22.real/M20.real
    e2 = M22.imag/M20.real
    whiskerLength = np.sqrt(np.abs(M22))
    lambdap = 0.5*(M20 + abs(M22))
    lambdam = 0.5*(M20 - abs(M22))
    fwhmw = np.sqrt(2.*np.log(2.))*(np.sqrt(lambdap)+np.sqrt(lambdam))
    return e1,e2,whiskerLength,fwhmw 


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
    return height,np.sqrt(width_x**2 + width_y**2)

def wr(x,y,xcen,ycen,sigma):
    """
    Returns a gaussian weight function with the given parameters
    """
    res=np.exp(-((x-xcen)**2+(y-ycen)**2)/(2.*sigma**2))/(2.*np.pi*sigma**2) 
    return res

def gauss_seeing(npix = None,fwhm=None,e1=None,e2=None):
    """
    generate a seeing PSF of given fwhm and e1 and e2
    fwhm in the unit of arcsec
    """
    fwhm = fwhm/scale
    M20 = (fwhm/2.35482)**2
    row,col = np.mgrid[-npix/2:npix/2,-npix/2:npix/2]
    rowc = row.mean()
    colc = col.mean()
    Mcc = 0.5*M20*(1+e1)
    Mrc = 0.5*e2*M20
    Mrr = 0.5*M20*(1-e1)
    rho = Mrc/np.sqrt(Mcc*Mrr)
    img = np.exp(-0.5/(1-rho**2)*(row**2/Mrr + col**2/Mcc - 2*rho*row*col/np.sqrt(Mrr*Mcc)))
    res = img/img.sum()
    return res


def moffat_seeing(npix = None, alpha=None,beta=None):
    row,col = np.mgrid[-npix/2:npix/2,-npix/2:npix/2]
    rowc = row.mean()
    colc = col.mean()
    img = (beta - 1)/(np.pi*alpha**2)/(1+((row**2+col**2)/alpha**2))**beta
    res = img/img.sum()
    return res

                 
   

def convolveH(image=None,kernel=None):
    """
    using the fft convolve full mode and then trim the size to the original size. This is because if use the same mode, the centroid will change due to edge effects.
    """
    res = sg.fftconvolve(image,kernel,mode='full')
    res = res/res.sum()
    rowM,colM = nd.center_of_mass(res)
    size = image.shape[0]/2
    resnew = res[int(rowM)-size:int(rowM)+size,int(colM)-size:int(colM)+size]
    return resnew


def adaptiveCentroid(data=None,sigma=None):
    """
    calculate the centroid using the adaptive approach
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/np.sum(IWrow)
        dcolmean = np.sum((colgrid-colmean)*IWcol)/np.sum(IWcol)
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break

    return rowmean,colmean

         
def des_image(exptime=100,mag=None, Nstar=1,ccd=None,seeing=[0.9,0.,0.],npix=npix,zenith=0,filter='g', theta=0., phi=0,corrector='corrector',x=None,y=None,z=None,suband=None,regular=False,setbkg=True):
    """
    This code generate a PSF star with seeing and sky background
    exptime is given in sec
    """
    gain = 0.21 # convert electrons to ADU
    zeropoint = 26.794176 # r band, from Nikolay
    objectphoton = exptime*10**(0.4*(zeropoint - mag))
    if setbkg == False:
        skyphoton = 0.
    else:
        skyphoton = 8.460140*exptime
    bkg = skyphoton*gain
    res = genImgV(filename=None,Nstar=Nstar,ccd=ccd,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, phi=phi,corrector=corrector,x=x,y=y,z=z,suband=suband,regular=regular)
    if Nstar == 1:
        psf = res[0][0][4:].reshape(npix,npix)
        img = (psf * objectphoton + skyphoton)*gain
        img = img + np.sqrt(img)
        img = rebin(img,(40,40))
        psf = rebin(psf,(40,40))
    if Nstar > 1:
        img = []
        psf = []
        for i in range(Nstar):
           psfi = res[0][i][4:].reshape(npix,npix)
           imgi = (psf * objectphoton + skyphoton)*gain
           imgi = imgi+np.sqrt(imgi)
           imgi = rebin(imgi,(40,40))
           psfi = rebin(psfi,(40,40))
           img.append(imgi = np.sqrt(imgi))
           psf.append(psfi)
    return img,bkg,psf


def des_psf_image(exptime=100,mag=None,seeing=[0.9,0.,0.],npix=npix,setbkg=True):
    
    """
    This code generate a PSF star with seeing and sky background (no optics psf)
    exptime is given in sec
    """
    gain = 0.21 # convert electrons to ADU
    zeropoint = 26.794176 # r band, from Nikolay
    objectphoton = exptime*10**(0.4*(zeropoint - mag))
    if setbkg == False:
        skyphoton = 0.
    else:
        skyphoton = 8.460140*exptime
    bkg = skyphoton*gain
    psf = gauss_seeing(npix,seeing[0],seeing[1],seeing[2])
    img = (psf * objectphoton + skyphoton)*gain
    img = img + np.sqrt(img)
    img = rebin(img,(40,40))
    psf = rebin(psf,(40,40))
    return img,bkg,psf


def complexMoments(data=None,sigma=None):
    """
    This one calcualte the 3 2nd moments and 4 thrid moments with the Gaussian weights.
    col : x direction
    row : y direction
    the centroid is using the adpative centroid.
    sigma is the stand deviation of the measurement kernel in pixel
    The output is in pixel coordinate
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        IWsum = IWmat.sum()
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/IWsum
        dcolmean = np.sum((colgrid-colmean)*IWcol)/IWsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mr = np.sum(rowgrid*IWrow)/IWsum
    Mc = np.sum(colgrid*IWcol)/IWsum
    Mrr = np.sum(rowgrid**2*IWrow)/IWsum
    Mcc = np.sum(colgrid**2*IWcol)/IWsum
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/IWsum
    Mrrr = np.sum(rowgrid**3*IWrow)/IWsum
    Mccc = np.sum(colgrid**3*IWcol)/IWsum
    Mrrc = np.sum(np.outer(rowgrid**2,colgrid)*IWmat)/IWsum
    Mrcc = np.sum(np.outer(rowgrid,colgrid**2)*IWmat)/IWsum
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    M31 = complex(3*Mc - (Mccc+Mrrc)/sigma**2, 3*Mr - (Mrcc + Mrrr)/sigma**2)
    M33 = complex(Mccc-3*Mrrc, 3.*Mrcc - Mrrr)
    return M20, M22, M31, M33



#------------------------------

def decamspot(xmm=None,ymm=None,seeing=[0.9,0.,0.],npix=None,zenith=0,filter='r', theta=0., phi=0,corrector='corrector',x=None,y=None,z=None,suband=None):
    #---generating the .par file------
    dir = os.getcwd()+'/'
    file = open(dir+'temp.par','w')
    file.write('RAYPATTERN '+str(raypattern) +'\n')
    file.write('NPIX '+str(npix) +'\n')
    file.write('SCALE '+str(scale)+'\n')
    file.write('FWHM '+str(diffusionfwhm)+'\n')
    file.write('ZENITH '+str(zenith)+'\n')
    file.write('FILTER '+filter +'\n')
    file.write('PHI '+corrector+' '+str(phi) +'\n')
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
    file.write('PHI '+corrector+' '+str(phi)+'\n')
    if x is not None:
        file.write('X '+corrector+' '+str(x)+'\n')
    if y is not None:
        file.write('Y '+corrector+' '+str(y)+'\n')
    if z is not None:
        file.write('Z '+corrector+' '+str(z)+'\n')
    file.write('OUTPUT '+dir+output+'\n')
    file.close()
    #---execute the raytrace code ------
    os.system(install_dir+'raytrace-3.13/decamspot '+dir+'temp.par')
    #---output the result as an image vector
    b=pf.getdata(dir+'temp.fit')
    if seeing != 0.:
        b=addseeingImgFFT(b,fwhm=seeing[0],e1=seeing[1],e2=seeing[2])
    hdr = pf.getheader(dir+'temp.fit')
    ypstamp,xpstamp = nd.center_of_mass(b) # y -> row, x-> col
    bb = b.reshape(npix*npix)
    pos = np.array([hdr['xcen'],hdr['ycen'],xpstamp,ypstamp])
    os.system('rm '+dir+'temp.fit temp.par')
    return np.concatenate((pos,bb)),hdr


def genImgV(filename=None,Nstar=None,ccd=None,seeing=[0.9,0.,0.],npix=None,zenith=0,filter='g', theta=0., phi=0,corrector='corrector',x=None,y=None,z=None,suband=None,regular=False):
    """
    seeing is the rms in arcseconds
    syntax: genImgV(filename=None,Nstar=None,ccd=None,seeing=0,npix=40,zenith=0,filter='g', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None)
    """
    if regular == True:
        Nstar = 16
        regridx,regridy = np.meshgrid(np.array([-13., -6.5, 6.5, 13]), np.array([-28., -14., 14., 28.]))
        regridx = regridx.flatten()
        regridy = regridy.flatten()
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
            if i == 0 and regular == False:
                xmm = ccd[1]
                ymm = ccd[2]
            else:
                if regular == False:
                    xmm = np.random.rand()*13*randfactor[np.random.randint(0,2)] + ccd[1]
                    ymm = np.random.rand()*28*randfactor[np.random.randint(0,2)] + ccd[2]
                else:
                    xmm = regridx[i]+ccd[1]
                    ymm = regridy[i]+ccd[2]
            res = decamspot(xmm=xmm,ymm=ymm,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, phi=phi,corrector=corrector,x=x,y=y,z=z,suband=suband)
            datalist.append(res[0])
            hdrlist.append(res[1])
        data = np.array(datalist)
    if filename is not None:
        hdu = pf.PrimaryHDU(data)
        hdu.header.set('RAYPATT',raypattern)
        hdu.header.set('NPIX',npix)
        hdu.header.set('SCALE',scale)
        hdu.header.set('FWHM',diffusionfwhm)
        hdu.header.set('ZENITH',zenith)
        hdu.header.set('FILTER',filter)
        hdu.header.set('THETA',theta)
        hdu.header.set('CORRT',corrector)
        hdu.header.set('PHI',phi)
        if x != None:
            hdu.header.set('X',x)
        if y != None:
            hdu.header.set('Y',y)
        if z != None:
            hdu.header.set('Z',z)
        if seeing != 0.:
            hdu.header.set('s_fwhm',seeing[0])
            hdu.header.set('e1',seeing[1])
            hdu.header.set('e2',seeing[2])
        if os.path.exists(filename):
            os.system('rm '+filename)
            hdu.writeto(filename)
        else:
            hdu.writeto(filename)        
    return data,hdrlist


def genImgVfixedPos(filename=None,seeing=0,npix=None,zenith=0,filter='r', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None):
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
        hdu.header.set('RAYPATT',raypattern)
        hdu.header.set('NPIX',npix)
        hdu.header.set('SCALE',scale)
        hdu.header.set('FWHM',fwhm)
        hdu.header.set('ZENITH',zenith)
        hdu.header.set('FILTER',filter)
        hdu.header.set('THETA',theta)
        hdu.header.set('CORRT',corrector)
        if x != None:
            hdu.header.set('X',x)
        if y != None:
            hdu.header.set('Y',y)
        if z != None:
            hdu.header.set('Z',z)
        if seeing != 0.:
            hdu.header.set('s_fwhm',seeing[0])
            hdu.header.set('e1',seeing[1])
            hdu.header.set('e2',seeing[2])
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
    size = np.sqrt(len(data[4:]))
    xmm = data[0]
    ymm = data[1]
    pl.matshow(data[4:].reshape(size,size),fignum=False)
    if colorbar == True:
        pl.colorbar()
    pl.xlim(0,size-1)
    pl.ylim(0,size-1)
    pl.xlabel('Pixels')
    pl.ylabel('Pixels')
    pl.grid(color='yellow')


def disImgAll(imgV=None):
    nrow,ncol = imgV.shape
    for i in range(nrow):
        x=imgV[i,0]*1000./15.
        y=imgV[i,1]*1000./15.
        size = np.sqrt(len(imgV[i,4:]))
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

  
def imgCCDctr(ccd=None,filter='r',seeing=[0.9,0.,0.],x=0,y=0,z=0.,theta=0.,contour=False,sigma=1.08/scale,npix=None,phi=0,zenith=0):
    xmm = ccd[1]
    ymm = ccd[2]
    res = genImgV(Nstar=1, ccd = ccd,seeing=seeing,theta=theta,x=x,y=y,z=z,npix=npix,zenith=zenith,phi=phi)
    data = res[0]
    img = data[0][4:].reshape(npix,npix)
    npix = npix/4
    scale = 0.27
    sigma= sigma/4.
    img = rebin(img,(npix,npix))
    mfit = mfwhm(img)
    gfit = gfwhm(img)
    s2fit = s2fwhm(img)
    g2dfit = g2dfwhm(img)
    wfit = wfwhm(img,sigma=sigma)
    pl.figure(figsize=(18,6))
    pl.subplot(1,3,1)
    pl.matshow(img,fignum=False)
    stampImg = img.copy()
    pl.xlabel('Pixel')
    pl.ylabel('Pixel')
    pl.grid(color='y')
    rowCen,colCen = adaptiveCentroid(data=img,sigma=sigma)
    M20, M22, M31, M33 =complexMoments(stampImg,sigma=sigma)
    pl.figtext(0.15,0.8, 'e1: '+str(round(wfit[0],3)) + ',  e2: '+str(round(wfit[1],3)), color='r')
    pl.figtext(0.15,0.75, 'rowCen: '+str(round(rowCen,4)) + ',  colCen: '+str(round(colCen,4)), color='r')
    pl.figtext(0.15,0.7, 'PSF whisker_Wmoments: '+str(round(wfit[2]*scale,4))+' [arcsec]', color='r')
    pl.figtext(0.15,0.65, 'PSF whisker_Amoments: '+str(round(g2dfit[2]*scale,4))+' [arcsec]', color='r')
    pl.figtext(0.15,0.6, 'CCD Position: '+ccd[0] +',   filter: '+filter, color='r')
    pl.subplot(1,3,2)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    radius = np.sqrt(row**2+col**2)
    img = img.flatten()
    radius = radius.flatten()
    idx = np.argsort(radius)
    img = img[idx]
    radius = radius[idx]
    rad,im,imerr=bp.bin_scatter(radius,img,binsize=1,fmt='bo',plot=False)
    halfmax = np.median(img[0:10])/2.
    pl.plot(radius,img,'k.')
    pl.grid(color='y')
    pl.hlines(halfmax,0,radius.max(),linestyle='solid',colors='b')
    pl.hlines(mfit[2]/2.,0,radius.max(),linestyle='solid',colors='r')
    pl.hlines(gfit[1]/2.,0,radius.max(),linestyle='solid',colors='g')
    pl.hlines(s2fit[1]/2.,0,radius.max(),linestyle='solid',colors='m')
    pl.hlines(g2dfit[0]/2.,0,radius.max(),linestyle='solid',colors='c',label='Adaptive Moments')
    #pl.vlines(fwhm/scale/2.,0, halfmax*4,linestyle='solid',colors='b',label='Weighted Moments')
    pl.vlines(wfit[3]/2.,0, halfmax*4,linestyle='solid',colors='b',label='Weighted Moments')
    pl.vlines(mfit[4]/2.,0, halfmax*4,linestyle='solid',colors='r')
    pl.vlines(gfit[3]/2.,0, halfmax*4,linestyle='solid',colors='g')
    pl.vlines(s2fit[3]/2.,0, halfmax*4,linestyle='solid',colors='m')
    pl.vlines(g2dfit[3]/2.,0, halfmax*4,linestyle='solid',colors='c')
    pl.plot(radius,mprofile(radius,mfit[0],mfit[1],mfit[2],mfit[3]),'r-',label='Moffat Fit')
    pl.plot(radius,gprofile(radius,gfit[0],gfit[1],gfit[2]),'g-',label='Gaussian Fit')
    pl.plot(radius,s2profile(radius,s2fit[0],s2fit[1],s2fit[2]),'m-',label='Sech2 Fit')
    pl.legend(loc='best')
    pl.ylim(0,halfmax*4)
    pl.xlim(0,npix/4.) 
    pl.xlabel('Radius [pixels]')
    pl.ylabel('Mean counts [ADU]')
    pl.title('Radial profile')
    pl.figtext(0.65,0.7,'Gaussian Weight '+r'$\sigma$: '+str(round(sigma*scale,3))+ ' arcsec',color='r')
    pl.figtext(0.65,0.6,'FWHM_Gaussian: '+str(round(gfit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.55,'FWHM_Moffat: '+str(round(mfit[4]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.5,'FWHM_Sech2: '+str(round(s2fit[3]*scale,3))+ ' arcsec')
    pl.figtext(0.65,0.45,'FWHM_Wmoments: '+str(round(wfit[3]*scale,3))+ ' arcsec') 
    pl.figtext(0.65,0.4,'FWHM_Amoments: '+str(round(g2dfit[3]*scale,3))+ ' arcsec') 
    pl.figtext(0.65,0.35,'M20: '+str(round(M20,5))+ ' pix')
    pl.figtext(0.65,0.3,'M22.real: '+str(round(M22.real,5))+ ' pix')
    pl.figtext(0.8,0.3,'M22.imag: '+str(round(M22.imag,5))+ ' pix')
    pl.figtext(0.65,0.25,'M31.real: '+str(round(M31.real,5))+ ' pix')
    pl.figtext(0.8,0.25,'M31.imag: '+str(round(M31.imag,5))+ ' pix')
    pl.figtext(0.65,0.2,'M33.real: '+str(round(M33.real,5))+ ' pix')
    pl.figtext(0.8,0.2,'M33.imag: '+str(round(M33.imag,5))+ ' pix')

    #return M20*scale**2,M22*scale**2,M31*scale**2,M33*scale**2
    return stampImg

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


def centroidChange(ccd=None, filter=None, suband=None):
    res = genImgV(Nstar=1,ccd=ccd,filter=filter,suband=suband)
    data = res[0]
    xcen = res[1][0]['xcen']
    ycen = res[1][0]['ycen']
    size = np.sqrt(len(data[0][4:]))
    img = data[0][4:].reshape(size,size)
    xmm = data[0][0]
    ymm = data[0][1]
    xcentroid,ycentroid = nd.center_of_mass(img)
    return xcentroid, ycentroid, xcen, ycen

def addseeingImgFFT(img = None,fwhm=1.,e1=0.,e2=0.):
    """
    fwhm input in the unit of the arcsec
    """
    kern = gauss_seeing(npix,fwhm=fwhm,e1=e1,e2=e2)
    img = img.astype('f') # required for the fftconvolve
    covimg = convolveH(img,kern)
    covimg = covimg/covimg.sum()
    return covimg

def addseeingImgFFTmoffat(img = None,alpha=None,beta=None):
    """
    fwhm input in the unit of the arcsec
    """
    kern = moffat_seeing(npix,alpha=alpha,beta=beta)
    img = img.astype('f') # required for the fftconvolve
    covimg = convolveH(img,kern)
    covimg = covimg/covimg.sum()
    return covimg



def addseeing(filename=None,fwhm=1.,e1=0.,e2=0.,fft=True):
    """
    fwhm is in the unit of the arcsec
    rebinfactor is an integer to indicate how much the image will be scaled down
    """
    hdu = pf.open(filename)
    n = len(hdu)
    hdu[0].header.set('s_fwhm',fwhm)
    hdu[0].header.set('e1',e1)
    hdu[0].header.set('e2',e2)
    kern = gauss_seeing(npix,fwhm=fwhm,e1=e1,e2=e2)
    for i in range(1,n):
        img=hdu[i].data[0][4:].reshape(npix,npix)
        img = img.astype('f')
        if fft == False:
            covimg = sg.convolve2d(img,kern,mode='same')
            newfname = filename.replace('_noseeing_','_withseeing_')+'_fwhm_'+str(fwhm)+'_e1_'+str(e1)+'_e2_'+str(e2)+'.fit'
        else:
            covimg = convolveH(img,kern)
            newfname = filename.replace('_noseeing_','_withseeing_')+'_fwhm_'+str(fwhm)+'_e1_'+str(e1)+'_e2_'+str(e2)+'_fftconvolve.fit'
        covimg = covimg/covimg.sum()
        hdu[i].data[0][4:] = covimg.flatten()
    hdu.writeto(newfname)
    #os.system('gzip '+newfname)
    return 'done'

def imageRebin(dir=None,rebinsize=(40,40)):
    """
    rebin the size of the image
    """
    if dir == None:
        dir = os.getcwd()+'/'
    filenameAll = gl.glob(dir+'*.fit*')
    filenameAll.sort()
    for filename in filenameAll:
        hdu = pf.open(filename)
        n = len(hdu)
        hdu[0].header.set('npix',rebinsize[0])
        hdu[0].header.set('scale',hdu[0].header.set('scale')*hdu[0].header.set('npix')/rebinsize[0])
        for i in range(1,n):
            img=hdu[i].data[0][4:].reshape(npix,npix)
            img = img.astype('f')
        if fft == False:
            covimg = sg.convolve2d(img,kern,mode='same')
            newfname = filename.replace('_noseeing_','_withseeing_')+'_fwhm_'+str(fwhm)+'_e1_'+str(e1)+'_e2_'+str(e2)+'.fit'
        else:
            covimg = sg.fftconvolve(img,kern,mode='same')
            newfname = filename.replace('_noseeing_','_withseeing_')+'_fwhm_'+str(fwhm)+'_e1_'+str(e1)+'_e2_'+str(e2)+'_fftconvolve.fit'
        covimg = covimg/covimg.sum()
        hdu[i].data[0][4:] = covimg.flatten()
    hdu.writeto(newfname)
    #os.system('gzip '+newfname)
    return 'done'



def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


def psfSizeCCD(ccd=None,filter='r',seeing=0,theta=0., zenith = 0., x=0., y=0.,z=0.):
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

def psfSizeAll(Nstar=None,filter='r',npix=None,seeing=0,theta=0., zenith = 0.,corrector='corrector', x=None, y=None,z=None):
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

def psfSizeAllZernike(Nstar=None,filter='r',npix=None,seeing=0,theta=0., zenith = 0.,corrector='corrector', x=None, y=None,z=None,rand=False):
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
        img = imgV[i,4:].reshape(npix,npix)
        coeff.append(mh.zernike.zernike_moments(img,30,degree = 3, cm=mh.center_of_mass(img)))
    return x, y, np.array(coeff)



def genPSFimage(filename=None):
    """
    convert the PSF image vector file to the set of PSF images
    """
    hdu=pf.open(filename)
    nn = len(hdu)
    for i in range(1,nn):
        img = hdu[i].data[0][4:].reshape(npix,npix)
        img = img/img.sum()
        hdu[i].data = img
    #hdu.scale('int16', '', bzero=32768)
    newfilename = filename[:-7]+'_stamp.fits'
    hdu.writeto(newfilename)
    os.system('gzip '+newfilename)    
        

def measureDataComplexM_multiext(filename,sigma = 1.1,scale=0.27):
    """
    Note that here, the sigma is not fwhm. Sigma is given in arcsec
    """
    hdu=pf.open(filename)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    sigma = sigma/scale
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        M31=np.zeros(Nobj).astype(complex)
        M33=np.zeros(Nobj).astype(complex)
        for i in range(Nobj):
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=hdui.data[i][4:].reshape(npix,npix),sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(M20), np.median(M22), np.median(M31), np.median(M33)])
    data=np.array(data)
    hp.mwrfits(filename[:-7]+'_complexMoments_gausswt_'+str(sigma*scale)+'.fit',data.T,colnames=colnames)
    return '---done !-----'

def wisker(data=None, sigma = None):
    """
    This code calculate the wisker length defined as sqrt(e1^2+e2^2)
    input: 
         data: 2d stamp image
         sigma: std of the Gaussian weight Kernel in pixel
    """
    M20, M22, M31, M33 = complexMoments(data=data,sigma=sigma)
    wisker_length = np.sqrt(M22.real**2+M22.imag**2)
    return wisker_length




def genImgVallCCD(filename=None,Nstar=None,seeing=[0.9,0.,0.],npix=None,zenith=0,filter='g', theta=0.,phi=0, corrector='corrector',x=0.,y=0.,z=0.,suband=None,regular=False):
    """
    Nstar is the number of stars on each CCD
    """
    hduList = pf.HDUList()
    hdu = pf.PrimaryHDU(np.array([0]))
    hdu.header.set('RAYPATT',raypattern)
    hdu.header.set('NPIX',npix)
    hdu.header.set('SCALE',scale)
    hdu.header.set('FWHM',diffusionfwhm)
    hdu.header.set('ZENITH',zenith)
    hdu.header.set('FILTER',filter)
    hdu.header.set('THETA',theta)
    hdu.header.set('CORRT',corrector)
    hdu.header.set('PHI',phi)
    if x != None:
        hdu.header.set('X',x)
    if y != None:
        hdu.header.set('Y',y)
    if z != None:
        hdu.header.set('Z',z)
    if seeing != 0.:
        hdu.header.set('s_fwhm',seeing[0])
        hdu.header.set('e1',seeing[1])
        hdu.header.set('e2',seeing[2])
    hduList.append(hdu)
    for ccd in N[1:]+S[1:]:
        print ccd
        res = genImgV(filename=filename,Nstar=Nstar,ccd=ccd,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, phi=phi,corrector=corrector,x=x,y=y,z=z,suband=suband,regular=regular)
        hdu = pf.PrimaryHDU(res[0])
        hdu.header.set('ccdPos',ccd[0])
        hdu.header.set('ccdXcen',ccd[1])
        hdu.header.set('ccdYcen',ccd[2])
        hduList.append(hdu)
    if filename != None:
        if os.path.exists(filename):
            os.system('rm '+filename)
            hduList.writeto(filename)
        else:
            hduList.writeto(filename)
        os.system('gzip '+filename)
    return hduList


def zernike_rad(m, n, rho):
    """
    Calculate the radial component of Zernike polynomial (m, n) 
    given a grid of radial coordinates rho.
    """
    if (n < 0 or m < 0 or abs(m) > n):
        raise ValueError
    if ((n-m) % 2):
        return rho*0.0
    pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
    return sum(pre_fac(k) * rho**(n-2.0*k) for k in xrange((n-m)/2+1))

def zernike(m, n, rho, phi):
    """
    Calculate Zernike polynomial (m, n) given a grid of radial
    coordinates rho and azimuthal coordinates phi.
    """
    if (m > 0): return zernike_rad(m, n, rho) * np.cos(m * phi)
    if (m < 0): return zernike_rad(-m, n, rho) * np.sin(-m * phi)
    return zernike_rad(0, n, rho)

def zernikel(j, rho, phi):
    """
    Calculate Zernike polynomial with Noll coordinate j given a grid of radial coordinates rho and azimuthal coordinates phi.
    """
    n = 0
    while (j > n):
        n += 1
        j -= n
    m = -n+2*j
    return zernike(m, n, rho, phi)


def zernikeFit(x, y, z,max_rad=225.,cm=[0,0],max_order=20):
    """
    Fit a set of x, y, z data to a zernike polynomial with the least square fitting. Note that here x, y, z are all 1 dim array. Here the max_rad is by default equal to 225 mm, the size of the decam focal plane.
    It will return the beta and the adjusted R2
    """
    x = x - cm[0]
    y = y - cm[1]
    n = len(x)
    p = max_order
    rho = np.sqrt(x**2+y**2)/max_rad #normalize to unit circle.
    phi = np.arctan2(y,x)
    dataX = []
    ok = rho <= 1.
    for j in range(max_order):
        dataX.append(zernikel(j,rho[ok],phi[ok]))
    dataX=np.array(dataX).T
    beta,SSE,rank,sing = np.linalg.lstsq(dataX,z[ok])# SSE is the residual sum square
    sigma = np.sqrt(SSE/(n-p))
    betaErr = sigma/np.dot(dataX.T,dataX).diagonal()
    SST = np.var(z[ok])*(len(z[ok])-1)# SST is the sum((z_i - mean(z))^2)
    R2 = 1 - SSE/SST
    R2adj = 1-(1-R2)*(len(z[ok])-1)/(len(z[ok])-max_order)# adjusted R2 for quality of fit.             
    return beta,betaErr, R2adj



def dispZernike(beta=1.,j=0,gridsize = 1, max_rad = 1):
    x,y = np.meshgrid(np.arange(-gridsize,gridsize,0.001),np.arange(-gridsize,gridsize,0.001))
    rho = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    ok = rho < max_rad
    znk = beta*zernikel(j,rho,phi)*ok
    pl.imshow(znk)
    return znk
    
def showZernike(beta=None,betaErr=None,gridsize = 1, max_rad = 1,significance=False):
    """
    significance shows how significant the coefficients are constrained. 
    """
    x,y = np.meshgrid(np.arange(-gridsize,gridsize,0.01),np.arange(-gridsize,gridsize,0.01))
    rho = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    ok = rho < max_rad
    if significance != False:
        sigIdx = np.abs(beta)/betaErr >= significance
        beta = beta[sigIdx]
    nn = len(beta)
    znk=0
    for j in range(nn):
        znk = znk + beta[j]*zernikel(j,rho,phi)*ok
    pl.imshow(znk)
    return znk


def zernike_diagnosis(Nstar=None,seeing=[0.9,0.,0.],npix=npix,zenith=0,filter='r', theta=0., phi=0,corrector='corrector',x=None,y=None,z=None,zernike_max_order=20,regular=False):
    """
    This function produce the zernike plots for a set of given parameters of the tilt/shift/defocus
    """
    hdu = genImgVallCCD(Nstar=Nstar,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta,phi=phi, corrector=corrector,x=x,y=y,z=z,regular=regular)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        M31=np.zeros(Nobj).astype(complex)
        M33=np.zeros(Nobj).astype(complex)
        #sigma = 1.1/0.27
        sigma = 1.08/scale
        for i in range(Nobj):
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=hdui.data[i][4:].reshape(npix,npix),sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(M20), np.median(M22), np.median(M31), np.median(M33)])
    data=np.array(data)
    pl.figure(figsize=(15,15))
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    pl.subplot(3,3,1)
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    znk=showZernike(beta=beta)
    pl.colorbar()
    pl.title(colnames[2])
    for i in range(3,6):
        pl.subplot(3,3,2*i-4)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        znk=showZernike(beta=beta)
        pl.colorbar()
        pl.title(colnames[i]+'_real')
        print '--- R2_adj of the fit is: '+str(R2_adj) +'---'
        pl.subplot(3,3,2*i-3)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        znk=showZernike(beta=beta)
        pl.colorbar()
        pl.title(colnames[i]+'_imag')
        print '--- R2_adj of the fit is: '+str(R2_adj) +'---'
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    return betaAll,betaErrAll, R2adjAll


def moments_display(Nstar=1,seeing=[0.9,0.,0.],npix=npix,zenith=0,filter='r', theta=0., phi=0,corrector='corrector',x=None,y=None,z=None,regular=False):
    """
    This function produce the zernike plots for a set of given parameters of the tilt/shift/defocus
    """
    hdu = genImgVallCCD(Nstar=Nstar,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta,phi=phi, corrector=corrector,x=x,y=y,z=z,regular=regular)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        M31=np.zeros(Nobj).astype(complex)
        M33=np.zeros(Nobj).astype(complex)
        for i in range(Nobj):
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=rebin(hdui.data[i][4:].reshape(npix,npix),(40,40)),sigma=2.)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(M20), np.median(M22), np.median(M31), np.median(M33)])
    data=np.array(data)
    datam = data.copy()
    #---add 20 percent noise --
    #data = addMomentsNoise(data,5.)
    # remove the mean for M31 and M33
    #data = subMeanM3x(data)
    pl.figure(figsize=(12,12))
    pl.subplot(2,2,1)
    phi22 = 0.5*np.arctan2(data[:,3].imag,data[:,3].real)
    x = data[:,0].real
    y = data[:,1].real
    phi22[x<0] = phi22+np.deg2rad(180)
    u = np.sqrt(np.abs(data[:,3]))*np.cos(phi22)
    v = np.sqrt(np.abs(data[:,3]))*np.sin(phi22)
    qvr = pl.quiver(x,y,u,v,width = 0.004, color='r',pivot='middle',headwidth=0.,headlength=0.,headaxislength=0.,scale_units='width')
    qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('X [mm]')
    pl.ylabel('Y [mm]')
    pl.title('M22')
    pl.subplot(2,2,2)
    phi31 = np.arctan2(data[:,4].imag,data[:,4].real)
    u = np.sqrt(np.abs(data[:,4]))*np.cos(phi31)
    v = np.sqrt(np.abs(data[:,4]))*np.sin(phi31)
    qvr=pl.quiver(x,y,u,v,width=0.003,color='r',pivot='middle',headwidth=4)
    qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('X [mm]')
    pl.ylabel('Y [mm]')
    pl.title('M31')
    pl.subplot(2,2,3)
    phi33 = np.arctan2(data[:,5].imag,data[:,5].real)/3.
    u = np.sqrt(np.abs(data[:,5]))*np.cos(phi33)
    v = np.sqrt(np.abs(data[:,5]))*np.sin(phi33)
    pl.quiver(x,y,u,v,width=0.003,color='r',headwidth=4)
    u = np.sqrt(np.abs(data[:,5]))*np.cos(phi33+np.deg2rad(120))
    v = np.sqrt(np.abs(data[:,5]))*np.sin(phi33+np.deg2rad(120))
    pl.quiver(x,y,u,v,width=0.003,color='r',headwidth=4)
    u = np.sqrt(np.abs(data[:,5]))*np.cos(phi33+np.deg2rad(240))
    v = np.sqrt(np.abs(data[:,5]))*np.sin(phi33+np.deg2rad(240))
    qvr=pl.quiver(x,y,u,v,width=0.003,color='r',headwidth=4)
    qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('X [mm]')
    pl.ylabel('Y [mm]')
    pl.title('M33')
    pl.subplot(2,2,4)
    u = np.sqrt(data[:,2].real) - np.sqrt(data[:,2].real).min()
    v = np.zeros(len(data[:,2].real))
    qvr = pl.quiver(x,y,u,v,width = 0.008, color='r',pivot='middle',headwidth=0.,headlength=0.,headaxislength=0.,scale_units='width')
    qk = pl.quiverkey(qvr, -150,-240,np.max(u),str(round(max(u),3))+' pix',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.grid(color='g')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.title('median '+r'$\sqrt{M20}$: '+str(round(np.median(scale*4*np.sqrt(data[:,2].real)),3))+' [arcsec]')
    return datam

def addMomentsNoise(data,percent):
    noise = data[:,2:]*percent
    noise = noise*np.random.randn(noise.shape[0],noise.shape[1])
    data[:,2:] = data[:,2:]+noise
    return data



def coeff_display(Nstar=1,seeing=[0.9,0.,0.],npix=npix,zenith=0,filter='r', theta=0., phi=0,corrector='corrector',x=0.,y=0.,z=0.,zernike_max_order=20,regular=False):
    """
    This function produce the zernike plots for a set of given parameters of the tilt/shift/defocus
    """
    hdu = genImgVallCCD(Nstar=Nstar,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta,phi=phi, corrector=corrector,x=x,y=y,z=z,regular=regular)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        M31=np.zeros(Nobj).astype(complex)
        M33=np.zeros(Nobj).astype(complex)
        #sigma = 1.1/0.27
        sigma = 1.08/scale
        for i in range(Nobj):
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=hdui.data[i][4:].reshape(npix,npix),sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(M20), np.median(M22), np.median(M31), np.median(M33)])
    data=np.array(data)
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    for i in range(3,6):
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    ind = np.arange(len(betaAll[0]))
    momname = ('M20','M22.Real','M22.imag','M31.real','M31.imag','M33.real','M33.imag')
    fmtarr = ['bo-','ro-','go-','co-','mo-','yo-','ko-']
    pl.figure(figsize=(17,13))
    for i in range(7):
        pl.subplot(7,1,i+1)
        pl.errorbar(ind,betaAll[i],yerr = betaErrAll[i],fmt=fmtarr[i])
        if i == 0:
            pl.title('x: '+str(hdu[0].header['x'])+'   y: '+str(hdu[0].header['y'])+'   z: '+str(hdu[0].header['z'])+'   tilt: '+str(hdu[0].header['theta'])+'   fwhm: '+str(hdu[0].header['s_fwhm'])+'   e1: '+str(hdu[0].header['e1'])+'   e2: '+str(hdu[0].header['e2']))
        pl.grid()
        pl.xlim(-1,21)
        if i ==0:
            pl.ylim(-10,65)
        elif i ==1:
            pl.ylim(-5,6)
        elif i ==2:
            pl.ylim(-5,6)
        elif i == 3:
            pl.ylim(-0.1,0.1)
        elif i == 4:
            pl.ylim(-0.1,0.1)
        elif i ==5:
            pl.ylim(-100,100)
        elif i == 6:
            pl.ylim(-100,100)
        pl.xticks(ind,('','','','','','','','','','','','','','','','','','','',''))
        pl.ylabel(momname[i])
    pl.xticks(ind,('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'))
    pl.xlabel('Zernike Coefficients')
    return betaAll,betaErrAll


    #pl.figure(figsize=(13,5))
    #ind = np.arange(len(betaAll[0]))
    #pl.pcolor(betaAll,vmin=-50,vmax=60)
    #pl.grid(color='red')
    #pl.yticks(np.arange(7)+0.5,('M20','M22.Real','M22.imag','M31.real','M31.imag','M33.real','M33.imag'))
    #pl.xticks(ind+0.5,('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'))
    #pl.xlabel('Zernike Coefficients')
    #pl.colorbar(shrink=1)
    #return betaAll,betaErrAll

def coeff_display_file(filename,zernike_max_order=20):
    """
    This function produce the zernike plots for a set of given parameters of the tilt/shift/defocus
    """
    hdu = pf.open(filename)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        M31=np.zeros(Nobj).astype(complex)
        M33=np.zeros(Nobj).astype(complex)
        #sigma = 1.1/0.27
        sigma = 1.08/scale
        for i in range(Nobj):
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=hdui.data[i][4:].reshape(npix,npix),sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(M20), np.median(M22), np.median(M31), np.median(M33)])
    data=np.array(data)
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    for i in range(3,6):
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    ind = np.arange(len(betaAll[0]))
    momname = ('M20','M22.Real','M22.imag','M31.real','M31.imag','M33.real','M33.imag')
    fmtarr = ['bo-','ro-','go-','co-','mo-','yo-','ko-']
    pl.figure(figsize=(17,13))
    for i in range(7):
        pl.subplot(7,1,i+1)
        pl.errorbar(ind,betaAll[i],yerr = betaErrAll[i],fmt=fmtarr[i])
        if i == 0:
            pl.title('x: '+str(hdu[0].header['x'])+'   y: '+str(hdu[0].header['y'])+'   z: '+str(hdu[0].header['z'])+'   tilt: '+str(hdu[0].header['theta'])+'   fwhm: '+str(hdu[0].header['s_fwhm'])+'   e1: '+str(hdu[0].header['e1'])+'   e2: '+str(hdu[0].header['e2']))
        pl.grid()
        pl.xlim(-1,21)
        if i ==0:
            pl.ylim(-50,100)
        elif i ==1:
            pl.ylim(-10,10)
        elif i ==2:
            pl.ylim(-10,10)
        elif i == 3:
            pl.ylim(-0.01,0.01)
        elif i == 4:
            pl.ylim(-0.01,0.01)
        elif i ==5:
            pl.ylim(-20,20)
        elif i == 6:
            pl.ylim(-20,20)
        pl.xticks(ind,('','','','','','','','','','','','','','','','','','','',''))
        pl.ylabel(momname[i])
    pl.xticks(ind,('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'))
    pl.xlabel('Zernike Coefficients')
    return betaAll,betaErrAll





def zernike_file(filename=None,zernike_max_order=20,significance=False):
    """
    This function produce the zernike plots for a set of given parameters of the tilt/shift/defocus
    """
    hdu = pf.open(filename)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        M31=np.zeros(Nobj).astype(complex)
        M33=np.zeros(Nobj).astype(complex)
        #sigma = 1.1/0.27
        sigma = 1.08/scale
        for i in range(Nobj):
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=hdui.data[i][4:].reshape(npix,npix),sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(M20), np.median(M22), np.median(M31), np.median(M33)])
    data=np.array(data)
    pl.figure(figsize=(10,10))
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    pl.subplot(3,3,1)
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    znk=showZernike(beta=beta)
    pl.colorbar()
    pl.title(colnames[2])
    for i in range(3,6):
        pl.subplot(3,3,2*i-4)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        znk=showZernike(beta=beta,betaErr=betaErr,significance=significance)
        pl.colorbar()
        pl.title(colnames[i]+'_real')
        print '--- R2_adj of the fit is: '+str(R2_adj) +'---'
        pl.subplot(3,3,2*i-3)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        znk=showZernike(beta=beta,betaErr=betaErr,significance=significance)
        pl.colorbar()
        pl.title(colnames[i]+'_imag')
        print '--- R2_adj of the fit is: '+str(R2_adj) +'---'
    pl.figtext(0.4,0.32,'------- optics/seeing parameters --------')
    pl.figtext(0.5,0.29,str(hdu[0].header.items()[6]))
    pl.figtext(0.5,0.27,str(hdu[0].header.items()[7]))
    pl.figtext(0.5,0.25,str(hdu[0].header.items()[10]))
    pl.figtext(0.5,0.23,str(hdu[0].header.items()[11]))
    pl.figtext(0.5,0.21,str(hdu[0].header.items()[13]))
    pl.figtext(0.5,0.19,str(hdu[0].header.items()[14]))
    pl.figtext(0.5,0.17,str(hdu[0].header.items()[15]))
    pl.figtext(0.5,0.15,str(hdu[0].header.items()[16]))
    pl.figtext(0.5,0.13,str(hdu[0].header.items()[17]))
    pl.figtext(0.5,0.11,str(hdu[0].header.items()[18]))
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    return betaAll,betaErrAll, R2adjAll




def zernike_diff(Nstar=1,seeing=[0.9,0.,0.],npix=npix,zenith=0,filter='g', theta=0., corrector='corrector',x=None,y=None,z=None,zernike_max_order=20,regular=False):
    """
    This function produce zernike coefficients and compare the difference
    """
    hdu = genImgVallCCD(Nstar=Nstar,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z,regular=regular)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        M20=np.zeros(Nobj)
        M22=np.zeros(Nobj).astype(complex)
        M31=np.zeros(Nobj).astype(complex)
        M33=np.zeros(Nobj).astype(complex)
        #sigma = 1.1/0.27
        sigma = 1.08/0.27
        for i in range(Nobj):
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=hdui.data[i][4:].reshape(npix,npix),sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(M20), np.median(M22), np.median(M31), np.median(M33)])
    data=np.array(data)
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    for i in range(3,6):
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    return betaAll,betaErrAll, R2adjAll

def comp_zernike(beta1=None,betaErr1=None,beta2=None,betaErr2=None):
    nmoments = 7
    chi2 = np.zeros(7)
    for i in range(nmoments):
        chi2[i] = np.sum((beta1[i,:] - beta2[i,:])**2/(betaErr1**2+betaErr2**2))
    return chi2



def dispZernike20():
    pl.figure(figsize=(20,16))
    for i in range(0,20):
        pl.subplot(4,5,i+1)
        dispZernike(j=i)
        pl.title('j = '+str(i))
    return 0

def zernike_coeff(filename=None,zernike_max_order=20):
    """
    measure the zernike coefficients and errors
    """
    hdu = pf.open(filename)
    nn = len(hdu)
    data = []
    colnames = ['x','y','M20','M22','M31','M33']
    sigma = 1.08/0.27
    for hdui in hdu[1:]:
        img = hdui.data[0][4:].reshape(npix,npix)
        img = rebin(img,(40,40))
        M20,M22,M31,M33=complexMoments(data=img,sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,M20,M22,M31,M33])
    data=np.array(data)
    betaAll=[]
    #betaErrAll=[]
    R2adjAll=[]
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    #betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    for i in range(3,6):
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        #betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        #betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
    betaAll = np.array(betaAll)
    #betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    x=hdu[0].header['x']
    y=hdu[0].header['y']
    z=hdu[0].header['z']
    phi = hdu[0].header['phi']
    theta=hdu[0].header['theta']
    s_fwhm=hdu[0].header['s_fwhm']
    e1=hdu[0].header['e1']
    e2=hdu[0].header['e2']
    return x,y,z,theta,phi,s_fwhm,e1,e2,betaAll,R2adjAll


def measure_zernike_coeff(filelist=None,computer=None):
    """
    In the save data file, the cols are:
    x,y,z,theta,phi,s_fwhm,e1,e2,betaAll
    """
    data=[]
    f = filelist
    n = len(f)
    for i in range(n):
        print i
        t = zernike_coeff(f[i])
        data.append(np.append(t[0:8],t[8].flatten()))
    data = np.array(data)
    if computer != None:
        datafile = open('/data/des09.a/data/jiangang_psf/finegrid/coefficients/zernike_coeff_finerGrid_'+computer+'.cp','w')
    else:
        datafile = open('zernike_coeff_data_matrix.cp','w')
    p.dump(data,datafile,2)
    return '---done--'

    
def rowcol2XY(row,col,CCD):
    """
    This code convert the row/col [in pixels] of a given CCD to the x, y
    of the Focal plane [in mm]. Assuming a constant pixel scale 0.015mm/pix
    Input: row, col of the object, the CCD name in a way (S1, S2, etc)
    The current convention is:
    each ccd, the origin of row and col is the south east corner.
    the direction row increase is West
    the direction col increase is North.
    In my Focal Plane definition file,
        positive X is South
        positive Y is East
    So, the row increase as -Y direction.
        the col increase as -X direction.
    """
    pixscale = 0.015 #mm/pix
    X = CCD[1]+1024*pixscale-(col*pixscale+pixscale/2.)
    Y = CCD[2]+2048*pixscale-(row*pixscale+pixscale/2.)
    return X,Y
    
def multimachine_addseeing(computer=None):
    machine = np.array(['des04','des05','des06','des07','des08','des09','des10'])
    #fwhm = [0.6, 0.8, 1.0, 1.2, 1.4]
    #e1 = [-0.08,-0.04,0,0.04,0.08]
    #e2 = [-0.08,-0.04,0,0.04,0.08]
    #fwhm = [0.8, 1.0,1.2]
    #e1 = [-0.08,0,0.08]
    #e2 = [-0.08,0,0.08]
    #allfile = gl.glob('/data/des07.b/data/jiangang/PSF_noseeing/lowres/*.gz')
    allfile=gl.glob('/data/des09.a/data/jiangang_psf/finegrid/psf_noseeing/*.fit')
    allfile=np.array(allfile)
    allfile.sort()
    nfile=len(allfile)
    nmachine = len(machine)
    machineIdx = np.arange(nmachine)
    mid = machineIdx[computer == machine]
    idx = np.arange(mid[0]*nfile/nmachine,(mid[0]+1)*nfile/nmachine)
    allfileComputer = allfile[idx]
    for fname in allfileComputer:
        fw=np.random.randint(800.,1400)*0.001
        e11 = np.random.randint(-800,800)*0.0001
        e22 = np.random.randint(-800,800)*0.0001
        t = addseeing(filename=fname,fwhm = fw,e1=e11,e2=e22)

def singlemachine_addseeing(dir=None):
    if dir == None:
        dir = os.getcwd()+'/'
    fwhm = [0.6, 0.8, 1.0, 1.2, 1.4]
    e1 = [-0.08,-0.04,0,0.04,0.08]
    e2 = [-0.08,-0.04,0,0.04,0.08]
    allfile = gl.glob(dir+'*.gz')
    allfile=np.array(allfile)
    allfile.sort()
    nfile=len(allfile)
    filecount = 0.
    for fname in allfile:
        filecount = filecount +1
        print 'file:'+str(filecount)
        for fw in fwhm:
            for e11 in e1:
                for e22 in e2:
                    t = addseeing(filename=fname,fwhm = fw,e1=e11,e2=e22)
    return ' ----- done ! -----'


def multimachine_psfgen(computer=None):
    machine = np.array(['des04','des05','des06','des07','des08','des09','des10'])
    #tiltrange = [-100,-80,-50,-20,0,50,80,100]
    #xshiftrange = [-1.0,-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0]
    #yshiftrange = [-1.0,-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0]
    #defocusrange = [-1.0,-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0]
    #nmachine = len(machine)
    #machineIdx = np.arange(nmachine)
    #mid = machineIdx[computer == machine]
    #tlt = tiltrange[mid]
    #tlt = 100.
    n=4000
    for i in range(n):
        print i
        xsft = np.random.randint(-1000,1000)*0.0001
        ysft = np.random.randint(-1000,1000)*0.0001
        tlt = np.random.randint(-4000,4000)*0.01
        defo=np.random.randint(-1000,1000)*0.0001
        phi = np.random.randint(0,18000)*0.01
        filename='/data/des09.a/data/jiangang_psf/finegrid/psf_noseeing/PSF_noseeing_theta'+str(tlt)+'_x_'+str(xsft)+'_y_'+str(ysft)+'_z_'+str(defo)+'_phi_'+str(phi)+'.fit'
                    #filename = '/data/des07.b/data/jiangang/PSF_noseeing/PSF_noseeing_theta'+str(tlt)+'_x_'+str(xsft)+'_y_'+str(ysft)+'_z_'+str(defo)+'.fit'
        t = genImgVallCCD(filename=filename,Nstar=1,seeing=0.,npix=npix,zenith=0,filter='r', theta=tlt,phi=phi, corrector='corrector',x=xsft,y=ysft,z=defo,suband=None,regular=False)
    return '----done!-----'




def multimachine_measure_zernike(machine=None):
    #machine = np.array(['des04','des05','des06','des07','des08','des09','des10'])
    #allfile=gl.glob('/data/des09.a/data/jiangang_psf/highres_small_psf_with_seeing/*.fit')
    allfile=gl.glob('/data/des09.a/data/jiangang_psf/finegrid/psf_withseeing/*.fit')
    allfile.sort()
    nfile=len(allfile)
    nmachine = 7
    mm = nfile/nmachine
    if machine == 'des04':
        files = allfile[0:mm]
        allfile=0
    if machine == 'des05':
        files = allfile[mm:2*mm]
        allfile=0
    if machine == 'des06':
        files = allfile[2*mm:3*mm]
        allfile=0
    if machine == 'des07':
        files = allfile[3*mm:4*mm]
        allfile=0
    if machine == 'des08':
        files = allfile[4*mm:5*mm]
        allfile=0
    if machine == 'des09':
        files = allfile[5*mm:6*mm]
        allfile=0
    if machine == 'des10':
        files = allfile[6*mm:7*mm]
        allfile=0
    t=measure_zernike_coeff(filelist=files,computer=computer)
    return '---done!----'
    


def getPCA(data):
    """
    row of the data matrix is observations, col is the variable
    """
    #covM = np.cov(data.T) #note that np.cov define row as variables, col as observations
    #corM = np.corrcoef(data.T) # we will use correlation matrix instead of cov.
    covM = np.cov(data.T)
    eigvalue,eigvector = np.linalg.eig(covM) # each col of the eigvector matrix corresponds to one eigenvalue. So, each col is the coeff of one component
    pca = np.dot(data,eigvector) # each col is one pca, each row is one obs in that pca. 
    return eigvalue,eigvector,pca
    
def SVMRegression(trainingObs,trainingParam,Obs):
    """
    using the SVM regression to get the parameter
    trainingObs: the zernike coefficients data matrix. each row is a new observation
    trainingParam: the hexapod configration and seeing, each row is a new obs.
    Obs: a given set of measured coefficients
    """
    svr = SVR(C=1.0, cache_size=2000, coef0=0.0, degree=3, epsilon=0.1, gamma=0.,  kernel='rbf', probability=False, shrinking=True, tol=0.00001)
    nparam=trainingParam.shape[1]
    vparam = []
    for i in range(0,nparam):
        print i
        svr.fit(trainingObs,trainingParam[:,i])
        vparam.append(svr.predict(Obs))
    vparam = np.array(vparam)
    return vparam.T
   


def KNeighborRegression(trainingObs,trainingParam,Obs,n_neighbors):
    """
    using the k nearest neighbor regression to get the parameter
    trainingObs: the zernike coefficients data matrix. each row is a new observation
    trainingParam: the hexapod configration and seeing, each row is a new obs.
    Obs: a given set of measured coefficients
    """
    #knn = nb.KNeighborsRegressor(algorithm='ball_tree',n_neighbors=n_neighbors,weights = 'distance')
    knn = nb.KNeighborsRegressor(algorithm='ball_tree',n_neighbors=n_neighbors)
    knn.fit(trainingObs,trainingParam)
    return knn.predict(Obs)
    

def standardizeData(tdata,vdata):
    """
    This code standardize the training data and validation data by the training data.
    """
    tmean = tdata.mean(axis=0)
    tstd = tdata.std(axis=0)
    tdataNew = (tdata - tmean)/tstd
    vdataNew = (vdata - tmean)/tstd
    return tdataNew, vdataNew




def validateFit(Tfile=None,Vfile=None,PCA=False,alg='NNR'):
    b=p.load(open(Tfile))
    vb = p.load(open(Vfile))
    nobs = len(b)
    tdata=b[:,8:]
    tpara=b[:,0:8]
    vdata = vb[:,8:]
    vparaTrue=vb[:,0:8]
    tdata,vdata = standardizeData(tdata,vdata)
    if PCA == True:
        evlue, eigvector,tdata=getPCA(tdata)
        vdata=np.dot(vdata,eigvector)
    if alg == 'NNR':
        vparaReg=KNeighborRegression(tdata,tpara,vdata,15)
    if alg == 'SVM':
        vparaReg = SVMRegression(tdata,tpara,vdata)
    pl.figure(figsize=(17,10))
    pl.subplot(2,3,1)
    bp.bin_scatter(vparaTrue[:,0],vparaReg[:,0],binsize=0.005,fmt='b.')
    pl.plot([vparaTrue[:,0].min(),vparaTrue[:,0].max()],[vparaTrue[:,0].min(),vparaTrue[:,0].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('x shift [mm]')
    pl.subplot(2,3,2)
    bp.bin_scatter(vparaTrue[:,1],vparaReg[:,1],binsize=0.005,fmt='b.')
    pl.plot([vparaTrue[:,1].min(),vparaTrue[:,1].max()],[vparaTrue[:,1].min(),vparaTrue[:,1].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('y shift [mm]')
    pl.subplot(2,3,3)
    bp.bin_scatter(vparaTrue[:,2],vparaReg[:,2],binsize=0.005,fmt='b.')
    pl.plot([vparaTrue[:,2].min(),vparaTrue[:,2].max()],[vparaTrue[:,2].min(),vparaTrue[:,2].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('defocus [mm]')
    pl.subplot(2,3,4)
    bp.bin_scatter(vparaTrue[:,3],vparaReg[:,3],binsize=3,fmt='b.')
    pl.plot([vparaTrue[:,3].min(),vparaTrue[:,3].max()],[vparaTrue[:,3].min(),vparaTrue[:,3].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('tilt angle [arcsec]')
    pl.subplot(2,3,5)
    bp.bin_scatter(vparaTrue[:,4],vparaReg[:,4],binsize=5,fmt='b.')
    pl.plot([vparaTrue[:,4].min(),vparaTrue[:,4].max()],[vparaTrue[:,4].min(),vparaTrue[:,4].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('phi angle [deg]')
    return vparaTrue,vparaReg

def remM3xZernike(tdata):
    """
    This code remove the 0th zernike coefficient for the M31, M33 from the training and validation data object. The data has 140 columns, from 0 - 59 are the 2nd moments coefficients. 60, 80, 100, 120 are the 0th coefficients for the 3rd moments. We remove them from the data structure.    
    """
    idx = np.concatenate((np.arange(0,60),np.arange(61,80),np.arange(81,100),np.arange(101,120),np.arange(121,140)))
    datanew = tdata[:,idx]
    return datanew


def subMeanM3x(data=None):
    """
    this code subtract the mean of the 3rd moments from the data. This is to remove the tracking errors. The data is a matrix with 6 columns. The columns are 
    x, y, M20, M22, M31, M33
    """
    datamean = data.mean(axis = 0)
    data[:,4:6] = data[:,4:6] - datamean[4:6]
    return data


def validateFitNew(Tfile=None,Vfile=None,PCA=False,alg='NNR',scatter=False):
    b=p.load(open(Tfile))
    vb = p.load(open(Vfile))
    nobs = len(b)
    tdata=b[:,8:28].copy()
    ttpara=b[:,0:5].copy()
    tpara = b[:,0:5].copy()
    vdata = vb[:,8:28].copy()
    vparaTrue=vb[:,0:5].copy()
    vvparaTrue=vb[:,0:5].copy()
    tpara[:,3] = ttpara[:,3]*np.cos(np.deg2rad(ttpara[:,4]))
    tpara[:,4] = ttpara[:,3]*np.sin(np.deg2rad(ttpara[:,4]))
    vparaTrue[:,3] = vvparaTrue[:,3]*np.cos(np.deg2rad(vvparaTrue[:,4]))
    vparaTrue[:,4] = vvparaTrue[:,3]*np.sin(np.deg2rad(vvparaTrue[:,4]))
    # remove the 0th coeff for 3rd moments
    #tdata = remM3xZernike(tdata)
    #vdata = remM3xZernike(vdata)
    tdata,vdata = standardizeData(tdata,vdata)
    if PCA == True:
        evlue, eigvector,tdata=getPCA(tdata)
        vdata=np.dot(vdata,eigvector)
    if alg == 'NNR':
        vparaReg=KNeighborRegression(tdata,tpara,vdata,10)
    if alg == 'SVM':
        vparaReg = SVMRegression(tdata,tpara,vdata)
    pl.figure(figsize=(17,10))
    pl.subplot(2,3,1)
    bp.bin_scatter(vparaTrue[:,0],vparaReg[:,0],binsize=0.01,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,0].min(),vparaTrue[:,0].max()],[vparaTrue[:,0].min(),vparaTrue[:,0].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('x shift [mm]')
    pl.xlim(-0.1,0.1)
    pl.ylim(-0.1,0.1)
    pl.subplot(2,3,2)
    bp.bin_scatter(vparaTrue[:,1],vparaReg[:,1],binsize=0.01,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,1].min(),vparaTrue[:,1].max()],[vparaTrue[:,1].min(),vparaTrue[:,1].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('y shift [mm]')
    pl.xlim(-0.1,0.1)
    pl.ylim(-0.1,0.1)
    pl.subplot(2,3,3)
    bp.bin_scatter(vparaTrue[:,2],vparaReg[:,2],binsize=0.01,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,2].min(),vparaTrue[:,2].max()],[vparaTrue[:,2].min(),vparaTrue[:,2].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('defocus [mm]')
    pl.xlim(-0.1,0.1)
    pl.ylim(-0.1,0.1)

    pl.subplot(2,3,4)
    bp.bin_scatter(vparaTrue[:,3],vparaReg[:,3],binsize=5,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,3].min(),vparaTrue[:,3].max()],[vparaTrue[:,3].min(),vparaTrue[:,3].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('tilt angle in x direction [arcsec]')
    pl.xlim(-40,40)
    pl.ylim(-40,40)

    pl.subplot(2,3,5)
    bp.bin_scatter(vparaTrue[:,4],vparaReg[:,4],binsize=5,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,4].min(),vparaTrue[:,4].max()],[vparaTrue[:,4].min(),vparaTrue[:,4].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('tilt angle in y direction [arcsec]')
    pl.xlim(-40,40)
    pl.ylim(-40,40)
    return vparaTrue,vparaReg


def genValidation(n=None):
    """
    generate the zernike coefficients files for validation
    """
    tilt = np.random.randint(-100,101,n)*1.#[-100,-80,-50,-20,0,50,80,100]
    xshift = np.random.randint(-1000,1001,n)*0.001#[-1.0,-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0]
    yshift = np.random.randint(-1000,1001,n)*0.001#[-1.0,-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0]
    defocus = np.random.randint(-1000,1001,n)*0.001#[-1.0,-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 1.0]
    fwhm = np.random.randint(-1400,1400,n)*0.001
    e1 = np.random.randint(-800,900,n)*0.0001
    e2 = np.random.randint(-800,900,n)*0.0001
    for i in range(n):
        print i
        filename='/home/jghao/research/decamFocus/psf_withseeing/validation_highres_small/randomseeing/PSF_validation_theta'+str(tilt[15])+'_x_'+str(xshift[15])+'_y_'+str(yshift[15])+'_z_'+str(defocus[15])+'_fwhm_'+str(fwhm[i])+'_e1_'+str(e1[i])+'_e2_'+str(e2[i])+'.fit'
        t = genImgVallCCD(filename=filename,Nstar=1,seeing=[fwhm[i],e1[i],e2[i]],npix=npix,zenith=0,filter='r', theta=tilt[15], corrector='corrector',x=xshift[15],y=yshift[15],z=defocus[15],suband=None,regular=False)
    
    filelist = gl.glob('/home/jghao/research/decamFocus/psf_withseeing/validation_highres_small/randomseeing/PSF_*.fit*')
    t=measure_zernike_coeff(filelist)
    return '-----done!-----'
 
def generateKNNobj():
    """
    This code generate the KNN object that fited by the training data using M20
    """
    Tfile='/home/jghao/research/decamFocus/psf_withseeing/finerGrid_coeff_matrix/zernike_coeff_finerGrid_training.cp'
    b=p.load(open(Tfile))
    nobs = len(b)
    #tdata=b[:,8:28].copy()
    tdata=b[:,9:28].copy() # remove the zero order zernike, i.e. remove the mean of the M20
    #-standardize the data. use this information in future validation data too.
    tmean = tdata.mean(axis=0)
    tstd = tdata.std(axis=0)
    tdata = (tdata - tmean)/tstd
    ttpara=b[:,0:5].copy()
    tpara = b[:,0:5].copy()
    tpara[:,3] = ttpara[:,3]*np.cos(np.deg2rad(ttpara[:,4]))
    tpara[:,4] = ttpara[:,3]*np.sin(np.deg2rad(ttpara[:,4]))
    knn = nb.KNeighborsRegressor(algorithm='ball_tree',n_neighbors=15)
    knn.fit(tdata,tpara)
    p.dump(knn,open('finerGridKnnObj_remMean.cp','w'),2)
    p.dump([tmean,tstd],open('finerGridStdConst_remMean.cp','w'),2)
    #np.savetxt('finerGridStdConst.txt',np.array([tmean,tstd]),fmt='%f10.5',delimiter = ',')
    return 'It is done !'
    

def validateFitKnnObj(Vfile=None,scatter=True):
    """
    This function validate the results using the saved KNN object. It use only M20. 
    """
    vb = p.load(open(Vfile))
    vdata = vb[:,9:28].copy() # remove the 0 order of the M20, i.e. remove the mean of M20
    vparaTrue=vb[:,0:5].copy()
    vvparaTrue=vb[:,0:5].copy()
    vparaTrue[:,3] = vvparaTrue[:,3]*np.cos(np.deg2rad(vvparaTrue[:,4]))
    vparaTrue[:,4] = vvparaTrue[:,3]*np.sin(np.deg2rad(vvparaTrue[:,4]))

    knn = p.load(open('finerGridKnnObj_remMean.cp','r'))
    tmean,tstd = p.load(open('finerGridStdConst_remMean.cp','r'))
    vdata = (vdata - tmean)/tstd
    vparaReg = knn.predict(vdata)
    pl.figure(figsize=(17,10))
    pl.subplot(2,3,1)
    bp.bin_scatter(vparaTrue[:,0],vparaReg[:,0],binsize=0.01,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,0].min(),vparaTrue[:,0].max()],[vparaTrue[:,0].min(),vparaTrue[:,0].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('x shift [mm]')
    pl.xlim(-0.1,0.1)
    pl.ylim(-0.1,0.1)
    pl.subplot(2,3,2)
    bp.bin_scatter(vparaTrue[:,1],vparaReg[:,1],binsize=0.01,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,1].min(),vparaTrue[:,1].max()],[vparaTrue[:,1].min(),vparaTrue[:,1].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('y shift [mm]')
    pl.xlim(-0.1,0.1)
    pl.ylim(-0.1,0.1)
    pl.subplot(2,3,3)
    bp.bin_scatter(vparaTrue[:,2],vparaReg[:,2],binsize=0.01,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,2].min(),vparaTrue[:,2].max()],[vparaTrue[:,2].min(),vparaTrue[:,2].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('defocus [mm]')
    pl.xlim(-0.1,0.1)
    pl.ylim(-0.1,0.1)

    pl.subplot(2,3,4)
    bp.bin_scatter(vparaTrue[:,3],vparaReg[:,3],binsize=5,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,3].min(),vparaTrue[:,3].max()],[vparaTrue[:,3].min(),vparaTrue[:,3].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('tilt angle in x direction [arcsec]')
    pl.xlim(-40,40)
    pl.ylim(-40,40)

    pl.subplot(2,3,5)
    bp.bin_scatter(vparaTrue[:,4],vparaReg[:,4],binsize=5,fmt='b.',scatter=scatter)
    pl.plot([vparaTrue[:,4].min(),vparaTrue[:,4].max()],[vparaTrue[:,4].min(),vparaTrue[:,4].max()],'r-')
    pl.xlabel('True Value')
    pl.ylabel('Regression Value')
    pl.title('tilt angle in y direction [arcsec]')
    pl.xlim(-40,40)
    pl.ylim(-40,40)
    return vparaTrue,vparaReg





if __name__ == '__main__':

    #genValidation(500)
    #t=singlemachine_addseeing()
    #computer = sys.argv[1]
    #multimachine_addseeing(computer)
    #multimachine_psfgen(computer)            
    #allfile=gl.glob('*fftconvolve.fit')
    #count = 0
    #for filename in allfile:
    #    count = count +1
    #    print count
    #    t = zernike_file(filename)
    #    pngname = filename[:-3]+'png'
    #    pl.savefig(pngname)
    #    pl.close()
    """
 
    Tfile='/home/jghao/research/decamFocus/psf_withseeing/finerGrid_coeff_matrix/zernike_coeff_finerGrid_training.cp'
    Vfile = '/home/jghao/research/decamFocus/psf_withseeing/finerGrid_coeff_matrix/zernike_coeff_finerGrid_validate.cp'

    # combine files ---
    f1=p.load(open('zernike_coeff_finerGrid_des04.cp','r'))
    f2=p.load(open('zernike_coeff_finerGrid_des05.cp','r'))
    f3=p.load(open('zernike_coeff_finerGrid_des06.cp','r'))
    f4=p.load(open('zernike_coeff_finerGrid_des07.cp','r'))
    f5=p.load(open('zernike_coeff_finerGrid_des08.cp','r'))
    f6=p.load(open('zernike_coeff_finerGrid_des09.cp','r'))
    f7=p.load(open('zernike_coeff_finerGrid_des10.cp','r'))
    f = np.concatenate((f1,f2,f3,f4,f5,f6,f7))
    p.dump(f,open('zernike_coeff_finerGrid_all.cp','w'),2)
    rnd = np.random.rand(len(f[:,1]))
    idx =np.argsort(rnd)
    validate = f[idx[0:500],]
    training = f[idx[500:],]
    p.dump(training,open('zernike_coeff_finerGrid_training.cp','w'),2)
    p.dump(validate,open('zernike_coeff_finerGrid_validate.cp','w'),2)
    """

    t= moments_display(Nstar=1,npix = npix)
    pl.savefig('moments_seeing0.9.png')
    pl.close()
    t= moments_display(Nstar=1,npix = npix,seeing=0)
    pl.savefig('moments_seeing0.png')
    pl.close()
    t= moments_display(Nstar=1,z=0.1,npix = npix)
    pl.savefig('moments_seeing0.9_z0.1.png')
    pl.close()
    t= moments_display(Nstar=1,z=-0.1,npix = npix)
    pl.savefig('moments_seeing0.9_z-0.1.png')
    pl.close()
    t= moments_display(Nstar=1,x=0.1,npix = npix)
    pl.savefig('moments_seeing0.9_x0.1.png')
    pl.close()
    t= moments_display(Nstar=1,x=-0.1,npix = npix)
    pl.savefig('moments_seeing0.9_x-0.1.png')
    pl.close()
    t= moments_display(Nstar=1,y=0.1,npix = npix)
    pl.savefig('moments_seeing0.9_y0.1.png')
    pl.close()
    t= moments_display(Nstar=1,y=-0.1,npix = npix)
    pl.savefig('moments_seeing0.9_y-0.1.png')
    pl.close()
    t= moments_display(Nstar=1,theta=30,phi=0,npix = npix)
    pl.savefig('moments_seeing0.9_theta30_phi0.png')
    pl.close()
    t= moments_display(Nstar=1,theta=-30,phi=0,npix = npix)
    pl.savefig('moments_seeing0.9_theta-30_phi0.png')
    pl.close()
    t= moments_display(Nstar=1,theta=30,phi=90,npix = npix)
    pl.savefig('moments_seeing0.9_theta30_phi90.png')
    pl.close()
    t= moments_display(Nstar=1,theta=-30,phi=90,npix = npix)
    pl.savefig('moments_seeing0.9_theta-30_phi90.png')
    pl.close()
    
    """
    #-----coeff ------
    t= coeff_display(Nstar=1,npix = npix)
    pl.savefig('coeff_seeing0.9.png')
    pl.close()
    t= coeff_display(Nstar=1,npix = npix)
    pl.savefig('coeff_seeing0.png')
    pl.close()
    t= coeff_display(Nstar=1,z=0.1,npix = npix)
    pl.savefig('coeff_seeing0.9_z0.1.png')
    pl.close()
    t= coeff_display(Nstar=1,z=-0.1,npix = npix)
    pl.savefig('coeff_seeing0.9_z-0.1.png')
    pl.close()
    t= coeff_display(Nstar=1,x=0.1,npix = npix)
    pl.savefig('coeff_seeing0.9_x0.1.png')
    pl.close()
    t= coeff_display(Nstar=1,x=-0.1,npix = npix)
    pl.savefig('coeff_seeing0.9_x-0.1.png')
    pl.close()
    t= coeff_display(Nstar=1,y=0.1,npix = npix)
    pl.savefig('coeff_seeing0.9_y0.1.png')
    pl.close()
    t= coeff_display(Nstar=1,y=-0.1,npix = npix)
    pl.savefig('coeff_seeing0.9_y-0.1.png')
    pl.close()
    t= coeff_display(Nstar=1,theta=30,phi=0,npix = npix)
    pl.savefig('coeff_seeing0.9_theta30_phi0.png')
    pl.close()
    t= coeff_display(Nstar=1,theta=-30,phi=0,npix = npix)
    pl.savefig('coeff_seeing0.9_theta-30_phi0.png')
    pl.close()
    t= coeff_display(Nstar=1,theta=30,phi=90,npix = npix)
    pl.savefig('coeff_seeing0.9_theta30_phi90.png')
    pl.close()
    t= coeff_display(Nstar=1,theta=-30,phi=90,npix = npix)
    pl.savefig('coeff_seeing0.9_theta-30_phi90.png')
    pl.close()
    """
