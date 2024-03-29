#! /usr/bin/env python
#---this is a python scrit that encapsulate the raytracing code from steve to make it more interactive/scriptable.
# the zernike polynomial definition codes is adopted from:
#http://www.staff.science.uu.nl/~werkh108/docs/teach/2011b_python/python102/examples/py102-example2-zernike.py
# J Hao 3/28/2012 @ FNAL

try:
    import numpy as np
    import pyfits as pf
    import pylab as pl
    import os
    from DECamCCD_def import *
    #from sklearn.decomposition import PCA
    #from scipy.optimize import leastsq
    #import mahotas as mh
    import scipy.ndimage as nd
    import healpy as hp
    import glob as gl
    from scipy.misc import factorial as fac
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

def wr(x,y,xcen,ycen,sigma):
    """
    Returns a gaussian weight function with the given parameters
    """
    res=np.exp(-((x-xcen)**2+(y-ycen)**2)/(2.*sigma**2))/(2.*np.pi*sigma**2) 
    return res

def gauss_seeing(npix = None,fwhm=None,e1=None,e2=None):
    """
    generate a seeing PSF of given fwhm and e1 and e2
    """
    scale = 0.27
    fwhm = fwhm/scale
    M20 = (fwhm/2.35482)**2
    row,col = np.mgrid[-npix/2:npix/2+1,-npix/2:npix/2+1]
    rowc = row.mean()
    colc = col.mean()
    Mcc = 0.5*M20*(1+e1)
    Mrc = 0.5*e2*M20
    Mrr = 0.5*M20*(1-e1)
    rho = Mrc/np.sqrt(Mcc*Mrr)
    img = np.exp(-0.5/(1-rho**2)*(row**2/Mrr + col**2/Mcc - 2*rho*row*col/np.sqrt(Mrr*Mcc)))
    res = img/img.sum()
    return res


                 
   
def amoments(data,rowmean=None,colmean=None,sigma=None):
    """
    This codes calcualte the moments with/without a Gaussian weight.
    the sigma is the weight the Gaussian weights. It is the sqrt(sigma_x**2+sigma_y**2)
    col : x direction
    row : y direction
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    IcolSum = np.sum(Icol)
    IrowSum = np.sum(Irow)
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    if rowmean == None:
        rowmean=np.sum(rowgrid*Irow)/IrowSum
        colmean=np.sum(colgrid*Icol)/IcolSum
    if sigma == None:
        rowvar = np.sum((rowgrid-rowmean)**2*Irow)/(IrowSum)
        colvar = np.sum((colgrid-colmean)**2*Icol)/(IcolSum)
    else:
        ROW,COL=np.indices((nrow,ncol))
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        wrcol = wrmat.sum(axis=0)
        wrrow = wrmat.sum(axis=1)
        wrcolsum = np.sum(Icol*wrcol)
        wrrowsum = np.sum(Irow*wrrow)
        rowvar = np.sum((rowgrid-rowmean)**2*Irow*wrrow)/(wrrowsum)
        colvar = np.sum((colgrid-colmean)**2*Icol*wrcol)/(wrcolsum)
        #rowcolcov = np.sum((rowgrid-rowmean)*(colgrid-colmean)*IWmat)/np.sum(IWmat)
        rowcolcov = np.sum(np.outer((rowgrid-rowmean),(colgrid-colmean))*IWmat)/np.sum(IWmat)
    return rowmean,colmean,rowvar,colvar,rowcolcov


def AdaptM(data=None,sigma=None):
    """
    Calculate the adaptive moements e1 and e2 and the variances in row and col
    """
    rowmean,colmean,rowvar,colvar,rowcolcov = amoments(data,sigma=sigma)
    #rowmean,colmean,rowvar,colvar,rowcolcov = amoments(data,rowmean,colmean,sigma=sigma)
    #rowmean,colmean,rowvar,colvar,rowcolcov = amoments(data,rowmean,colmean,sigma=sigma)
    mrrcc = rowvar + colvar
    me1 = (colvar - rowvar)/mrrcc
    me2 = 2.*rowcolcov/mrrcc
    scale = 0.27
    fwhm = np.sqrt(rowvar+colvar)*2.35482*0.27
    return me1,me2,fwhm

def moments7(data=None,sigma=None,rowmean=None,colmean=None):
    """
    This one calcualte the 3 2nd moments and 4 thrid moments with the Gaussian weights.
    col : x direction
    row : y direction
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    IcolSum = np.sum(Icol)
    IrowSum = np.sum(Irow)
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    if rowmean == None:
        rowmean=np.sum(rowgrid*Irow)/IrowSum
        colmean=np.sum(colgrid*Icol)/IcolSum
    if sigma == None:
        rowvar = np.sum((rowgrid-rowmean)**2*Irow)/(IrowSum)
        colvar = np.sum((colgrid-colmean)**2*Icol)/(IcolSum)
    else:
        ROW,COL=np.indices((nrow,ncol))
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        wrcol = wrmat.sum(axis=0)
        wrrow = wrmat.sum(axis=1)
        wrcolsum = np.sum(Icol*wrcol)
        wrrowsum = np.sum(Irow*wrrow)
        Mrr = np.sum((rowgrid-rowmean)**2*Irow*wrrow)/(wrrowsum)
        Mcc = np.sum((colgrid-colmean)**2*Icol*wrcol)/(wrcolsum)
        Mrc = np.sum(np.outer((rowgrid-rowmean),(colgrid-colmean))*IWmat)/np.sum(IWmat)
        Mrrr = np.sum((rowgrid-rowmean)**3*Irow*wrrow)/(wrrowsum)
        Mccc = np.sum((colgrid-colmean)**3*Icol*wrcol)/(wrcolsum)
        Mrrc = np.sum(np.outer((rowgrid-rowmean)**2,(colgrid-colmean))*IWmat)/np.sum(IWmat)
        Mrcc = np.sum(np.outer((rowgrid-rowmean),(colgrid-colmean)**2)*IWmat)/np.sum(IWmat)
        return Mrr, Mcc, Mrc, Mrrr, Mccc, Mrrc, Mrcc

def complexMoments(data=None,sigma=None):
    """
    This one calcualte the 3 2nd moments and 4 thrid moments with the Gaussian weights.
    col : x direction
    row : y direction
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    IcolSum = np.sum(Icol)
    IrowSum = np.sum(Irow)
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    ROW,COL=np.indices((nrow,ncol))
    rowmean=np.sum(rowgrid*Irow)/IrowSum
    colmean=np.sum(colgrid*Icol)/IcolSum
    wrmat = wr(ROW,COL,rowmean,colmean,sigma)
    IWmat = data*wrmat
    wrcol = wrmat.sum(axis=0)
    wrrow = wrmat.sum(axis=1)
    wrcolsum = np.sum(Icol*wrcol)
    wrrowsum = np.sum(Irow*wrrow)
    rowmean = np.sum(rowgrid*Irow*wrrow)/wrrowsum
    colmean = np.sum(colgrid*Icol*wrcol)/wrcolsum
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mr = np.sum(rowgrid*Irow*wrrow)/wrrowsum
    Mc = np.sum(colgrid*Icol*wrcol)/wrcolsum
    Mrr = np.sum(rowgrid**2*Irow*wrrow)/(wrrowsum)
    Mcc = np.sum(colgrid**2*Icol*wrcol)/(wrcolsum)
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/np.sum(IWmat)
    Mrrr = np.sum(rowgrid**3*Irow*wrrow)/(wrrowsum)
    Mccc = np.sum(colgrid**3*Icol*wrcol)/(wrcolsum)
    Mrrc = np.sum(np.outer(rowgrid**2,colgrid)*IWmat)/np.sum(IWmat)
    Mrcc = np.sum(np.outer(rowgrid,colgrid**2)*IWmat)/np.sum(IWmat)
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    M31 = complex(3*Mc - (Mccc+Mrrc)/sigma**2, 3*Mr - (Mrcc + Mrrr)/sigma**2)
    M33 = complex(Mccc-3*Mrrc, 3.*Mrcc - Mrrr)
    return M20, M22, M31, M33



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
seeing=0.9  # in arcseconds, fwhm
#------------------------------

def decamspot(xmm=None,ymm=None,seeing=0.9,npix=40,zenith=0,filter='r', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None):
    #---generating the .par file------
    dir ='/home/jghao/research/ggsvn/decam-fermi/pyRaytrace/'
    file = open(dir+'temp.par','w')
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
    os.system(dir+'raytrace-3.13/decamspot '+dir+'temp.par')
    #---output the result as an image vector
    b=pf.getdata(dir+'temp.fit')
    if seeing !=0:
        seeing = seeing/scale # convert arcseconds to pixel
        sgm = seeing/2.35482  # convert fwhm to sigma
        sgmx = np.sqrt(sgm**2/2.)
        sgmy = sgmx
        b=nd.filters.gaussian_filter(b,(sgmx,sgmy))
    hdr = pf.getheader('temp.fit')
    ypstamp,xpstamp = nd.center_of_mass(b) # y -> row, x-> col
    bb = b.reshape(npix*npix)
    pos = np.array([hdr['xcen'],hdr['ycen'],xpstamp,ypstamp])
    os.system('rm '+dir+'temp.fit temp.par')
    return np.concatenate((pos,bb)),hdr


def genImgV(filename=None,Nstar=None,ccd=None,seeing=0,npix=40,zenith=0,filter='g', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None,regular=False):
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
            res = decamspot(xmm=xmm,ymm=ymm,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z,suband=suband)
            datalist.append(res[0])
            hdrlist.append(res[1])
        data = np.array(datalist)
    if filename is not None:
        hdu = pf.PrimaryHDU(data)
        hdu.header.update('RAYPATT',raypattern)
        hdu.header.update('NPIX',npix)
        hdu.header.update('SCALE',scale)
        hdu.header.update('FWHM',fwhm)
        hdu.header.update('ZENITH',zenith)
        hdu.header.update('FILTER',filter)
        hdu.header.update('THETA',theta)
        hdu.header.update('CORRT',corrector)
        if x != None:
            hdu.header.update('X',x)
        if y != None:
            hdu.header.update('Y',y)
        if z != None:
            hdu.header.update('Z',z)
        hdu.header.update('Seeing',seeing)
        if os.path.exists(filename):
            os.system('rm '+filename)
            hdu.writeto(filename)
        else:
            hdu.writeto(filename)        
    return data,hdrlist


def genImgVfixedPos(filename=None,seeing=0,npix=40,zenith=0,filter='r', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None):
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
        hdu.header.update('RAYPATT',raypattern)
        hdu.header.update('NPIX',npix)
        hdu.header.update('SCALE',scale)
        hdu.header.update('FWHM',fwhm)
        hdu.header.update('ZENITH',zenith)
        hdu.header.update('FILTER',filter)
        hdu.header.update('THETA',theta)
        hdu.header.update('CORRT',corrector)
        if x != None:
            hdu.header.update('X',x)
        if y != None:
            hdu.header.update('Y',y)
        if z != None:
            hdu.header.update('Z',z)
        hdu.header.update('Seeing',seeing)
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

  
def imgCCDctr(ccd=None,filter='g',seeing=0,z=0.,theta=0.,contour=False):
    xmm = ccd[1]
    ymm = ccd[2]
    res = genImgV(Nstar=1, ccd = ccd,seeing=seeing,theta=theta,z=z)
    data = res[0]
    disImgCCD(data,ccd)
    if contour is True:
        pl.contourf(data[0][2:].reshape(npix,npix),n=100)
    e1,e2, rowvar,colvar =AdaptM(data[0][2:].reshape(npix,npix),sigma=1.1)
    xcen = res[1][0]['xcen']
    ycen = res[1][0]['ycen']
    pl.figtext(0.2,0.84,'CCD: '+ccd[0] +',   '+'Filter: '+filter, color='r')
    pl.figtext(0.2,0.8, 'e1: '+str(round(e1,3)) + ',  e2: '+str(round(e2,3)), color='r')
    pl.figtext(0.2,0.75, 'xcen: '+str(xcen) + ',  ycen: '+str(ycen), color='r')
    pl.figtext(0.2,0.7, 'seeing: '+str(seeing)+'" (fwhm)', color='r')
    if z:
        pl.figtext(0.2,0.65, 'defocus: '+str(z)+' (mm)', color='r')
    if theta:
        pl.figtext(0.2,0.6,'tilt: '+str(theta)+' (arcsec)', color='r')
    return e1,e2, rowvar,colvar

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

    
def decompZernike(ccd=None,seeing=0.,theta=0.,x=0.,y=0.,z=0.,filter='g'):
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

def psfSizeAll(Nstar=None,filter='r',npix=40,seeing=0,theta=0., zenith = 0.,corrector='corrector', x=None, y=None,z=None):
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

def psfSizeAllZernike(Nstar=None,filter='r',npix=40,seeing=0,theta=0., zenith = 0.,corrector='corrector', x=None, y=None,z=None,rand=False):
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


def genPSFimage(filename=None,dir=None):
    """
    convert the PSF image vector file to the set of PSF images
    """
    b=pf.getdata(filename)
    nn = len(b)
    npix = int(np.sqrt(len(b[0][2:])))
    for i in range(nn):
        img = b[i][2:].reshape(npix,npix)
        img = img/img.sum()
        h = pf.PrimaryHDU(img)
        h.header.update('xmm',b[i][0])
        h.header.update('ymm',b[i][1])
        h.scale('int16', '', bzero=32768)
        h.writeto(dir+'psf_'+str(i)+'.fit')


def measuredata(filename):
    b=pf.getdata(filename)
    Nobj = len(b)
    x=np.zeros(Nobj)
    y=np.zeros(Nobj)
    e1=np.zeros(Nobj)
    e2=np.zeros(Nobj)
    rowvar=np.zeros(Nobj)
    colvar=np.zeros(Nobj)
    colnames = ['x','y','e1','e2','rowvar','colvar']
    sigma = 1.1/0.27
    for i in range(Nobj):
        e1[i],e2[i],rowvar[i],colvar[i]=AdaptM(b[i][2:].reshape(40,40),sigma=sigma)
        x[i]=b[i][0]
        y[i]=b[i][1]
        data = [x,y,e1,e2,rowvar,colvar]
    hp.mwrfits(filename[:-4]+'_moments_gausswt_11.fit',data,colnames=colnames)
    return 0


def measuredataM7(filename):
    b=pf.getdata(filename)
    Nobj = len(b)
    x=np.zeros(Nobj)
    y=np.zeros(Nobj)
    Mrr=np.zeros(Nobj)
    Mcc=np.zeros(Nobj)
    Mrc=np.zeros(Nobj)
    Mrrr=np.zeros(Nobj)
    Mccc=np.zeros(Nobj)
    Mrrc=np.zeros(Nobj)
    Mrcc=np.zeros(Nobj)
    
    colnames = ['x','y','Mrr','Mcc','Mrc','Mrrr','Mccc','Mrrc','Mrcc']
    sigma = 1.1/0.27
    for i in range(Nobj):
        Mrr[i],Mcc[i],Mrc[i],Mrrr[i],Mccc[i],Mrrc[i],Mrcc[i]=moments7(b[i][4:].reshape(40,40),sigma=sigma)
        x[i]=b[i][0]
        y[i]=b[i][1]
        data = [x,y,Mrr, Mcc, Mrc, Mrrr, Mccc, Mrrc, Mrcc]
    hp.mwrfits(filename[:-4]+'_moments7_gausswt_11.fit',data,colnames=colnames)
    return 0


def measureDataM7_multiext(filename):
    hdu=pf.open(filename)
    nn = len(hdu)
    data = []
    colnames = ['x','y','Myy','Mxx','Myx','Myyy','Mxxx','Myyx','Myxx']
    #colnames = ['x','y','Mrr','Mcc','Mrc','Mrrr','Mccc','Mrrc','Mrcc']
    for hdui in hdu[1:]:
        Nobj = hdui.data.shape[0]
        Mrr=np.zeros(Nobj)
        Mcc=np.zeros(Nobj)
        Mrc=np.zeros(Nobj)
        Mrrr=np.zeros(Nobj)
        Mccc=np.zeros(Nobj)
        Mrrc=np.zeros(Nobj)
        Mrcc=np.zeros(Nobj)
        sigma = 1.1/0.27
        for i in range(Nobj):
            Mrr[i],Mcc[i],Mrc[i],Mrrr[i],Mccc[i],Mrrc[i],Mrcc[i]=moments7(hdui.data[i][4:].reshape(40,40),sigma=sigma)
        x=hdui.header['ccdXcen']
        y=hdui.header['ccdYcen']
        data.append([x,y,np.median(Mrr), np.median(Mcc), np.median(Mrc), np.median(Mrrr), np.median(Mccc),np.median(Mrrc), np.median(Mrcc)])
    data=np.array(data)
    hp.mwrfits(filename[:-7]+'_moments7_gausswt_11.fit',data.T,colnames=colnames)
    return '---done !-----'


def measureDataComplexM_multiext(filename):
    hdu=pf.open(filename)
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
    hp.mwrfits(filename[:-7]+'_complexMoments_gausswt_'+str(sigma*0.27)+'.fit',data.T,colnames=colnames)
    return '---done !-----'





def genImgVallCCD(filename=None,Nstar=None,seeing=0,npix=40,zenith=0,filter='g', theta=0., corrector='corrector',x=0.,y=0.,z=0.,suband=None,regular=False):
    """
    Nstar is the number of stars on each CCD
    """
    hduList = pf.HDUList()
    hdu = pf.PrimaryHDU(np.array([0]))
    hdu.header.update('RAYPATT',raypattern)
    hdu.header.update('NPIX',npix)
    hdu.header.update('SCALE',scale)
    hdu.header.update('FWHM',fwhm)
    hdu.header.update('ZENITH',zenith)
    hdu.header.update('FILTER',filter)
    hdu.header.update('THETA',theta)
    hdu.header.update('CORRT',corrector)
    if x != None:
        hdu.header.update('X',x)
    if y != None:
        hdu.header.update('Y',y)
    if z != None:
        hdu.header.update('Z',z)
    hdu.header.update('Seeing',seeing)
    hduList.append(hdu)
    for ccd in N[1:]+S[1:]:
        print ccd
        res = genImgV(filename=filename,Nstar=Nstar,ccd=ccd,seeing=seeing,npix=npix,zenith=zenith,filter=filter, theta=theta, corrector=corrector,x=x,y=y,z=z,suband=suband,regular=regular)
        hdu = pf.PrimaryHDU(res[0])
        hdu.header.update('ccdPos',ccd[0])
        hdu.header.update('ccdXcen',ccd[1])
        hdu.header.update('ccdYcen',ccd[2])
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



def dispZernike(beta=1.,j=0,gridsize = 10, max_rad = 10):
    x,y = np.meshgrid(np.arange(-gridsize,gridsize,0.01),np.arange(-gridsize,gridsize,0.01))
    rho = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    ok = rho < max_rad
    znk = beta*zernikel(j,rho,phi)*ok
    pl.imshow(znk)
    return znk
    
def showZernike(beta=None,gridsize = 1, max_rad = 1):
    x,y = np.meshgrid(np.arange(-gridsize,gridsize,0.01),np.arange(-gridsize,gridsize,0.01))
    rho = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    ok = rho < max_rad
    nn = len(beta)
    znk=0
    for j in range(nn):
        znk = znk + beta[j]*zernikel(j,rho,phi)*ok
    pl.imshow(znk)
    return znk


def zernike_diagnosis(Nstar=None,seeing=0,npix=40,zenith=0,filter='g', theta=0., corrector='corrector',x=None,y=None,z=None,zernike_max_order=20,regular=False):
    """
    This function produce the zernike plots for a set of given parameters of the tilt/shift/defocus
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

#----generate the data ------
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar1000_notilt_nodefocus.fit',Nstar=1000,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_notilt_z0.01mm.fit',Nstar=1000,z=0.01,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_notilt_z0.1mm.fit',Nstar=1000,z=0.1,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_tilt_10arcsec_nodefocus.fit',Nstar=1000,theta=10.,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_tilt_100arcsec_nodefocus.fit',Nstar=1000,theta=100.,seeing=0.9)

    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_tilt_xshift_0.01mm_nodefocus_notilt.fit',Nstar=1000,theta=0.,x=0.01,seeing=0.9)
    genImgV(filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_nstar_1000_tilt_xshift_0.1mm_nodefocus_notilt.fit',Nstar=1000,theta=0.,x=0.1,seeing=0.9)
  

    
    #-----analyzing the generated data -------------
 
 
    filenameAll = gl.glob('/home/jghao/research/decamFocus/PSF_seeing_*xshift*.fit')
    Nfile=len(filenameAll)

    for j in range(Nfile):
        measuredataM7(filenameAll[j])
    
    """

    #--generate multi extension data -----
    filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_rband_nstar20_notilt_nodefocus_multi_ext.fit'
    genImgVallCCD(filename=filename,Nstar=20,seeing=0.9,npix=40,zenith=0,filter='r', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None)
    filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_rband_nstar20_tilt_100arecse_nodefocus_multi_ext.fit'
    genImgVallCCD(filename=filename,Nstar=20,seeing=0.9,npix=40,zenith=0,filter='r', theta=100., corrector='corrector',x=None,y=None,z=None,suband=None)
    #-------one psf per ccd-----
    filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_rband_nstar1_notilt_nodefocus_multi_ext.fit'
    genImgVallCCD(filename=filename,Nstar=1,seeing=0.9,npix=40,zenith=0,filter='r', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None)
    filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_rband_nstar1_tilt_100arecse_nodefocus_multi_ext.fit'
    genImgVallCCD(filename=filename,Nstar=1,seeing=0.9,npix=40,zenith=0,filter='r', theta=100., corrector='corrector',x=None,y=None,z=None,suband=None)

    #--generate multi extension data with regular grid on each CCD -----
    filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_rband_nstar16_regular_notilt_nodefocus_multi_ext.fit'
    genImgVallCCD(filename=filename,Nstar=16,seeing=0.9,npix=40,zenith=0,filter='r', theta=0., corrector='corrector',x=None,y=None,z=None,suband=None,regular=True)
    filename='/home/jghao/research/decamFocus/PSF_seeing_0.9_rband_nstar16_regular_tilt_100arecse_nodefocus_multi_ext.fit'
    genImgVallCCD(filename=filename,Nstar=16,seeing=0.9,npix=40,zenith=0,filter='r', theta=100., corrector='corrector',x=None,y=None,z=None,suband=None,regular=True)
