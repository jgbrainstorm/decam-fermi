# this codes compare the PTC for several measurements

from DECamCCD import *

def getFileNameList(expid_start,expid_end,directory,skip=1):
    expid = np.arange(expid_start,expid_end+1,skip)
    fileNameList=[]
    for eid in expid:
        fileNameList.append(directory+'DECam_00'+str(eid)+'.fits.fz')
    return fileNameList



def ptc(NameFits,NameBias,Channel,shift=None,left=None,overplot=False,pch='o',color='blue',label='date of ptc'):
    detector = pf.open(NameBias)[Channel].header['DETPOS']   
    if detector[0]=='N':
        if left == None or left == 1:
            colmin=300
            colmax=800
            rowmin=3946 
            rowmax=4046
        else:
            colmin=1360
            colmax=1860
            rowmin=3946
            rowmax=4046
    if detector[0]=='S':
        if left == None or left == 1:
            colmin=1360
            colmax=1860
            rowmin=100  
            rowmax=200
        else:
            colmin=300
            colmax=800
            rowmin=100
            rowmax=200 
    if detector[0]=='F':
        if detector[-1]=='N':
            if left == None or left == 1:
                colmin = 300
                colmax = 800
                rowmin = 1850
                rowmax = 1950
            else:
                colmin = 1360
                colmax = 1860
                rowmin = 1850
                rowmax = 1950
        if detector[-1]=='S':
            if left == None or left == 1:
                colmin = 1360
                colmax = 1860
                rowmin=100
                rowmax=200
            else:
                colmin = 300
                colmax = 800
                rowmin = 100
                rowmax = 200
    if shift == 1:
        colshift = np.random.random_integers(0,500,1)
        rowshift = np.random.random_integers(0,3000,1)
        rowmin = rowmin + rowshift
        rowmax = rowmax + rowshift
        colmin = colmin + colshift
        colmax = colmax + colshift
    NFile = len(NameFits)
    num = int(round((NFile-1)/2.))
    mean_b = np.zeros(num)
    var_b = np.zeros(num)
    exptime=np.zeros(num)
    for i in range(0,num):
        bias,biashdr = readCCDFits(NameBias,Channel)
        bias = bias[rowmin:rowmax,colmin:colmax]
        ba,hdra = readCCDFits(NameFits[i*2],Channel)
        ba = ba[rowmin:rowmax,colmin:colmax]
        ba = ba - bias
        bb,hdrb = readCCDFits(NameFits[i*2+1],Channel)
        bb = bb[rowmin:rowmax,colmin:colmax]
        bb = bb - bias
        uni = np.unique(ba)
        add_b = ba + bb
        diff_b= ba - bb
        mean_b[i] = robust_mean(add_b)/2.
        var_b[i] = robust_var(diff_b)/2.
        exptime[i]=pf.open(NameFits[i*2])[0].header['expreq']        
    sidx=np.argsort(exptime)
    exptime=exptime[sidx]
    mean_b=mean_b[sidx]
    var_b=var_b[sidx]
    ok = (mean_b > 5000)*(mean_b <20000)*(var_b < 8000)          
    if len(mean_b[ok]) > 2:
        (a,b,SEa,SEb,R2) = linefit(mean_b[ok],var_b[ok])
        gain = b
        (ta,tb,tSEa,tSEb,tR2) = linefit(exptime[ok],mean_b[ok])           
    else:
        gain=999
        a=0
        b=0
        SEb=999
        ta=0
        tb=0
        tSEa=999
        tSEb=999       
    if overplot == False:
        pl.figure(figsize=(10,10))
    if left == None:
        pl.title('Position:' +detector+'     Channel: '+str(Channel))
    if left == 1:
        pl.title('Position:' +detector+'     Channel: '+str(Channel)+' left')
    if left == 0:
        pl.title('Position:' +detector+'     Channel: '+str(Channel)+' right')
    dff=(mean_b-(tb*exptime+ta))/(tb*exptime+ta)
    uuu=mean_b[(abs(dff) > 0.01)*(mean_b > 10000.)]
    if uuu.any() == True:
        ffwid=np.argwhere(mean_b==min(uuu))[0][0]
        ffwid=ffwid-1
        ffw=mean_b[ffwid]
    else:
        ffw=999
        ffwid=0
    pl.plot(mean_b,dff,pch,color=color,label=label)
    pl.plot(mean_b[ffwid],dff[ffwid],pch,color='red')
    pl.hlines(0.01,min(mean_b),max(mean_b),colors='r',linestyles='dashed')
    pl.hlines(-0.01,min(mean_b),max(mean_b),colors='r',linestyles='dashed')
    pl.ylim(-0.05,0.05)
    pl.xlim(min(mean_b),max(mean_b))
    pl.xlabel('Mean ADU)')
    pl.ylabel('Fractional deviation from the fitted line')
    if type(ffw)==float:
        pl.figtext(0.1,0.9,'Fullwell:'+str(round(ffw))+' (ADU)')
    return 0
        

if __name__ == "__main__":
    
   dir_9_27_12 = '/blue-orch/DTS/src/20120927/src/'
   namebias_9_27_12=dir_9_27_12+'DECam_00136715.fits.fz'
   namefits_9_27_12 = getFileNameList(136717,136858,dir_9_27_12)

   dir_10_2_12 = '/blue-orch/DTS/src/10121002/src/'
   namebias_10_2_12=dir_10_2_12+'DECam_00137808.fits.fz'
   namefits_10_2_12 = getFileNameList(137811,137880,dir_10_2_12)

   dir_10_19_12 = '/blue-orch/DTS/src/20121019/src/'
   namebias_10_19_12=dir_10_19_12+'DECam_00140391.fits.fz'
   namefits_10_19_12 = getFileNameList(140394,140515,dir_10_19_12)


   t=ptc(namefits_9_27_12,namebias_9_27_12,1,left=1,label='9/27/2012')
   t=ptc(namefits_10_2_12,namebias_10_2_12,1,left=1,overplot=True,pch='*',color='blue',label='10/2/2012')
   t=ptc(namefits_10_19_12,namebias_10_19_12,1,left=1,overplot=True,pch='+',color='blue',label='10/19/2012')
