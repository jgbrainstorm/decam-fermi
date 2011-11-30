"""
This code define the CCD coordinate

"""
import os
from DECamCCD_def import *

def moveto(CCD,xoffset=None,yoffset=None):

    if xoffset == None and yoffset == None:
        xoffset = 0.0
        yoffset = 0.0

#-----South side------------
    S1=["S1",16.908,191.670]
    S2=["S2",16.908,127.780]
    S3=["S3",16.908,63.890]
    S4=["S4",16.908,0.]
    S5=["S5",16.908,-63.890]
    S6=["S6",16.908,-127.780]
    S7=["S7",16.908,-191.670]
    
    S8=["S8",50.724,159.725]
    S9=["S9",50.724,95.835]
    S10=["S10",50.724,31.945]
    S11=["S11",50.724,-31.945]
    S12=["S12",50.724,-95.835]
    S13=["S13",50.724,-159.725]
    
    S14=["S14",84.540,159.725]
    S15=["S15",84.540,95.835]
    S16=["S16",84.540,31.945]
    S17=["S17",84.540,-31.945]
    S18=["S18",84.540,-95.835]
    S19=["S19",84.540,-159.725]
    
    S20=["S20",118.356,127.780]
    S21=["S21",118.356,63.890]
    S22=["S22",118.356,0.]
    S23=["S23",118.356,-63.890]
    S24=["S24",118.356,-127.780]
    
    S25=["S25",152.172,95.835]
    S26=["S26",152.172,31.945]
    S27=["S27",152.172,-31.945]
    S28=["S28",152.172,-95.835]
    
    S29=["S29",185.988,63.890]
    S30=["S30",185.988,0.]
    S31=["S31",185.988,-63.890]
    
    
#-----North side------------

    N1=["N1",-16.908,191.670]
    N2=["N2",-16.908,127.780]
    N3=["N3",-16.908,63.890]
    N4=["N4",-16.908,0.]
    N5=["N5",-16.908,-63.890]
    N6=["N6",-16.908,-127.780]
    N7=["N7",-16.908,-191.670]

    N8=["N8",-50.724,159.725]
    N9=["N9",-50.724,95.835]
    N10=["N10",-50.724,31.945]
    N11=["N11",-50.724,-31.945]
    N12=["N12",-50.724,-95.835]
    N13=["N13",-50.724,-159.725]

    N14=["N14",-84.540,159.725]
    N15=["N15",-84.540,95.835]
    N16=["N16",-84.540,31.945]
    N17=["N17",-84.540,-31.945]
    N18=["N18",-84.540,-95.835]
    N19=["N19",-84.540,-159.725]

    N20=["N20",-118.356,127.780]
    N21=["N21",-118.356,63.890]
    N22=["N22",-118.356,0.]
    N23=["N23",-118.356,-63.890]
    N24=["N24",-118.356,-127.780]

    N25=["N25",-152.172,95.835]
    N26=["N26",-152.172,31.945]
    N27=["N27",-152.172,-31.945]
    N28=["N28",-152.172,-95.835]

    N29=["N29",-185.988,63.890]
    N30=["N30",-185.988,0.]
    N31=["N31",-185.988,-63.890]

#------------Guiding CCD-----------
    GS2=["G2S",185.988,111.8075]
    GS1=["G1S",152.172,161.7525]
    GN2=["G2N",-185.988,111.8075]
    GN1=["G1N",-152.172,161.7525]

    S=["S",S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23,S24,S25,S26,S27,S28,S29,S30,S31]
    N=["N",N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31]

    if CCD=='g1n' or CCD=='G1N':
        x_target = GN1[2]+xoffset
        y_target = GN1[1]+yoffset
    elif CCD=='g2n' or CCD=='G2N':
        x_target = GN2[2]+xoffset
        y_target = GN2[1]+yoffset
    elif CCD=='g1s' or CCD=='G1S':
        x_target = GS1[2]+xoffset
        y_target = GS1[1]+yoffset
    elif CCD=='g2s' or CCD=='G2S':
        x_target = GS2[2]+xoffset
        y_target = GS2[1]+yoffset
    elif CCD[0]=='s' or CCD[0]=='S':
        idx=int(CCD[1:])
        x_target = S[idx][2]+xoffset
        y_target = S[idx][1]+yoffset
    elif CCD[0]=='n' or CCD[0]=='N':
        idx=int(CCD[1:])
        x_target = N[idx][2]+xoffset
        y_target = N[idx][1]+yoffset 
    else:
        print "----Not in the range! Try again!-----"

    os.system('/usr/remote/user/sispi/jiangang/decam-fermi/sendsockcmd -h 139.229.14.66 -p 2055 -t 1000000 "ma,'+str(x_target)+','+str(y_target)+'"')
    return(x_target,y_target)    
   
#------move to left channel-----

def moveto_l(CCD,xoffset=None,yoffset=None):
    if xoffset == None and yoffset == None:
        xoffset = 0.0
        yoffset = 0.0

    (x_target,y_target) = moveto(CCD,xoffset,yoffset)
    if CCD[0]=='s' or CCD[0]=='S':        
        y_target = y_target - 9.0
    else:
        y_target = y_target + 9.0

    os.system('/usr/remote/user/sispi/jiangang/decam-fermi/sendsockcmd -h 139.229.14.66 -p 2055 -t 1000000 "ma,'+str(x_target)+','+str(y_target)+'"')
    return(x_target,y_target)    

  

#-------move to right channel----

def moveto_r(CCD,xoffset=None,yoffset=None):
    if xoffset == None and yoffset == None:
        xoffset = 0.0
        yoffset = 0.0

    (x_target,y_target) = moveto(CCD,xoffset,yoffset)
    if CCD[0]=='s' or CCD[0]=='S':        
        y_target = y_target + 9.0
    else:
        y_target = y_target - 9.0

    os.system('/usr/remote/user/sispi/jiangang/decam-fermi/sendsockcmd -h 139.229.14.66 -p 2055 -t 1000000 "ma,'+str(x_target)+','+str(y_target)+'"')
   
    return(x_target,y_target)      

#-------move to right channel----

def moveto_origin():
    os.system('/usr/remote/usr/sispi/jiangang/decam-fermi/sendsockcmd -h 139.229.14.66 -p 2055 -t 1000000 "ma,0,0"')
    return("moved to origin")


def moveto_xy(x,y):
    os.system('/usr/remote/user/sispi/jiangang/decam-fermi/sendsockcmd -h 139.229.14.66 -p 2055 -t 1000000 "ma,'+str(x)+','+str(y)+'"')
    return("moved to the position you specified")
