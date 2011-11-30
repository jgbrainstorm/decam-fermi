#!/usr/bin/env python

# this script run the test motion of the xy stage for the flatness measurement
# choose the correct region and change the position of the projector accordingly.
#For questions, contact J. Hao: jghao@fnal.gov

import sys
import os
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
import xymove as mv


#---upper left--------

#ccd=['s30','s31','s28','s27','s22','s23','s24','s19','s18','s17','s11','s12','s13','s7','s6','s5','s4','n4','n5','n6','n7']

#---lower left--------

ccd=['s7','s6','s5','s4','n4','n5','n6','n7','n13','n12','n11','n17','n18','n19','n24','n23','n22','n27','n28','n31','n30']


#----upper right ----

#ccd=['s30','s29','s25','s26','s22','s21','s20','s14','s15','s16','s10','s9','s8','s1','s2','s3','s4','n4','n3','n2','n1']

#----lower right ----

#ccd=['s1','s2','s3','s4','n4','n3','n2','n1','n8','n9','n10','n16','n15','n14','n20','n21','n22','n26','n25','n29','n30']

from DECamCCD import *

for position in ccd:
    print '----now at: '+position+'  -----'
    i=20
    mv.moveto(position)
    os.system('sleep 5')
       

