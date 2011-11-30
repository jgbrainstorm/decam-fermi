#!/usr/bin/env python

# This script run the flatness measurement. You need to choose one quandrant each time and mount the projector at appropriate place accordingly.

#For questions, contact J. Hao: jghao@fnal.gov

import sys
import os
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
import xymove as mv

#---upper left--------

#ccd=['s30','s31','s28','s27','s22','s23','s24','s19','s18','s17','s11','s12','s13','s7','s6','s5','s4','n4','n5','n6','n7']

#---lower left--------

#ccd=['s7','s6','s5','s4','n4','n5','n6','n7','n13','n12','n11','n17','n18','n19','n24','n23','n22','n27','n28','n31','n30']


#----upper right ----

ccd=['s30','s29','s25','s26','s22','s21','s20','s14','s15','s16','s10','s9','s8','s1','s2','s3','s4','n4','n3','n2','n1']

#----lower right ----

#ccd=['s1','s2','s3','s4','n4','n3','n2','n1','n8','n9','n10','n16','n15','n14','n20','n21','n22','n26','n25','n29','n30']




import sys, time
import PML.core
import SISPIlib.logger as Log
from fermiMCCD import *
try:
    # Connect to OCS
    ocs = PML.core.PML_Connection('OCS','OCS')

    # Use the get state function to retrieve the state of the OCS state machine
    # Configure if OCS is in state 'ALIVE'
    if ocs('get','state') == 'INITIALIZED':
        ocs('configure','')
        # make sure the OCS is done
        while ocs('get','state') != 'READY':
          time.sleep(1)

    # Are we ready to run our tests?
    reported_ocs_state = ocs('get','state')
    if reported_ocs_state != 'READY':
        # No - something is wrong - we better exit.
        Log.error("OCS is not ready (reported state: %s" % reported_ocs_state)
        sys.exit()
    else:
        Log.info("OCS is ready.")

    # Now we can do something.
    # For example take 10 exposures with increasing exposure time
    #for i in range(1,11):
        # Use the expose function with wait=True to wait until the exposure is completely done and the OCS is
        # back in READY. Note that all PML arguments have to be passed as string
    #---take one bias image--------
    i=20
    response = ocs('expose','exptype=dark, wait=True,exptime = %d' % i)
    os.system('sleep 5')
    Log.info(" Exposure %d: %s" %(i,str(response)))
    os.system('sleep 5')
    for position in ccd:
        print '----now at: '+position+'  -----'
        i=20
        mv.moveto(position)
        os.system('sleep 20')
        response = ocs('expose','wait=True, exptime = %d' % i)
        os.system('sleep 5')
        Log.info(" Exposure %d: %s" %(i,str(response)))
        os.system('sleep 5') 
        
    # You don't have to use the SISPI logger. print works just fine.
    # However, using the logger gives you formatted message including a time stamp and an
    # archive, ie a log file. Log files are kept in your instance directory under logs.
    print "Done"
except KeyboardInterrupt, msg:
    pass
