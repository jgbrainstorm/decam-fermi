set startdir [file dirname [info script]]

source $startdir/optlist.tcl
source $startdir/optplot.tcl
source $startdir/poly.tcl
source $startdir/scalecomp.tcl
source $startdir/setup.tcl
source $startdir/distort.tcl
source $startdir/tilt.tcl
source $startdir/opticplot.tcl
source $startdir/baffle.tcl
source $startdir/lstsq.tcl
source $startdir/plot.tcl
source $startdir/read.tcl
source $startdir/glass.tcl
source $startdir/catalog.tcl
source $startdir/maksutov.tcl
source $startdir/diffract.tcl
source $startdir/fft.tcl
source $startdir/snap.tcl
source $startdir/stop.tcl
source $startdir/ray.tcl
source $startdir/waveFront.tcl
source $startdir/intercept.tcl
source $startdir/ee.tcl
source $startdir/region.tcl
source $startdir/fork.tcl
source $startdir/parallel.tcl
source $startdir/tolerance.tcl
source $startdir/lensIO.tcl
source $startdir/lens.tcl
source $startdir/paramRead.tcl
source $startdir/util.tcl
source $startdir/aberr.tcl
source $startdir/ghost.tcl
source $startdir/geom.tcl
source $startdir/refract.tcl
source $startdir/zern.tcl
source $startdir/zlist.tcl
source $startdir/atm.tcl
source $startdir/oslo.tcl
source $startdir/pupil.tcl
source $startdir/syminv.tcl
source $startdir/gradient.tcl
source $startdir/quaternion.tcl

#Handy proc to reload startup scripts.
proc reload {} [list uplevel #0 source $startdir/rayStartup.tcl]

#If ccx program is in path, load the finite element analysis code interfaces.
if {[auto_execok ccx] != "" && [auto_execok triangle] != ""} {
   echo Loading finite element analysis code.
   source $startdir/fea/init.tcl
} else {
#   echo Cannot find ccx and triangle programs in PATH, fea code not loaded.
   }

#Load tkTable.tcl by hand because otherwise Tcl mechanisms are too intricate
#to be sensible
if {[info exists env(TKTABLE_DIR)] && [info command winfo] != ""} {
   source $env(TKTABLE_DIR)/etc/tkTable.tcl
   if {[info command .r] == ""} {
	plotInit a
	lower .r
	}
   }

#Set glassmode - I will use dispersion mode from now on!!!
#Options are "glass" and "disp".  "glass" uses interpolation in opus tables.
global GLASSMODE
set GLASSMODE disp

#Procedure to check or set verbosity.
proc verbose {{val ""}} {
   global VERBOSE
   if {![info exists VERBOSE]} {set VERBOSE 0}
   if {$val == ""} {return $VERBOSE}
   if {![ctype digit $val]} {return 0}
   set VERBOSE $val
   return $VERBOSE
   }

################################
#Procedure to return version of cray
proc version {} {
   global env
   if {[info exists env(SETUP_CRAY)]} {
	set version [lindex $env(SETUP_CRAY) 1]
   } else {
	set version unknown
	}
   return $version
   }

verbose 1

echo Version: [version]

#Lower the priority.  Often, if I start a long CPU-bound activity, it is nice
#to other people on the machine if my process is "niced".

exec renice 10 [id process]
source $startdir/ile.tcl

