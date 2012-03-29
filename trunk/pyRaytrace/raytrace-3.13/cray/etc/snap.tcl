#Construct snap design with full focal plane, starting from initial snap.dat
#file

##################################################################
proc filterZ {n} {
   set z [expr pow(1.17,$n-1) - 1.]
   return $z
   }

#################################################################
proc snapDesign {optic} {

#Arbitrary central wavelength of shortest wave filter.
   set wave0 .44
   setFocal $optic 1 wave .44

#Flatten focal plane.
   setSurf $optic 5 a2 0

#Insert stops
#Struts are 40 mm wide!
#Well, my old strutAdd actually took half-width as input, but now it takes
#full width.
   strutAdd $optic 0 -2000 40 0
   strutAdd $optic 0 -2000 40 120
   strutAdd $optic 0 -2000 40 240

#Clean up the optical design - I entered the wrong conic constant for the
#primary!
   setSurf $optic 5 a6 0
   setSurf $optic 6 a6 0
   setSurf $optic 7 a6 0
   setSurf $optic 5 ccon -0.981128
   setSurf $optic 6 z -2000
   setSurf $optic 7 z 1780
   setSurf $optic 8 z 10

#Fill in filters 1-9 with full field and full SNAP filter complement.
   set scale [showFocal $optic 1 scale]
   set scale -9.54902625
   setFocal $optic 1 scale $scale
   set wave(1) $wave0
   loop i 2 10 {
	set wave($i) [format %.2f [expr .44*(1.+[filterZ $i])]]
	waveAdd $optic $wave($i)
	setFocal $optic $i scale $scale
	}

#Set filter index for each focal plane position by hand - too messy for
#an equation!
   set ifil(11) 1
   set ifil(12) 2
   set ifil(15) 3
   set ifil(16) 4
   set ifil(19) 5
   set ifil(20) 6
   set ifil(13) 2
   set ifil(14) 3
   set ifil(17) 4
   set ifil(18) 5
   set ifil(21) 6
   set ifil(22) 1
   set ifil(23) 3
   set ifil(24) 4
   set ifil(27) 5
   set ifil(28) 6
   set ifil(31) 1
   set ifil(32) 2
   set ifil(25) 4
   set ifil(26) 5
   set ifil(29) 6
   set ifil(30) 1
   set ifil(33) 2
   set ifil(34) 3
   set ifil(35) 5
   set ifil(36) 6
   set ifil(39) 1
   set ifil(40) 2
   set ifil(43) 3
   set ifil(44) 4
   set ifil(37) 6
   set ifil(38) 1
   set ifil(41) 2
   set ifil(42) 3
   set ifil(45) 4
   set ifil(46) 5
   set ifil(47) 7
   set ifil(48) 8
   set ifil(49) 9
   set ifil(50) 9
   set ifil(51) 7
   set ifil(52) 8
   set ifil(53) 8
   set ifil(54) 9
   set ifil(55) 7

#Add the individual detectors in filters 11 to 55
#Number of surfaces now - must fill in optical index for new filters
   set nsurf [exprGet $optic.nsurf]

#Reset id of last surface to its index.  There is no high-level wrapper
#for this (yet).
   handleSet $optic.optic<$nsurf>->idint $nsurf
   set z [showSurf $optic $nsurf z]

#Optical CCDs
   set ccdwidth [expr 287.1e-3*128.]
   set xw1 [expr 12.*287.1e-3]
   set xw2 [expr 64.*287.1e-3]
   set yw [expr 14.*287.1e-3]

echo ccdwidth [format %.1f $ccdwidth]
echo xw1 [format %.1f $xw1]
echo xw2 [format %.1f $xw2]
echo yw [format %.1f $yw]
   set icolor [expr [exprGet $optic.ncolor] + 1]
   loop y 0 3 {
	loop x 0 3 {
	   loop iy 0 2 {
		loop ix 0 2 {
		   set xc [expr (-1.)*(0.5*$xw2 + 3.*$ccdwidth + 2.*$xw1) + \
			$x*($xw1 + $ccdwidth) + ($ix+0.5)*$ccdwidth/2.]
		   set yc [expr 129. + (3.*$ccdwidth + 2.*$yw) - \
			$y*($yw + $ccdwidth) - ($iy+0.5)*$ccdwidth/2.]
		   set phi [expr atan2($yc, $xc)]
		   incr icolor
		   set surf $icolor
		   setSurf $optic $icolor x $xc
		   setSurf $optic $icolor y $yc
		   setSurf $optic $icolor z $z
		   setSurf $optic $icolor phi $phi
		   setIndex $optic $icolor $icolor [showIndex $optic \
			[surfId $optic $nsurf] 1]

#Add wavelength by hand in focal plane - too much custom stuff for waveAdd
#The next statement should be doable via setSurf
		   setIndex $optic 0 $icolor 1
		   setFocal $optic $icolor wave $wave($ifil($icolor))

#Position on sky.  Don't know if my signs are correct yet!
		   set xoff [expr $xc*$scale/60.]
		   set yoff [expr $yc*$scale/60.]
		   setFocal $optic $icolor xoff $xoff
		   setFocal $optic $icolor yoff $yoff
		   setFocal $optic $icolor xsize [expr -$ccdwidth/4.]
		   setFocal $optic $icolor ysize [expr -$ccdwidth/4.]
		   setFocal $optic $icolor scale $scale
		   setFocal $optic $icolor dist 0
		   setFocal $optic $icolor rot 0

#Add new color to existing surfaces (except last one!)
		   loop i 1 $nsurf {
			set isurf [surfId $optic $i]
			set index [showIndex $optic $isurf 1]
			setIndex $optic $isurf $icolor $index
			}
		   opticInfo $optic $icolor
		   }
		}
	   }
	}

#IR arrays
   loop y 0 3 {
	loop x 0 3 {
	   set xc [expr .5*$xw2 + $x*($ccdwidth+$xw1) + 0.5*$ccdwidth]
	   set yc [expr 129. + (3.*$ccdwidth + 2.*$yw) - $y*($ccdwidth+$yw) - \
		0.5*$ccdwidth]
	   set phi [expr atan2($yc, $xc)]
	   incr icolor
	   set surf $icolor
	   setSurf $optic $icolor x $xc
	   setSurf $optic $icolor y $yc
	   setSurf $optic $icolor z $z
	   setSurf $optic $icolor phi $phi
	   setIndex $optic $icolor $icolor [showIndex $optic \
			[surfId $optic $nsurf] 1]

#Add wavelength by hand in focal plane - too much custom stuff for waveAdd
#The next statement should be doable via setSurf
	   setIndex $optic 0 $icolor 1
	   setFocal $optic $icolor wave $wave($ifil($icolor))

#Position on sky.  Don't know if my signs are correct yet!
	   set scale [showFocal $optic 1 scale]
	   set xoff [expr $xc*$scale/60.]
	   set yoff [expr $yc*$scale/60.]
	   setFocal $optic $icolor xoff $xoff
	   setFocal $optic $icolor yoff $yoff
	   setFocal $optic $icolor xsize [expr -$ccdwidth/2.]
	   setFocal $optic $icolor ysize [expr -$ccdwidth/2.]
	   setFocal $optic $icolor scale $scale
	   setFocal $optic $icolor dist 0
	   setFocal $optic $icolor rot 0

#Add new color to existing surfaces (except last one!)
	   loop i 1 $nsurf {
		set isurf [surfId $optic $i]
		set index [showIndex $optic $isurf 1]
		setIndex $optic $isurf $icolor $index
		}
	   opticInfo $optic $icolor
	   }
	}

   colorcount $optic
   stopcomp $optic
   opticinc $optic 1

#Reset field centers for off-axis detectors.  Do a least squares fit.
#   echo Resetting CCD tilts
#   loop i 11 56 {
#	setflag $optic 0 $i 1

#Refocus and tilt?  I will NOT for now so I really do match the SNAP design
#	setflag $optic $i 5 1
#	setflag $optic $i 7 1
#	}
#   lstsq $optic 0
   echo Resetting scale, field centers
   loop i 11 56 {
	setColorFlag $optic $i 1
	setFocalFlag $optic $i 1
	setFocalFlag $optic $i scale 1
	setFocalFlag $optic $i xoff 1
	setFocalFlag $optic $i rot 1
	}
   plstsq $optic 1

   echo Resetting scale for first 9 filters
   loop i 1 10 {
	setflag $optic 0 $i 1
	setFocalFlag $optic $i scale 1
	}
   lstsq $optic 1

   stopcomp $optic
   opticinc $optic 1
   return
   }

###########################################################################
#Run a trace and return the mean incidence angle.
proc incidenceReport {optic xmm ymm icolor} {
   rtrace $optic $xmm $ymm $icolor 0
   set ang [exprGet $optic.diagram->incidence]
   set ang [format %02f [expr $ang*57.3]]
   return $ang
   }

######################################################################
proc snell {ang n} {
   set ang2 [format %.3f [expr 57.3*sin($ang/57.3)*(1./$n)]]
   return $ang2
   }

proc dz {ang n} {
   set ang2 [snell $ang $n]
   set dz [format %.3f [expr 1./cos($ang2/57.3) - 1.]]
   return $dz
   }

#####################################################################
#Layout the spot pattern used by Scholl for tolerance studies

proc snapSpot {} {
   global _optic
   set _optic(spots) ""
   foreach rad ".46 .73 1" {
	lappend _optic(spots) [list $rad 0 1]
	lappend _optic(spots) [list -$rad 0 1]
	lappend _optic(spots) [list 0 $rad 1]
	lappend _optic(spots) [list 0 -$rad 1]
	}
   set _optic(niter) 1
   return
   }
