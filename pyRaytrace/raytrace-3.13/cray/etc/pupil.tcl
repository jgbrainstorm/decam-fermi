#Compute tweaks to ray aiming parameters to correct for displacement
#and distortion of entrance pupil.

############################################################################
#Get aperture stop radius

proc apRadius {hndl ifil} {
   set surfids [surfIdsGet $hndl $ifil]

#Get the aperture stop
   set iapp -1
   loop i 0 [llength $surfids] {
	set surfid [lindex $surfids $i]
	if {[showSurf $hndl $surfid stoptype] == 2} {
	   set iapp $i
	   set apid $surfid
	   break
	   }
	}
   if {$iapp < 0} {
	error "No aperture stop found for filter $ifil"
	}
   set astop [showSurf $hndl $apid outstop]
   return $astop
   }

############################################################################
#All tweaks can be accommodated by adjusting xfract, yfract in "ray".

#The nominal position on the entract pupil is xfract * outstop (aperture stop)
# * EMAG, etc.

#The most useful remapping probably involves a recenter, a scaling,
#an ellipticity, and a position angle (5 numbers).

proc entranceMap {hndl xmm ymm ifil} {
   set surfids [surfIdsGet $hndl $ifil]

#Get the aperture stop
   set iapp -1
   loop i 0 [llength $surfids] {
	set surfid [lindex $surfids $i]
	if {[showSurf $hndl $surfid stoptype] == 2} {
	   set iapp $i
	   set apid $surfid
	   break
	   }
	}
   if {$iapp < 0} {
	error "No aperture stop found for filter $ifil"
	}
   set astop [showSurf $hndl $apid outstop]
   if {$astop == 0.} {
	error "Aperture stop has radius 0."
	}
   set emag [showFocal $hndl $ifil emag]

#First, see if we can run a ray and hit the aperture stop.
#Search in radius in units of fractions.
   set xcen 0.
   set ycen 0.
   set rad [expr sqrt(pow($xmm,2) + pow($ymm,2))]
   if {$rad == 0} {
	set rad 1
	}
   set xrad [expr $xmm/$rad]
   set yrad [expr $ymm/$rad]
   set flag 0
   set r2min 1.e20
   set xcenmin 0.
   set ycenmin 0.
   set ipass 0
   for {set fract 0} {$fract <= 1.e3} {set fract [expr $fract + .5]} {
	set xcen [expr $fract*$xrad]
	set ycen [expr $fract*$yrad]
	set flag 0
	if {[ray $hndl $xmm $ymm [expr 0.-$xcen] [expr 0.-$ycen] $ifil 0] \
	   == 1} {
	   set flag 1
	   set ipass 1
	   set xapp [exprGet $hndl.diagram->xlens<$iapp>]
	   set yapp [exprGet $hndl.diagram->ylens<$iapp>]
	   set r2 [expr pow($xapp,2) + pow($yapp,2)]
	   if {$r2 < $r2min} {
		set r2min $r2
		set xcenmin $xcen
		set ycenmin $ycen
	   } else {
		break
		}
	   }
	set xcen [expr -1.*$fract*$xrad]
	set ycen [expr -1.*$fract*$yrad]
	if {[ray $hndl $xmm $ymm [expr 0.-$xcen] [expr 0.-$ycen] $ifil 0] \
	   == 1} {
	   set flag 1
	   set ipass 1
	   set xapp [exprGet $hndl.diagram->xlens<$iapp>]
	   set yapp [exprGet $hndl.diagram->ylens<$iapp>]
	   set r2 [expr pow($xapp,2) + pow($yapp,2)]
	   if {$r2 < $r2min} {
		set r2min $r2
		set xcenmin $xcen
		set ycenmin $ycen
	   } else {
		break
		}
	   }
	if {$flag == 0 && $ipass == 1} break
	}
   if {$ipass == 0} {
	error "Chief ray does not reach aperture stop"
	}

   set xcen $xcenmin
   set ycen $ycenmin
   set xcen0 $xcenmin
   set ycen0 $ycenmin
   set xapp [exprGet $hndl.diagram->xlens<$iapp>]
   set yapp [exprGet $hndl.diagram->ylens<$iapp>]

#Run a chief ray and see where it lands.  xcen and ycen will be subtracted
#from nominal xfract, yfract to get actual values needed to land in desired
#place on aperture stop.
   set ntry 0
   set fract 1.
   while {1} {
	if {[ray $hndl $xmm $ymm [expr 0.-$xcen] [expr 0.-$ycen] $ifil 0] \
	   < 1} {
	   set fract [expr $fract/2.]
#echo Shrink fract $fract
	   set xcen $xcen0
	   set ycen $ycen0
	   incr ntry
	   if {$ntry >= 100} {
		error "Unable to zero in chief ray"
		}
	   continue
	   }
	set xapp0 $xapp
	set yapp0 $yapp
	set xapp [exprGet $hndl.diagram->xlens<$iapp>]
	set yapp [exprGet $hndl.diagram->ylens<$iapp>]
	if {abs($xapp) < 1.e-4 && abs($yapp) < 1.e-4} break

#For now, assume that desired xcen, ycen are in same direction as
#xapp, yapp.  I would need to run more rays to get a better feel for
#orientation.

	if {$xcen != $xcen0 && $xapp != $xapp0} {
	   set slope [expr -($xcen-$xcen0)/($xapp-$xapp0)]
	} else {
	   set slope [expr $emag / $astop]
	   }
	set xcen0 $xcen
	set ycen0 $ycen
	set xcen [expr $xcen0 + $fract * $slope * $xapp]
	set ycen [expr $ycen0 + $fract * $slope * $yapp]
	incr ntry
	if {$ntry >= 100} break
	}
   if {$ntry == 100} {
	echo No convergence in entranceMap
	}
   return [list $xcen $ycen]
   }

###########################################################################
#Similar proc to "ray" but fine-tune the entrance pupil before running the
#ray.

proc rayPupil {hndl xmm ymm xfract yfract ifil stopcheck} {
   set list [entranceMap $hndl $xmm $ymm $ifil]
   set xcen [lindex $list 0]
   set ycen [lindex $list 1]
   return [ray $hndl $xmm $ymm [expr $xfract - $xcen] [expr $yfract - $ycen] \
	$ifil $stopcheck]
   }

#######################################################################
#Inverse of entranceMap.
#I will create a new optic structure with just the surfaces in front of the
#entrance pupil and run rays from the aperture stop back to the entrance.
#I will create a fictitious "focal plane" if there is no finite object
#distance.
#This routine presumes the aperture stop is on the z axis and that
#there are no tilts.
#
#This is now the preferred way to map entrance pupil aberrations.
	
proc apertureMap {hndl ifil} {
   set surfids [surfIdsGet $hndl $ifil]

#Get the aperture stop
   set iapp -1
   loop i 0 [llength $surfids] {
	set surfid [lindex $surfids $i]
	if {[showSurf $hndl $surfid stoptype] == 2} {
	   set iapp $i
	   set apid $surfid
	   break
	   }
	}
   if {$iapp < 0} {
	error "No aperture stop found for filter $ifil"
	}

   set xapp [showSurf $hndl $apid x]
   set yapp [showSurf $hndl $apid y]
   set zapp [showSurf $hndl $apid z]
   set theta [showSurf $hndl $apid theta]
   if {$xapp != 0 || $yapp != 0 || $theta != 0} {
	echo Aperture stop is not centered/aligned with optical axis
	return
	}

   set hndl1 [opticNew]

#Copy full optic structure, then delete unwanted stuff.
   opticCopy $hndl $hndl1
   set nsurf [exprGet $hndl1.nsurf]
   for {set i 1} {$i <= $nsurf} {incr i} {
	set surfid [surfId $hndl1 1]
	opticRemove $hndl1 $surfid
	}

#Now copy active surfaces in inverse order
   setSurf $hndl1 0 z $zapp

#Can I keep refractive indices with same sign?  That would be an interesting
#trick!
   setIndex $hndl1 0 $ifil [showIndex $hndl $apid $ifil]
   set surfids [lrange $surfids 1 $iapp]
   foreach surfid $surfids {
	opticInsert $hndl1 0
	opticSurfCopy $hndl $surfid $hndl1 1
	}

#Create a fictitious "focal plane"
   set lastsurf [exprGet $hndl1.nsurf]
   set focal [opticInsert $hndl1 $lastsurf]
   setIndex $hndl1 $focal $ifil [showIndex $hndl 0 $ifil]
   setGlass $hndl1 $focal air
   set zfocal [showSurf $hndl 0 z]
   if {abs($zfocal) >= 1.e10} {
	setSurf $hndl1 $focal z -1.e4
	}


#Create a fictitious "entrance pupil".  In this case it is the location and
#size of the focal plane as viewed from the aperture stop.  It is used
#by "ray" to set the angles of rays.  For now, I will set the z distance
#to be the radius of the aperture stop. An xfract=1 ray goes out at 45 deg.
   setFocal $hndl1 $ifil emag 1.
   set zapp [showSurf $hndl1 1 z]
   set stop [showSurf $hndl1 1 outstop]
   set zpupil [expr $zapp + $stop]
   setFocal $hndl1 $ifil entrance $zpupil

#Need to shuffle refractive indexes.

   set nsurf [exprGet $hndl1.nsurf]
#   echo nsurf $nsurf
   for {set i 1} {$i < $nsurf} {incr i} {
	setIndex $hndl1 $i $ifil [showIndex $hndl1 [expr $i+1] $ifil]
	setGlass $hndl1 $i [showGlass $hndl1 [expr $i+1]]
	}

#Clean out any old "mapping" flags
   handleSet $hndl1.pupil<$ifil>->mode 0

#Now I can start mapping the exit pupil.  I want the slope and position
#of the last ray, then map the z-axis intercept to the angle
   set zintercepts ""
   set angles ""
   set a2list ""

#If I really wanted to, I could stop at the point where "ang" is greater
#than the maximum defined by "xrad" and "yrad".

#Keep track of last angle
   set lastang 0.

#Max size of field
   set angMax [fieldAngMax $hndl $ifil]

   loop i 0 10000 {
	set xfract [expr $i/20.]
	if {![ray $hndl1 0 0 $xfract 0 $ifil 0]} {
	   break
	   }
	set np [exprGet $hndl1.diagram->np]

#Index of "focal plane", last surface
	set ifocal [expr $np-1]
	set isurf [expr $np-2]
	set xfocal [exprGet $hndl1.diagram->xray<$ifocal>]
	set zfocal [exprGet $hndl1.diagram->zray<$ifocal>]
	set xsurf [exprGet $hndl1.diagram->xray<$isurf>]
	set zsurf [exprGet $hndl1.diagram->zray<$isurf>]
	set dx [expr 1.*($xfocal-$xsurf)]
	set dz [expr 1.*($zfocal-$zsurf)]
	set ang [expr 57.29578*atan2(-$dx,-$dz)]

#break if ang changes sign
	if {$ang*$lastang < 0} break
	set lastang $ang
	if {$dx == 0.} continue
	if {$dx > 0.} {
	   set ang [expr $ang + 180.]
	   }

	if {abs($ang) > $angMax} break
	set zintercept [expr $zsurf - $xsurf * ($dz/$dx)]
	lappend angles $ang

#For least squares, use square of angle in radians
	set a2 [expr pow($ang/57.29578,2)]
	lappend a2list $a2
	lappend zintercepts $zintercept
	}
   set plot 0
   if {$plot} {
	plotInit a
	pgEnv [eval min $angles] [eval max $angles] \
	   [eval min $zintercepts] [eval max $zintercepts] 0 0
	pgLine $angles $zintercepts
	pgSci 2
	loop i 0 [llength $angles] {
	   pgPoint [lindex $angles $i] [lindex $zintercepts $i] 3
	   }

#List trace of last successful ray.
	set xfract [expr ($i-2.)/1000.]
	ray $hndl1 0 0 $xfract 0 $ifil 0
	rayList $hndl1
	}

   opticDel $hndl1

#Make least squares fit
   for {set n 1} {$n <= 10} {incr n} {
	set coeffs [polyfit $a2list $zintercepts $n]
	set reslist [polycomp $a2list $zintercepts $coeffs]
	set rms [polyrms $reslist]
	if {$rms < .0005} break
	}

#Print warning if fit is out of limits:
   if {$n > 10 || $rms > .02} {
	echo Polynomial fit n = $n rms = $rms may be poor
	}
   return $coeffs
   }

#######################################################################
#In wide-angle lenses, the exit pupil is both displaced and distorted.
#Let's compute the distortion, once the displacement has been pinned down.
#What I actually return is the ratio of the actual size of the pupil to
#the paraxial predicted size.  This is then a multiplier to EMAG to get the
#actual pupil magnification.
#It is different in the tangential and sagittal directions.

proc pupilDistort {hndl ifil} {

#I assume apertureMap has been called already.

#First, set defaults
   handleSet $hndl.pupil<$ifil>->ntcoeff 0
   handleSet $hndl.pupil<$ifil>->nscoeff 0

#Get the aperture stop
   set surfids [surfIdsGet $hndl $ifil]
   set iapp -1
   loop i 0 [llength $surfids] {
	set surfid [lindex $surfids $i]
	if {[showSurf $hndl $surfid stoptype] == 2} {
	   set iapp $i
	   set apid $surfid
	   break
	   }
	}
   if {$iapp < 0} {
	error "No aperture stop found for filter $ifil"
	}

   set astop [showSurf $hndl $apid outstop]

#Loop through a bunch of positions.  Let's plot what the distortions look
#like.

   set tmags ""
   set smags ""
   set xrad [showFocal $hndl $ifil xrad]

   set angs ""
   set angrads ""
   loop i 0 100 {
	set xmm [expr $xrad*1.*$i/100.]
	set list [focaltosky $hndl $xmm 0 $ifil]
	set xang [expr [lindex $list 0]/60.]
	if {![ray $hndl $xmm 0 1 0 $ifil 0]} continue
	set x1 [exprGet $hndl.diagram->xlens<$iapp>]
	if {![ray $hndl $xmm 0 -1 0 $ifil 0]} continue
	set x2 [exprGet $hndl.diagram->xlens<$iapp>]

#Actual pupil radius
	set radius [expr abs($x2-$x1)/2.]
	if {$radius == 0} {
	   echo In pupilDistort, bizarre - radius = 0
	   continue
	   }
	set etmag [expr $astop/$radius]

#Now y direction
	if {![ray $hndl $xmm 0 0 1 $ifil 0]} continue
	set y1 [exprGet $hndl.diagram->ylens<$iapp>]

#Actual pupil radius
	set radius [expr abs($y1)]
	set esmag [expr $astop/$radius]
	lappend angs [expr pow($xang,2)]
	lappend angrads [expr pow($xang/57.29578,2)]
	lappend tmags $etmag
	lappend smags $esmag
	}
   set angmin 0
   set angmax [eval max $angs]
   set magmin [eval min $tmags $smags]
   set magmax [eval max $tmags $smags]
   if {$magmin == $magmax} {
	set magmax [expr $magmin+.5]
	set magmin [expr $magmin-.5]
	}
#   plotInit a
#   pgEnv $angmin $angmax 0. $magmax 0 0
#   pgSci 2
#   pgLine $angs $tmags
#   pgText 5 .8 Tangential
#   pgSci 5
#   pgLine $angs $smags
#   pgText 5 .5 Sagittal

#Make least squares fits
   for {set n 1} {$n <= 10} {incr n} {
	set tcoeffs [polyfit $angrads $tmags $n]
	set reslist [polycomp $angrads $tmags $tcoeffs]
	set rms [polyrms $reslist]
	if {$rms < .005} break
	}
   if {$n > 10 || $rms > .02} {
	echo Tangential: n = $n, rms = $rms fit may be poor
	}

   for {set n 1} {$n <= 10} {incr n} {
	set scoeffs [polyfit $angrads $smags $n]
	set reslist [polycomp $angrads $smags $scoeffs]
	set rms [polyrms $reslist]
	if {$rms < .005} break
	}
   if {$n > 10 || $rms > .02} {
	echo Saggital: n = $n, rms = $rms fit may be poor
	}

   return [list $tcoeffs $scoeffs]
   }

########################################################################
#Reset mode and ncoeff for pupil info
proc pupilInfoInit {hndl {ifils ""}} {
   if {$ifils == ""} {
	set ncolor [exprGet $hndl.ncolor]
	set ifils [range 1-$ncolor]
	}
   foreach ifil $ifils {
	handleSet $hndl.pupil<$ifil>->mode 0
	handleSet $hndl.pupil<$ifil>->ncoeff 0
	handleSet $hndl.pupil<$ifil>->ntcoeff 0
	handleSet $hndl.pupil<$ifil>->nscoeff 0
	}
   return
   }

########################################################################
#Compute entrance pupil info.
#This routine is optional, but if run, must be done after opticInfo.

proc pupilInfo {hndl {ifils ""}} {
   if {$ifils == ""} {
	set ncolor [exprGet $hndl.ncolor]
	set ifils [range 1-$ncolor]
	}
   foreach ifil $ifils {
	set coeffs [apertureMap $hndl $ifil]
	if {[llength $coeffs] == 0} continue

#apertureMap returns a set of polynomial coeffs for computing the z position
#of the exit pupil vs. incidence angle (squared) of a chief ray.  I will
#subtract off the nominal entrance pupil z position from the 0'th order
#coeff in order to make this a correction to the nominal z position.
#
#The following code should probably be in its own wrapped procedure, but
#for now, we will show it explicitly here.
	set znom [showFocal $hndl $ifil entrance]
	handleSet $hndl.pupil<$ifil>->mode 1
	handleSet $hndl.pupil<$ifil>->ncoeff [llength $coeffs]
	loop i 0 [llength $coeffs] {
	   set coeff [lindex $coeffs $i]
	   if {$i == 0} {set coeff [expr $coeff - $znom]}
	   handleSet $hndl.pupil<$ifil>->coeff<$i> $coeff
	   }
#Reset defaults for distortion
	handleSet $hndl.pupil<$ifil>->ntcoeff 0
	handleSet $hndl.pupil<$ifil>->nscoeff 0


#pupilDistort returns two sets of polynomial coeffs for computing the
#pupil distortion vs. incidence angle (squared, in radians) of the
#tangential and sagittal directions.
#
#The following code should probably be in its own wrapped procedure, but
#for now, we will show it explicitly here.
	set list [pupilDistort $hndl $ifil]
	set tcoeffs [lindex $list 0]

#If tangential coeffs didn't work, sagittal probably didn't either
	if {[llength $tcoeffs] == 0} continue
	handleSet $hndl.pupil<$ifil>->ntcoeff [llength $tcoeffs]
	loop i 0 [llength $tcoeffs] {
	   set coeff [lindex $tcoeffs $i]
	   handleSet $hndl.pupil<$ifil>->tcoeff<$i> $coeff
	   }

	set scoeffs [lindex $list 1]
	if {[llength $scoeffs] == 0} continue
	handleSet $hndl.pupil<$ifil>->nscoeff [llength $scoeffs]
	loop i 0 [llength $scoeffs] {
	   set coeff [lindex $scoeffs $i]
	   handleSet $hndl.pupil<$ifil>->scoeff<$i> $coeff
	   }
	}
   return
   }

#########################################################################
#List entrance pupil information

proc pupilList {hndl} {
   set ncolor [exprGet $hndl.ncolor]
   set format "%6s %8s %4s %s"
   set format2 "%6d %8.3f %4d"
   echo [format $format Filter z(nom) Mode "Pupil Center Coeffs"]
   for  {set ifil 1} {$ifil <= $ncolor} {incr ifil} {
	set znom [showFocal $hndl $ifil entrance]
	set mode [exprGet $hndl.pupil<$ifil>->mode]
	set ncoeff [exprGet $hndl.pupil<$ifil>->ncoeff]
	puts -nonewline [format $format2 $ifil $znom $mode]
	if {$mode > 0} {
	   loop j 0 $ncoeff {
		set coeff [exprGet $hndl.pupil<$ifil>->coeff<$j>]
		puts -nonewline [format " %5.2f" $coeff]
		}
	   }
	puts ""
	}

   echo [format $format Filter Emag Mode "Tangential Mag. Coeffs"]
   for  {set ifil 1} {$ifil <= $ncolor} {incr ifil} {
	set enom [showFocal $hndl $ifil emag]
	set mode [exprGet $hndl.pupil<$ifil>->mode]
	set ntcoeff [exprGet $hndl.pupil<$ifil>->ntcoeff]
	puts -nonewline [format $format2 $ifil $enom $mode]
	if {$mode > 0} {
	   loop j 0 $ntcoeff {
		set coeff [exprGet $hndl.pupil<$ifil>->tcoeff<$j>]
		puts -nonewline [format " %5.2f" $coeff]
		}
	   }
	puts ""
	}

   echo [format $format Filter Emag Mode "Sagittal Mag. Coeffs"]
   for  {set ifil 1} {$ifil <= $ncolor} {incr ifil} {
	set enom [showFocal $hndl $ifil emag]
	set mode [exprGet $hndl.pupil<$ifil>->mode]
	set nscoeff [exprGet $hndl.pupil<$ifil>->nscoeff]
	puts -nonewline [format $format2 $ifil $enom $mode]
	if {$mode > 0} {
	   loop j 0 $nscoeff {
		set coeff [exprGet $hndl.pupil<$ifil>->scoeff<$j>]
		puts -nonewline [format " %5.2f" $coeff]
		}
	   }
	puts ""
	}
   return
   }

######################################################################
#Plot a single ray on an exist opticPlot

proc rayPlot {hndl} {
   set np [exprGet $hndl.diagram->np]
   set xpts ""
   set zpts ""
   loop i 0 $np {
	lappend xpts [exprGet $hndl.diagram->xray<$i>]
	lappend zpts [exprGet $hndl.diagram->zray<$i>]
	}
   pgLine $xpts $zpts
   return
   }

#######################################################################
#Compute area of entrance pupil.
proc pupilArea {hndl xmm ymm ifil} {

#Radius of aperture stop
   set astop [apRadius $hndl $ifil]
   set emag [showFocal $hndl $ifil emag]
   set area [expr 3.141593*pow($astop*$emag,2)]
   set mode [exprGet $hndl.pupil<$ifil>->mode]
   if {$mode > 0} {
	set list [focaltosky $hndl $xmm $ymm $ifil]
	set xang [expr [lindex $list 0]/(60.*57.296)]
	set yang [expr [lindex $list 1]/(60.*57.296)]
	set ang [expr sqrt($xang*$xang + $yang*$yang)]
	set tmag 1.
	set smag 1.
	set ntcoeff [exprGet $hndl.pupil<$ifil>->ntcoeff]
	if {$ntcoeff > 0} {
	   loop i 0 $ntcoeff {
		set tmag [expr $tmag + [exprGet \
		   $hndl.pupil<$ifil>->tcoeff<$i>] * pow($ang, 2*$i)]
		}
	   }
	set nscoeff [exprGet $hndl.pupil<$ifil>->nscoeff]
	if {$nscoeff > 0} {
	   loop i 0 $nscoeff {
		set smag [expr $smag + [exprGet \
		   $hndl.pupil<$ifil>->scoeff<$i>] * pow($ang, 2*$i)]
		}
	   }
	set area [expr $area*$tmag*$smag]
	}
   return $area
   }
