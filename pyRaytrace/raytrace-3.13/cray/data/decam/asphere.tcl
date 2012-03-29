#Create an optical design from reading gladder's designs

#decamSpot
#opticGhost $hndl "3 4 5 6 7 8 9 10 11 12 13 14" 15 2 128 12
#ghostAdd (in kentools)
#tolExpand
#tol1Flags
#asphTol1
# ... 2 3 4
#tolRMS		compute rms in tolerance Table.
#flags	Least squares flags for tweaking.
#imgSummary
#tolCarlo <hndl> <niter>
#tempTol	Tolerance temperature
#adcTest <hndl> filter zenith	Compute image width for given SDSS filter, z
#######################################################################
#Define a spot pattern where the last spot has lscale = 1

proc decamSpot {} {
   clearSpot
   addSpot 0 0 .17
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8

#Let's weight things the way we report them
   addSpot 1.0 0 1.
   return
   }

#######################################################################
#Spot pattern for FOM.  Need a 4-fold pattern.

proc decamFomSpot {} {
   clearSpot

#Middle is weighted more but contributes to 4 quadrants.
   addSpot 0 0 .34
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8
   addSpot 1.0 1.
   addSpot -.2 0 .2
   addSpot -.4 0 .4
   addSpot -.6 0 .6
   addSpot -.8 0 .8
   addSpot -1.0 0 1.
   addSpot 0 .2 .2
   addSpot 0 .4 .4
   addSpot 0 .6 .6
   addSpot 0 .8 .8
   addSpot 0 1.0 1.
   addSpot 0 -.2 .2
   addSpot 0 -.4 .4
   addSpot 0 -.6 .6
   addSpot 0 -.8 .8
   addSpot 0 -1.0 1.
   return
   }

#######################################################################
#Define a spot pattern for use with tolerancing of tilts and decenters.

proc tiltSpot {} {
   clearSpot
   addSpot 0 0 .17
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8

#Let's weight things the way we report them
   addSpot 1.0 0 1.
   addSpot -.2 0 .2
   addSpot -.4 0 .4
   addSpot -.6 0 .6
   addSpot -.8 0 .8

#Let's weight things the way we report them
   addSpot -1.0 0 1.
   return
   }

#######################################################################
#For monte carlo, we want all secondary variables to be twiddled.

proc tolExpand {} {
   global tolTable

   foreach surf "2.02 2.03 2.04" {
	foreach param "x z theta" {
	   if {![info exists tolTable(2.01,$param,tol)]} continue
	   foreach type "fomlim tol type weight" {
		set tolTable($surf,$param,$type) $tolTable(2.01,$param,$type)
		}
	   set combo [list $surf $param]
	   if {[lsearch $tolTable() $combo] < 0} {
		lappend tolTable() $combo
		}
	   }
	}

#Construct list of surfaces.
   set surfs ""
   foreach combo $tolTable() {
	set surf [lindex $combo 0]
	if {[lsearch $surfs $surf] < 0} {lappend surfs $surf}
	}
   foreach surf $surfs {
	if {[info exists tolTable($surf,x,tol)]} {
	   foreach param "fomlim tol type weight" {
		set tolTable($surf,y,$param) $tolTable($surf,x,$param)
		}
	   set combo [list $surf y]
	   if {[lsearch $tolTable() $combo] < 0} {
		lappend tolTable() $combo
		}
	   }
	if {[info exists tolTable($surf,theta,tol)]} {
	   foreach param "fomlim tol type weight" {
		set tolTable($surf,phi,$param) $tolTable($surf,theta,$param)
		}
	   set tolTable($surf,phi,tol) 3.14
	   set combo [list $surf phi]
	   if {[lsearch $tolTable() $combo] < 0} {
		lappend tolTable() $combo
		}
#Increase tolerance on theta
	   set tolTable($surf,theta,tol) [expr $tolTable($surf,theta,tol) \
		* 1.414]
	   }
	}
   return
   }

#######################################################################
#Flags for use when tolerancing radii.
#Allow recentering of all lenses except front window.
#Allow refocus of primary

proc tol1Flags {hndl} {
   clearFlags
   clearSpot

   addSpot 0 0 .1
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8

#Keep scale fixed?  Maybe not - at least for tolerancing studies.
   addSpot 1. 0 3

   linkColorFlag $hndl "1 2" 2
   linkColorFlag $hndl "3 4" 4
   linkColorFlag $hndl "5 6" 6
   linkColorFlag $hndl "7 8" 8

   setSurfFlag $hndl 2.01 z 1
   setSurfFlag $hndl 2.02 z 1
   setSurfFlag $hndl 2.03 z 1
   setSurfFlag $hndl 2.04 z 1

   linkSurfFlag $hndl "3 4" z 2
   linkSurfFlag $hndl "5 6" z 3
   linkSurfFlag $hndl "7 8" z 4
   linkSurfFlag $hndl "11 12" z 5
#   linkSurfFlag $hndl "13 14" z 6

#Set other parameters so we can run fastLstsq right away.
   global _optic
   set _optic(lscale) ""
   set _optic(fomType) spot
   set _optic(niter) 1
   set _optic(flagcache) $_optic(flags)
   return
   }

#######################################################################
#Generic procedure to tolerance correctors.  Build in a little intelligence.
#Find all surfaces that are not stops.  Check curvature of a surface; if 0,
#don't tolerance ccon or x, y.
#For primary mirror, don't tolerance the shape.

#This version is used with tol1Flags.  Tolerance radii of curvature assuming
#we can recenter lenses.

#I provide increment in F.O.M. per degree of freedom as the target.

#My standard tolTarget is 1.e-4.

proc asphTol1 {hndl tolTarget} {
   global tolTable
   tolInit
   set nsurf [exprGet $hndl.nsurf]

#Tolerance the radii of all lenses.  Tolerance the refraction indices of
#glasses
   loop i 1 $nsurf {
	set surfid [surfId $hndl $i]
	set index [showIndex $hndl $surfid 1]
	if {$index == 0} continue
	set stoptype [showSurf $hndl $surfid stoptype]
	if {$stoptype != 0} continue
	set curv [showSurf $hndl $surfid curv]
	if {$curv != 0} {
	   tolTarget $surfid curv $tolTarget
	   }

#If surface is followed by glass, tolerance the refraction indices and the
#z position; we will reposition the lenses in the compensation step, so
#this does give a true measure of the tolerance on thickness.
	if {abs($index) > 1} {
	   tolTarget $surfid 1 $tolTarget
	   tolTarget $surfid z $tolTarget
	   }
	}

   foreach combo $tolTable() {
	echo $combo
	}
   tolLimit $hndl 1 0 .02
   return
   }

#######################################################################
#Flags for use when tolerancing full system.
#I refocus the primary mirror - that's all

proc tol2Flags {hndl} {
   clearFlags
   clearSpot

   addSpot 0 0 .1
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8
   addSpot 1. 0 3

   linkColorFlag $hndl "1 2" 2
   linkColorFlag $hndl "3 4" 4
   linkColorFlag $hndl "5 6" 6
   linkColorFlag $hndl "7 8" 8

#Refocus primary.
   setSurfFlag $hndl 2.01 z 1
   setSurfFlag $hndl 2.02 z 1
   setSurfFlag $hndl 2.03 z 1
   setSurfFlag $hndl 2.04 z 1

#Set other parameters so we can run fastLstsq right away.
   global _optic
   set _optic(lscale) ""
   set _optic(fomType) spot
   set _optic(niter) 1
   set _optic(flagcache) $_optic(flags)
   return
   }

#######################################################################
#Generic procedure to tolerance correctors.  Build in a little intelligence.
#Find all surfaces that are not stops.  Check curvature of a surface; if 0,
#don't tolerance ccon or x, y.  I will otherwise tolerance ccon just to
#use it as a surrogate for a surface shape parameter.
#For primary mirror, don't tolerance the shape.
#In fact, because I cannot link surfaces for tolerancing purposes yet, I
#will ignore primary here and deal with it in a separate routine completely.

proc asphTol2 {hndl tolTarg} {
   global tolTable
   tolInit
   set nsurf [exprGet $hndl.nsurf]

#Tolerance z of focal plane.
   for {set i 1} {$i <= $nsurf} {incr i} {
	set surfid [surfId $hndl $i]
	set index [showIndex $hndl $surfid 1]
	if {$index == 0} continue
	set stoptype [showSurf $hndl $surfid stoptype]
	if {$stoptype != 0} continue
	set params "z {theta phi}"
	set curv [showSurf $hndl $surfid curv]
	if {$curv != 0.} {
	   lappend params {x y} curv
	} else {
	   lappend params curv
	   }
	foreach param $params {
	   tolTarget $surfid $param $tolTarg
	   }

#Link refreaction indices
	if {abs($index) > 1} {
	   tolTarget $surfid [list "1 2 3 4 5 6 7 8"] $tolTarg
	   }
	}

#Treat decentering of the primary special.
   set params [list "x y" "theta phi"]
   foreach param $params {
	tolTarget "2.01 2.02 2.03 2.04" $param $tolTarg
	}
   foreach combo $tolTable() {
	echo $combo
	}
   tolLimit $hndl 100 0 .02
   return
   }

#######################################################################
#Flags for use when tolerancing inhomogeneities in glass.
#No refocus.

proc tol3Flags {hndl} {
   clearFlags
   clearSpot

   addSpot 0 0 .1
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8
   addSpot 1. 0 3

   linkColorFlag $hndl "1 2" 2
   linkColorFlag $hndl "3 4" 4
   linkColorFlag $hndl "5 6" 6
   linkColorFlag $hndl "7 8" 8

#Set other parameters so we can run fastLstsq right away.
   global _optic
   set _optic(lscale) ""
   set _optic(fomType) spot
   set _optic(niter) 1
   set _optic(flagcache) $_optic(flags)
   return
   }

#######################################################################
#Tolerance refraction index without refocus.  This sort of sets a limit
#on inhomogeneities.

proc asphTol3 {hndl tolTarg} {
   global tolTable
   tolInit
   set nsurf [exprGet $hndl.nsurf]

#Tolerance z of focal plane.
   for {set i 1} {$i <= $nsurf} {incr i} {
	set surfid [surfId $hndl $i]
	set index [showIndex $hndl $surfid 1]
	if {$index == 0} continue
	set stoptype [showSurf $hndl $surfid stoptype]
	if {$stoptype != 0} continue

#Link refraction indices
	if {abs($index) > 1} {
	   tolTarget $surfid [list "1 2 3 4 5 6 7 8"] $tolTarg
	   }
	}

   foreach combo $tolTable() {
	echo $combo
	}
   tolLimit $hndl 100 0 .02
   return
   }

#######################################################################
#Flags for use when tolerancing decentering of primary.  Assumption is that
#tilt primary to recollimate as best as possible.

proc tol4Flags {hndl} {
   clearFlags
   clearSpot

   addSpot 0 0 .1
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8
   addSpot 1. 0 3

   linkColorFlag $hndl "1 2" 2
   linkColorFlag $hndl "3 4" 4
   linkColorFlag $hndl "5 6" 6
   linkColorFlag $hndl "7 8" 8

   setSurfFlag $hndl 15 z 1
   linkSurfFlag $hndl "2.01 2.02 2.03 2.04" theta 2

#Set other parameters so we can run fastLstsq right away.
   global _optic
   set _optic(lscale) ""
   set _optic(fomType) spot
   set _optic(niter) 1
   set _optic(flagcache) $_optic(flags)
   return
   }

#######################################################################
#Tolerance refraction index without refocus.  This sort of sets a limit
#on inhomogeneities.

proc asphTol4 {hndl tolTarg} {
   global tolTable
   tolInit
   set nsurf [exprGet $hndl.nsurf]

   tolTarget "2.01 2.02 2.03 2.04" "x 2" $tolTarg

   foreach combo $tolTable() {
	echo $combo
	}
   tolLimit $hndl 100 0 .02
   return
   }

#######################################################################
#Solve for common scale factor
proc scaleFlags {hndl} {
   clearFlags
   linkColorFlag $hndl "1 2" 2
   linkColorFlag $hndl "3 4" 4
   linkColorFlag $hndl "5 6" 6
   linkColorFlag $hndl "7 8" 8

   linkFocalFlag $hndl "1 2 3 4 5 6 7 8" scale 10

   clearSpot
   addSpot 0 0 1 1
   addSpot 0 1 1 1
   return
   }

#######################################################################
#Compute RMS error in tolerance table

proc tolRMS {} {
   global tolTable
   set fom $tolTable(fom)
   set rms 0
   foreach combo $tolTable() {
	set surf [lindex $combo 0]
	set param [lindex $combo 1]
	set fomlim $tolTable($surf,$param,fomlim)
	set weight $tolTable($surf,$param,weight)
	set rms [expr $rms + $weight * (pow($fom+$fomlim,2)-pow($fom,2))]
	}
   set rms [expr sqrt($rms)]
   return $rms
   }

#######################################################################
#This routine runs in kentools.  Combine images from running:
#   opticGhost $hndl "3 4 5 6 7 8 9 10 11 12 13 14" 15 2 64
#If defaults are blank, I have hardwired the apr6 design.  Otherwise it
#is generic.

proc ghostAdd {{surfs ""} {filter ""} {factor ""}} {
   if {$surfs == ""} {set surfs "3 4 5 6 7 8 9 10 11 12 13 14"}
   if {$filter == ""} {set filter 2}
   foreach surf $surfs {
	set hndl($surf) [openimage ghost-$surf-$filter.fit]
	set nrow [exprGet $hndl($surf).nrow]
	set ncol [exprGet $hndl($surf).ncol]
	}

#Surfaces with standard MgF2
   set mgf2surfs "3 14"
   foreach surf $surfs {
	if {$factor == ""} {
	   if {[lsearch $mgf2surfs $surf]} {
		set reflect .016
	   } else {
#No more sol-gel
#		set reflect .008
		set reflect .016
		}
	   set ccd .15
	   set scale [expr $ccd * $reflect]
	} else {
	   set scale $factor
	   }
	*con $hndl($surf) $scale
	}
   set out [regNew -type FL32 $nrow $ncol]
   regSetWithDbl $out 0
   foreach surf $surfs {
	+file $out $hndl($surf)
	}
   imgSmooth $out
   imgSmooth $out
   imgSmooth $out
   echo Output file is out.fit
   regWriteAsFits $out out.fit
   foreach surf $surfs {
	regDel $hndl($surf)
	}
   regDel $out
   return
   }

####################################################################
#Compute maximum gradient in ghost image and compare with DES specification.

proc ghostSummary {{file out.fit}} {
   global img

   set reg [openimage $file]
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]
   set rcen [expr $nrow/2]
   set ccen [expr $ncol/2]
   set img(sky) 0
   axSet row $rcen
   axSet col $ccen
   *con $reg 100.

#Don't go to corner - that is off the focal plane!
   set rmax [expr round($rcen*1)]
   set profile [annulus $reg $rmax]

#Units of flux are percent of incident sky.
   set n 0
   chainEach hndl $profile {
	set rad($n) [exprGet $hndl.rad]
	set cnts($n) [exprGet $hndl.cnts]
	incr n
	}
   set maxgrad 0

#Scale in mm per pixel.
   set scale [expr 225./$rcen]

   set rmax 0
   loop i 1 $n {
	set i1 [expr $i-1]
	set dr [expr ($rad($i) - $rad($i1)) * $scale]
	set grad [expr abs($cnts($i) - $cnts($i1))/$dr]
	if {$grad > $maxgrad} {
	   set maxgrad $grad
	   set rmax $rad($i)
	   }
	}

#Specification is that gradient shall not exceed .05% per mm
   set mmmax [expr $rmax*$scale]
   set mmmax [format %.0f $mmmax]
   echo Maximum gradient is [format %.4f $maxgrad] percent per mm at r = \
	$rmax pixels = $mmmax mm
   echo Requirement is no bigger than .05 percent per mm.
   closeimage $reg
   return
   }

####################################################################
#Scratch routine for setting flags for least squares.
#Keep this oriented towards may11

proc flags {hndl} {
   clearFlags
   stopcomp $hndl

#Focus?
#   focusFlag $hndl
   foreach i "1 2 4" {
	setSurfFlag $hndl 2.0$i z 1
	}

#Colors
   colorLink $hndl
   setFocal $hndl 1 weight 1.5
   setFocal $hndl 2 weight 1.5
   setFocal $hndl 3 weight 1
   setFocal $hndl 4 weight 1
   setFocal $hndl 5 weight 1
   setFocal $hndl 6 weight 1
   setFocal $hndl 7 weight 2
   setFocal $hndl 8 weight 2

#Let's resolve for curvatures and adjust spacings of all elements except 1st
#and last

#Keep certain curvatures fixes.
   foreach i "3 4 5 6 8 11 12 13" {setSurfFlag $hndl $i curv 1}

   set asurfs "11"

   foreach asurf $asurfs {
	setSurfFlag $hndl $asurf a4 1
	setSurfFlag $hndl $asurf a6 1
	setSurfFlag $hndl $asurf a8 1
	setSurf $hndl $asurf a10 0.
	}


#Clear out old asphere
   set unspheres "12"
   foreach unsphere $unspheres {
	setSurf $hndl $unsphere a4 0
	setSurf $hndl $unsphere a6 0
	setSurf $hndl $unsphere a8 0
	}

   clearSpot
   decamSpot

   return
   }

####################################################################
#Resolve blanco-2605 with different flags

proc blancoFlags {hndl} {
   clearFlags
   stopcomp $hndl

#Colors
   colorLink $hndl
   setFocal $hndl 1 weight 1
   setFocal $hndl 2 weight 1
   setFocal $hndl 3 weight 1
   setFocal $hndl 4 weight 1
   setFocal $hndl 5 weight 1.5
   setFocal $hndl 6 weight 1.5
   setFocal $hndl 7 weight 2
   setFocal $hndl 8 weight 2

#Focus
   focusFlag $hndl

#Let's resolve for curvatures and adjust spacings of all elements except 1st
#and last

#Keep certain curvatures fixes.
#Keep surface 12 fixed; allow window to flex on both sides.
   foreach i "3 5 7 11 12 13" {setSurfFlag $hndl $i curv 1}

#Asphere surface(s).  11 and 12 are for blanco-2605
#Keep 11 fixed.

   set asurfs "6"

   foreach asurf $asurfs {
	setSurfFlag $hndl $asurf a4 1
	setSurfFlag $hndl $asurf a6 1
	setSurfFlag $hndl $asurf a8 1
	setSurf $hndl $asurf a10 0.
	}

#Use a conic on surface 6?
   setSurf $hndl 5 ccon 0
   setSurf $hndl 6 ccon 0
   setSurf $hndl 7 ccon 0

#Clear out old asphere
   set unspheres "5 12"
   foreach unsphere $unspheres {
	setSurf $hndl $unsphere a4 0
	setSurf $hndl $unsphere a6 0
	setSurf $hndl $unsphere a8 0
	}

   decamSpot
   return
   }

###########################################################################
#Some analysis of images

proc imgSummary {hndl} {

#Check D80 in arcsec
   set verbose [verbose]
   verbose 0
   echo Requirement for FWHM: maximum size is .4 arcsec.
   set format "%3s %10s %10s %10s %10s"
   echo [format $format XMM g r i z]
   set format "%3.0f %10.2f %10.2f %10.2f %10.2f"
   set filters "g r i z"
   set family vista
   echo Using vista filters, since des not installed yet.
   loop i 0 [llength $filters] {
	set filter [lindex $filters $i]
	set icolor1 [expr 2*$i+1]
	set vfil ${filter}prime
	set design($filter) [filterSample $hndl $icolor1 $family $vfil]
	set sum($filter) 0.
	}
   set wsum 0.
   set xrad [showFocal $hndl 1 xrad]
   foreach xfract "0 .2 .4 .6 .8 1." {
	set xmm [expr $xfract*$xrad]
	set ymm 0
	set weight $xfract
	set wsum [expr $wsum + $weight]
	foreach filter "g r i z" {
	   set ncolor [exprGet $design($filter).ncolor]
	   set icolors ""
	   for {set i 1} {$i <= $ncolor} {incr i} {
		lappend icolors $i
		}

#Convert D80 to FWHM
	   set dsec($filter) [expr [raySum $design($filter) $xmm $ymm \
		$icolors]/1.53]
	   set sum($filter) [expr $sum($filter) + \
		pow($dsec($filter),2)*$weight]
	   }
	echo [format $format $xmm $dsec(g) $dsec(r) $dsec(i) $dsec(z)]
	}

   foreach filter "g r i z" {
	opticDel $design($filter)
	set rms($filter) [expr sqrt($sum($filter)/$wsum)]
	}

   set format "%4s%10.2f %10.2f %10.2f %10.2f"
   echo
   echo [format $format RMS $rms(g) $rms(r) $rms(i) $rms(z)]

#Targets
   set targtot .27 ;# From R. Bernstein
   set targi .22 ;# From R. Bernstein

#Alignment tolerance
   set targalign [expr sqrt(pow($targtot,2) - pow($targi,2))]

   set targr [expr sqrt(pow($targtot * pow((.69+.82)/(.56+.68),.4),2) - \
	pow($targalign,2))]

   set targg [expr sqrt(pow($targtot * pow((.69+.82)/(.39+.54),.4),2) - \
	pow($targalign,2))]

   set targu [expr sqrt(pow($targtot * pow((.69+.82)/(.35+.36),.4),2) - \
	pow($targalign,2))]

   set targz $targi
   echo [format $format TARG $targg $targr $targi $targz]
   echo

#Grand RMS - useful for doing global comparisons.
#   set rms [expr sqrt($sum/$n)]
#   echo Global RMS D80 = [format %.2f $rms]

#Check D80 for u filter.
   set design(u) [filterSample $hndl 1 $family uprime]

#Refocus
   clearFlags
   set ncolor [exprGet $design(u).ncolor]
   set colors ""
   loop i 0 $ncolor {
	lappend colors [expr $i+1]
	}
   linkColorFlag $design(u) $colors 2
   setSurfFlag $design(u) 2.01 z 1
   decamSpot
   fastLstsqInit
   fastLstsq $design(u)

   set format "%3s %10s"
   echo [format $format XMM u]
   set format "%3.0f %10.2f"
   set sum(u) 0.

   foreach xfract "0 .2 .4 .6 .8 1" {
	set xmm [expr $xfract*$xrad]
	set weight $xfract
	set ymm 0
	foreach icolor "1" {
	   set dsec(u) [expr [raySum $design(u) $xmm $ymm $colors]/1.53]
	   set sum(u) [expr $sum(u) + pow($dsec(u),2)*$weight]
	   }
	echo [format $format $xmm $dsec(u)]
	}

#Grand RMS - useful for doing global comparisons.
   set rms(u) [expr sqrt($sum(u)/$wsum)]
   echo "RMS " = [format %.2f $rms(u)]
   echo TARG = [format %.2f $targu]
   opticDel $design(u)
   verbose $verbose
   return

#Check spacing for filter:
   stopcomp $hndl
   set radius [showSurf $hndl 9 outstop]
   set diam [expr 2.*$radius]
   set z8 [showSurf $hndl 8 z]
   set z9 [showSurf $hndl 9 z]
   set diff [expr $z8 - $z9]
   set target [expr $diam + 51.]
   echo Air gap for filter: Minimum is [format %.0f $target], actual is \
[format %.0f $diff]

# Filter - lens separation
   set z10 [showSurf $hndl 10 z]
   set z11 [showSurf $hndl 11 z]
   set diff [expr $z10 - $z11]
   echo Filter - lens separation: spec is 10, actual is [format %.0f $diff]

#Check window-focal plane separation
   set z14 [showSurf $hndl 14 z]
   set z15 [showSurf $hndl 15 z]
   set diff [expr $z14-$z15]
   echo Dewar Window separation: spec is 40, actual is [format %.0f $diff]

#Check shutter space
   set z12c [showSurf $hndl 12 z]
   set stop12 [showSurf $hndl 12 outstop]
   set isurf [surfIndex $hndl 12]
   set z12e [expr [zsurf $hndl $isurf $stop12 0] + $z12c]
   set z12min [expr min($z12c,$z12e)]

   set z13c [showSurf $hndl 13 z]
   set stop13 [showSurf $hndl 13 outstop]
   set isurf [surfIndex $hndl 13]
   set z13e [expr [zsurf $hndl $isurf $stop13 0] + $z13c]
   set z13max [expr max($z13c,$z13e)]
   set diff [expr $z12min - $z13max]
   echo Shutter gap: minimum is 51, actual is [format %.0f $diff]

#Total length of assembly
   set z3 [showSurf $hndl 3 z]
   set diff [expr $z3 - $z15]
   echo Total length of corrector assembly: maximum is 2792, \
	   actual is [format %.0f $diff]

#Scale factor
   set scale [showScale $hndl 1]
   echo Scale factor: spec is 17.7 - 18 arcsec/mm, \
	actual is [format %.1f $scale]

#Field of view
   set fov [expr 2.*[showFocal $hndl 1 xsize]]
   echo Field diameter: spec is 451 mm, actual is [format %.1f $fov]

#Thickness of filter
   set z9 [showSurf $hndl 9 z]
   set z10 [showSurf $hndl 10 z]
   set diff [expr $z9 - $z10]
   echo Filter thickness: spec is 15, actual is [format %.1f $diff]

#Window max/min thickness
   set z14c [showSurf $hndl 14 z]
   set stop14 [showSurf $hndl 14 outstop]
   set isurf [surfIndex $hndl 14]
   set z14e [expr [zsurf $hndl $isurf $stop14 0] + $z14c]

   set z13c [showSurf $hndl 13 z]
   set stop13 [showSurf $hndl 13 outstop]
   set isurf [surfIndex $hndl 13]
   set z13e [expr [zsurf $hndl $isurf $stop13 0] + $z13c]

   set cthick [expr $z13c - $z14c]
   set ethick [expr $z13e - $z14e]
   set minthick [expr min($cthick, $ethick)]

   echo Window min thickness: minimum is 25, actual is [format %.1f $minthick]
   verbose $verbose
   return
   }

#######################################################################
#General procedure to create 4 surfaces for a primary mirror for my 4
#filters.  This is useful when inheriting ZEMAX designs that have different
#focus for each filter.

#I specify the surfId of the primary mirror.

proc primaryDup {hndl primsurf} {

#Index of primsurf
   set index1 [surfIndex $hndl $primsurf]
   setSurfName $hndl $index1 $primsurf.01
   opticInsert $hndl $primsurf.01
   opticInsert $hndl $primsurf.01
   opticInsert $hndl $primsurf.01
   set index2 [expr $index1+1]
   set index3 [expr $index1+2]
   set index4 [expr $index1+3]

   setSurfName $hndl $index2 $primsurf.02
   setSurfName $hndl $index3 $primsurf.03
   setSurfName $hndl $index4 $primsurf.04

#Copy mirror parameters
   set vars "curv ccon x y z a4 a6 a8 instop outstop stoptype"
   foreach var $vars {
	set $var [showSurf $hndl $primsurf.01 $var]
	}
   foreach i "2 3 4" {
	foreach var $vars {
	   setSurf $hndl $primsurf.0$i $var [set $var]
	   }
	}

#Refractive indices
   foreach ifil "3 4 5 6 7 8" {
	setIndex $hndl $primsurf.01 $ifil 0
	}
   foreach ifil "3 4" {
	setIndex $hndl $primsurf.02 $ifil -1
	}
   foreach ifil "5 6" {
	setIndex $hndl $primsurf.03 $ifil -1
	}
   foreach ifil "7 8" {
	setIndex $hndl $primsurf.04 $ifil -1
	}

   set name [showName $hndl $primsurf.01]
   foreach i "2 3 4" {
	setName $hndl $primsurf.0$i $name
	}
   return
   }

#######################################################################
#Monte Carlo.
#I assume I have already run 4 tolerancings and saved them to a file:
#   design1tol.tcl	Tolerance accuracy of curvatures and indexes,
#			adjust spacings of elements
#   design2tol.tcl	Mega tolerancing.  Tolerances on curvature and index
#			are for measurement.
#   design4tol.tcl	Primary mirror displacement and tilt.  (I should
#			merge this into 2, because I adjust tilt in my
#			current scheme as a compensation variable.)
#   design3tol.tcl	Refraction index homogeneity.  No refocus here.

proc tolCarlo {design niter} {
   global tolTable
   set hndl [lensRead ${design}.len]
   set file1 ${design}tol1.tcl
   set file2 ${design}tol2.tcl
   set file3 ${design}tol3.tcl
   set file4 ${design}tol4.tcl
   global FORK

   if {[info exists FORK]} {
	set command forkLoop
   } else {
	set command loop
	}
   $command i 0 $niter {
	source $file1
	tolToOptic
#Need to randomize the random number generator!
        loop j 0 [expr $i*[llength $tolTable()]] {
           expr rand()
           }
	set hndl1 [monteCarlo $hndl]
	source $file2
	tolExpand
	tolToOptic
#Need to randomize the random number generator!
        loop j 0 [expr $i*[llength $tolTable()]] {
           expr rand()
           }
	set hndl2 [monteCarlo $hndl1]
	source $file4
	tolToOptic
#Need to randomize the random number generator!
        loop j 0 [expr $i*[llength $tolTable()]] {
           expr rand()
           }
	set hndl4 [monteCarlo $hndl2]
	source $file3
	tolToOptic
#Need to randomize the random number generator!
        loop j 0 [expr $i*[llength $tolTable()]] {
           expr rand()
           }
	set hndl3 [monteCarlo $hndl4]
	set fom [fom $hndl3]
	set file mc$i
	set fid [open $file w]
	puts $fid "$fom"
	close $fid
	opticDel $hndl1
	opticDel $hndl2
	opticDel $hndl3
	opticDel $hndl4
	}

   set fomlist ""
   loop i 0 $niter {
        set file mc$i
        set fid [open $file]
        set fom [gets $fid]
        close $fid
        exec rm $file
        lappend fomList $fom
        }
   set tolTable(monteCarlo) $fomList

#Recompute target
   set rms 0
   foreach file "$file1 $file2 $file3 $file4" {
	source $file

#Recompute target fom.  I've edited all tolerances
	tolToOptic
	tolLimit $hndl 1

#Write out result.
	arrayWrite tolTable $file
	set rms [expr $rms + pow([tolRMS],2)]
	}
   set fomtarg [expr sqrt(pow($tolTable(fom),2) + $rms)]
   set tolTable(fomtarg) $fomtarg
   set tolTable(monteCarlo) $fomList
   opticDel $hndl

   return
   }

#############################################################################
#Compute effects of temperature on design.  dt is temp change in deg. K
#Copy optics structure and return one for new temperature.

proc tempTol {hndl dt} {
   set hndl1 [opticNew]
   opticCopy $hndl $hndl1

#Need to change index of refraction and spacing
#Change in fused silica index with temperature
   set dndt 1.e-5

#Note: BK7 has dndt that depends on wavelength 2e-6 (blue) to 1e-6 (red)

#Tempco of steel.
   set ppmsteel 1.2e-5

#Tempco of aluminum
   set ppmsteel 2.3e-5
   set ppmfs .5e-6

#Reposition  lenses, change indices for corrector and focal plane.
#Primary will be refocus, so don't bother with it.

   foreach lens [list "3 4" "5 6" "7 8" "9 10" "11 12" "13 14"] {
	set surfid1 [lindex $lens 0]
	set surfid2 [lindex $lens 1]
	set z1 [showSurf $hndl1 $surfid1 z]
	set z2 [showSurf $hndl1 $surfid2 z]
	set thick [expr $z2 - $z1]

#Rescale z for back of lens using steel; rescale thickness of front using
#FS:
	set z1 [expr $z1 * (1. + $ppmsteel*$dt)]
	set z2 [expr $z1 + $thick * (1. + $ppmfs*$dt)]
	setSurf $hndl1 $surfid1 z $z1
	setSurf $hndl1 $surfid2 z $z2
	foreach icolor "1 2 3 4 5 6 7 8" {
	   set index [showIndex $hndl1 $surfid1 $icolor]
	   set index [expr $index * (1. + $dndt*$dt)]
	   setIndex $hndl1 $surfid1 $icolor $index
	   }
	}

#Focal plane
   set z1 [showSurf $hndl1 15 z]
   set z1 [expr $z1 * (1. + $ppmsteel*$dt)]
   setSurf $hndl1 15 z $z1

#Curvatures.  For this, we include the primary.
#Whoops, the primary is cervit, so we DON'T want to include.
#Cervit has CTE 0 +/- .03e-6
   foreach surf "3 4 5 6 7 8 9 10 11 12 13 14" {
	set curv [showSurf $hndl1 $surf curv]
	set curv [expr $curv/(1. + $ppmfs*$dt)]
	setSurf $hndl1 $surf curv $curv
	}
   return $hndl1
   }

#########################################################################
#List out tolerances in a form that they can be easily imported into
#scalc.

proc tolSpread {{file ""}} {
   if {$file == ""} {
	set stdout stdout
   } else {
	set stdout [open $file w]
	}
   global tolTable
   set fom $tolTable(fom)
   puts $stdout "surfs,params,tolerance,rms,weight"
   foreach combo $tolTable() {
	set surfs [lindex $combo 0]
	set surf [lindex $surfs 0]
	set params [lindex $combo 1]
	set param [lindex $params 0]
	set tol $tolTable($surfs,$params,tol)
	set weight $tolTable($surfs,$params,weight)
	set fomlim $tolTable($surfs,$params,fomlim)
	set rms [expr sqrt(pow($fom+$fomlim,2) - pow($fom,2))]
	puts $stdout [format "%s,%s,%.4g,%.4g,%d" $surf $param $tol $rms \
	   $weight]
	}
   if {$file != ""} {close $stdout}
   return
   }

##########################################################################
#Test differential refraction effects.

proc adcTest {optic filter zenith} {

#I will create 4 new colors, offset the field center by the differential
#refraction, compute the new D80 for the combined set, then reset ncolor
#to effectively delete the new colors.

   if {$filter == "g"} {
	set ifil1 1
	set ifil2 2
   } elseif {$filter == "r"} {
	set ifil1 3
	set ifil2 4
   } elseif {$filter == "i"} {
	set ifil1 5
	set ifil2 6
   } else {
	set ifil1 7
	set ifil2 8
	}
   set wave1 [showWave $optic $ifil1]
   set wave2 [showWave $optic $ifil2]
   loop i 0 4 {
	set wave [expr $wave1 + ($wave2-$wave1)*$i/4.]
	set icolor($i) [waveAdd $optic $wave $ifil1]
	foreach id "3 5 7 9 11 13" {
	   setIndex $optic $id $icolor($i) [expr -1.*[glass SIO2 $wave]]
	   }
	setFocal $optic $icolor($i) xoff 0
	}
   set xmm 0.
   set ymm 225.5
#   spotPlot $optic $xmm $ymm "$icolor(0) $icolor(1) $icolor(2) $icolor(3)"
   set old [raySum $optic $xmm $ymm \
	"$icolor(0) $icolor(1) $icolor(2) $icolor(3)"]

   loop i 0 4 {
	set wave [expr $wave1 + ($wave2-$wave1)*$i/4.]
	set offset [diffRefract $wave1 $wave $zenith]
	echo Offset $offset
	setFocal $optic $icolor($i) xoff [expr $offset/60.]
	}

   set new [raySum $optic $xmm $ymm \
	"$icolor(0) $icolor(1) $icolor(2) $icolor(3)"]
#   spotPlot $optic $xmm $ymm "$icolor(0) $icolor(1) $icolor(2) $icolor(3)"
   set diff [expr pow($new,2) - pow($old,2)]
   if {$diff < 0.} {set diff 0} else {set diff [expr sqrt($diff)]}
   loop i 0 4 {
	setIndex $optic 0 $icolor($i) 0
	}
   colorcount $optic
   return [format %.2f $diff]
   }

########################################################################
proc adcPlot {optic} {
   plotInit a
   pgEnv 0 60 0 1 0 0
   pgLabel "Zenith Angle" "Incremental D80 (arcsec)" "Differential Refraction"
   loop ifil 0 4 {
 	set filter [lindex "g r i z" $ifil]
	pgSci [lindex "3 2 6 1" $ifil]
	set xlist ""
	set ylist ""
	foreach zenith "0 15. 30. 45. 60." {
	   if {$zenith == 0.} {
		set offset 0.
	   } else {
		set offset [adcTest $optic $filter $zenith]
		}
	   pgPoint $zenith $offset 4
	   lappend xlist $zenith
	   lappend ylist $offset
	   update
	   }
	pgLine $xlist $ylist
	update
	}
   return
   }

#######################################################################
#Insert adc into may11.

proc adcInsert {hndl} {
   set z1 -10250.
   set s1 8
   set thick 25.
   surfInsert $hndl [expr $s1] [expr $z1] LLF6
   surfInsert $hndl [expr $s1+1] [expr $z1-$thick] UBK7
   surfInsert $hndl [expr $s1+2] [expr $z1-2.*$thick] air
   surfInsert $hndl [expr $s1+3] [expr $z1-2.5*$thick] LLF6
   surfInsert $hndl [expr $s1+4] [expr $z1-3.5*$thick] UBK7
   surfInsert $hndl [expr $s1+5] [expr $z1-4.5*$thick] air
   return
   }

####################################################################
#Colors
proc colorLink {hndl} {
   linkColorFlag $hndl "1 2" 2
   linkColorFlag $hndl "3 4" 4
   linkColorFlag $hndl "5 6" 6
   linkColorFlag $hndl "7 8" 8
   return
   }

######################################################################
#Focus flags
proc focusFlag {hndl} {
   foreach i "1 2 3 4" {
	setSurfFlag $hndl 2.0$i z 1
	}
   return
   }

######################################################################
#Add uv filter to end of an optics structure.

proc uvAdd {hndl} {
   set icolor0 [exprGet $hndl.ncolor]
   set wave1 .33
   set wave2 .39
   set icolor1 [waveAdd $hndl $wave1 $icolor0]
   set icolor2 [waveAdd $hndl $wave2 $icolor0]
   set i [expr $icolor0/2]
   set oldsurf 2.0$i
echo oldsurf $oldsurf
   set inssurf 3.0$i
echo inserted surface $inssurf
   set ii [expr $icolor2/2]
   set newsurf 2.0$ii
echo newsurf $newsurf

   set glass [showGlass $hndl $oldsurf]
echo glass is $glass
   set zpos [showSurf $hndl $oldsurf z]
   surfInsert $hndl $oldsurf $zpos $glass

#Change name of new surface
   set index [surfIndex $hndl $inssurf]
   setSurfName $hndl $index $newsurf


#Most parameters are unset for new surface
   foreach var "curv ccon stoptype" {
	setSurf $hndl $newsurf $var [showSurf $hndl $oldsurf $var]
	}
   setName $hndl $newsurf PRIMARY

   for {set i 1} {$i <= $icolor0} {incr i} {
	setIndex $hndl $newsurf $i 0
	}
   setIndex $hndl $newsurf $icolor1 [showIndex $hndl $oldsurf $icolor1]
   setIndex $hndl $newsurf $icolor2 [showIndex $hndl $oldsurf $icolor2]
   setIndex $hndl $oldsurf $icolor1 0
   setIndex $hndl $oldsurf $icolor2 0

   clearFlags
   linkColorFlag $hndl "$icolor1 $icolor2" 2
   setSurfFlag $hndl $newsurf z 1
   clearSpot
   addSpot 0 0 1
   addSpot .5 0 2
   addSpot 1 0 3
   fastLstsqInit
   fastLstsq $hndl
   return
   }

#########################################################################
#Evaluate stellar halo ghosts.
#Use ghostDup to create a new optic structure, run edge-rays for on-axis
#spot position.
#Use CCD and reflection efficiencies as for ghostAdd

proc haloGhost {optic {ccd .15} {reflect .016}} {
   set surfs [surfIdsGet $optic 1]
   set flag 0
   loop i 0 [llength $surfs] {
	set surf [lindex $surfs $i]
	set name [showName $optic $surf]
	if {$name == "FOCAL"} {
	   set flag 1
	   break
	   }
	}
   if {$flag == 0} {
	echo No FOCAL plane found
	return
	}
   set focal [lindex $surfs $i]
   set window [lindex $surfs [expr $i-1]]
   set ghost [ghostDup $optic $window $focal 1]
   ray $ghost 0 0 1 0 1 0
   set np [exprGet $ghost.diagram->np]
   set xmax [exprGet $ghost.diagram->xray<[expr $np-1]>]
   opticDel $ghost

#Convert to arcsec
   set scale [showScale $optic 1]
   set rad [expr $xmax*$scale]
   set magdiff [expr -2.5*log10($ccd*$reflect/($rad*$rad*3.1416))]
   set magdiff [format %.2f $magdiff]
   return $magdiff
   }

##########################################################################
#Compute sample wavelengths for a filter given the min, max limits.
#I assume that these are hard edges of filters.  I divide filter into
#n equal bands and compute the center of each band.

#I will weight the edges by .5.
proc filterWaves {wmin wmax {n 6}} {
   set dw [expr ($wmax-$wmin)/(1.*($n-1))]
   loop i 0 $n {
	set w [expr $wmin + $i*$dw]
	if {$i == 0 || $i == $n-1} {
	   set weight .5
	} else {
	   set weight 1.
	   }
	echo Wavelength [format %.2f $w] weight $weight
	}
   return
   }

########################################################################
#Insert reflection coeffs into design

proc reflect {optic} {
   foreach surf "3 4 5 6 7 8 9 10 11 12 13 14" {
	setSurf $optic $surf reflect .01
	}
   foreach surf "2.01 2.02 2.03 2.04" {
	setSurf $optic $surf reflect 1.
	}
   setSurf $optic 15 reflect .15
   return
   }

##########################################################################
#Basic set of flags for doing figure of merit analysis.  Include refocus.

proc fomFlags {hndl} {
   decamSpot
   clearFlags
   colorLink $hndl

#Adjust focus?
   set focus 1
   if {$focus} {
	foreach i "1 2 3 4" {
	   setSurfFlag $hndl 2.0$i z 1
	   }
	}
   return
   }

