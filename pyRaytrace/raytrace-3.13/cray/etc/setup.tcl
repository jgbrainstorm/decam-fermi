#Initialize sdss telescope from internal design
#basicInit	Simple lens
#rcInit
#cassInit

#Manipulate wavelength and indices
#waveAdd
#waveCopy
#waveSwitch

#Get all surfaces for a given color
#surfGet

#Guess pixel scale and set for all filters
#scaleGuess
#scaleEdge

#####################################################
proc initset {} {
   set hndl [opticNew]
   setup $hndl
   colorcount $hndl
   stopcomp $hndl
   opticinc $hndl 1
   return $hndl
   }

##########################################################
#Create a basic layout with one lens and a focal plane.
#I will have 4 wavelengths.

proc basicInit {diam fratio} {

#Now initialize surfaces.

#Aperture stop radius
   set radius [expr $diam/2.]

#Lens thickness
   set thick [expr $diam/5.]

#Half-width of focal plane in mm

   set xrad $radius

#Focal length
   set focal [expr $diam*$fratio]

#Glass type
   set glass BK7

#Scale in arcsec/mm
   set scale [expr 206265./$focal]

#Curvature of back of lens
   set curv [expr 1./($focal*([glass $glass .4] -1.))]

#Location of focal plane
   set zfocal [expr $focal + $thick/2.]

#Wavelengths
   set waves [list .39 .55 .68 .9]

#Create internal data structures.
   set optic [opticNew]
   setDesign $optic Basic

   for {set i 1} {$i <= [llength $waves]} {incr i} {
	set wave [lindex $waves [expr $i-1]]
	setFocal $optic $i xpos 0
	setFocal $optic $i ypos 0
	setFocal $optic $i xrad $xrad
	setFocal $optic $i yrad $xrad
	setFocal $optic $i scale $scale
	setFocal $optic $i wave $wave
	setFocal $optic $i dist 0
	setFocal $optic $i rot 0
	setIndex $optic 0 $i 1
	}
   colorcount $optic

#Now surfaces
#Don't set optic.nsurf - setparam automatically creates new surfaces as
#needed.

   setSurf $optic 0 z -1.e14
   setGlass $optic 0 air

#lensInsert and surfInsert need to find an aperture stop.  I will insert
#a fake surface here just to do that, then delete once the lens is in place.
   opticInsert $optic 0
   setSurf $optic 1 stoptype 2
   setSurf $optic 1 outstop $radius
   for {set i 1} {$i <= [llength $waves]} {incr i} {
	setIndex $optic 1 $i 1
	}

#Lens
   lensInsert $optic 1 0 $thick $glass

   opticRemove $optic 1

   setSurf $optic 1 stoptype 2
   setSurf $optic 1 outstop $radius
   setName $optic 1 C1
   setSurf $optic 1 curv [expr $curv/2.]
   setSurf $optic 2 curv [expr -$curv/2.]
   setName $optic 2 C1

#Focal plane.
   surfInsert $optic 2 $zfocal air
   setName $optic 3 FOCAL

   rayPattern $optic 6 1

#Compute focal length and exit pupil
   for {set i 1} {$i <= [llength $waves]} {incr i} {
	opticInfo $optic $i
	}

   stopComp $optic
   opticinc $optic 1
   return $optic
   }

##########################################################
#Create a new RC design
proc rcInit {} {
   echo You will be prompted for the following:
   echo "   " primary mirror diameter
   echo "   " primary mirror focal ratio
   echo "   " final focal ratio
#   echo "   " back focal distance as a fraction of primary mirror diameter
   echo "   " back focal distance
   echo
   while {1} {
	puts stdout "Primary mirror diameter: " nonewline
	flush stdout
	set D [gets stdin]
	if {![catch {format %f $D}]} break
	echo Try again ...
	}
   while {1} {
	puts stdout "Primary mirror focal ratio: " nonewline
	flush stdout
	set f1 [gets stdin]
	if {![catch {format %f $f1}]} break
	echo Try again ...
	}
   while {1} {
	puts stdout "Final focal ratio: " nonewline
	flush stdout
	set f2 [gets stdin]
	if {![catch {format %f $f2}]} break
	echo Try again ...
	}
   while {1} {
	puts stdout "Back-focal distance: " nonewline
	flush stdout
	set back [gets stdin]
	if {![catch {format %f $back}]} break
	echo Try again ...
	}

#Derived parameters
#Primary mirror focal length
   set f1 [expr $f1]
   set F1 [expr $f1*$D]
#Primary mirror radius of curvature
   set R1 [expr -2.*$F1]
   set curv1 [expr 1./$R1]
#Magnification
   set M [expr $f2/$f1]
#Back focal distance ratio (in terms of primary focal length)
   set beta [expr $back/$F1]
#Secondary mirror focal length
   set F2 [expr $F1*$M*(1.+$beta)/(pow($M,2)-1.)]
#Radius of curvature of secondary
   set R2 [expr -2.*$F2]
   set curv2 [expr 1./$R2]
#Delta
   set delta [expr ($beta+1.)/($M+1.)]
#Distance of secondary from focal point of primary
   set d [expr $delta*$F1]
#Mirror separation
   set sep [expr $F1-$d]
#Final focal length
   set F3 [expr $D*$f2]
#Primary conic constant
   set ck1 [expr -1. -2. * (1. + $beta) / \
        (pow($M,2) * ($M-$beta))]
#Secondary conic constant
   set ck2 [expr -pow(($M + 1.) / ($M - 1.),2) \
        + pow($M,3) * (1. + $M) * ($ck1 + 1.) / ((1. + $beta) * \
	pow($M - 1.,3))]
#Focal plane curvature
   set fmed [expr (2. / ($M * $R1)) * (((pow($M,2) - 2.) * ($M - $beta) \
        + $M * ($M + 1.)) / ($M *(1. + $beta)) - $M * pow($M - $beta,2) \
	* ($ck1 + 1.) / (2. * pow(1. + $beta,2)))]
   set curv3 [expr 1.*$fmed]
#Petzval curvature
   set fpetz [expr (2. / $R1) * (($M * ($M - $beta) - ($M + 1.)) / \
        ($M * (1. + $beta)))]

#Now create internal data structures.
   set optic [opticNew]
   handleSet $optic.tel->diam $D
   handleSet $optic.tel->fr1 $f1
   handleSet $optic.tel->fr2 $f2
   handleSet $optic.tel->back $back
   handleSet $optic.tel->fl1 $F1
   handleSet $optic.tel->ck1 $ck1
   handleSet $optic.tel->rad1 $R1
   handleSet $optic.tel->finner .25
   handleSet $optic.tel->mag $M
   handleSet $optic.tel->beta $beta
   handleSet $optic.tel->fl2 $F2
   handleSet $optic.tel->ck2 $ck2
   handleSet $optic.tel->rad2 $R2
   handleSet $optic.tel->delta $delta
   handleSet $optic.tel->D $d
   handleSet $optic.tel->sep $sep
   handleSet $optic.tel->f3 $F3
   handleSet $optic.tel->fmed $fmed
   handleSet $optic.tel->fpetz $fpetz

#Now initialize surfaces.
#Size of focal plane: Use 40 arcminute diameter
   set size [expr $F3*(20./60./57.3)]
#Scale in arcsec/mm
   set scale [expr 206265./$F3]

   setFocal $optic 1 xpos 0
   setFocal $optic 1 ypos 0
   setFocal $optic 1 xrad $size
   setFocal $optic 1 yrad $size
   setFocal $optic 1 scale $scale
   setFocal $optic 1 wave .58
   setFocal $optic 1 dist 0
   setFocal $optic 1 rot 0

#Now surfaces
#Don't set optic.nsurf - setparam automatically creates new surfaces as
#needed.
   setDesign $optic "Ritchey Cretien"
   handleSet $optic.ncolor 1

   setSurf $optic 0 z -1.e14
   setIndex $optic 0 1 1

   setSurf $optic 1 z 0
   setSurf $optic 1 curv $curv1
   setSurf $optic 1 ccon $ck1
   setSurf $optic 1 stoptype 2
   setSurf $optic 1 outstop [format %.1f [expr $D/2.]]
   setIndex $optic 1 1 -1

   setSurf $optic 2 z [expr -1.*$sep]
   setSurf $optic 2 curv $curv2
   setSurf $optic 2 ccon $ck2
   setIndex $optic 2 1 1

   setSurf $optic 3 z [expr $back]
   setSurf $optic 3 a2 [expr .5*$curv3]
   setIndex $optic 3 1 1

   rayPattern $optic 6 1
   colorcount $optic

#Compute focal length and exit pupil
   opticInfo $optic 1
   stopComp $optic
   opticinc $optic 1
   return $optic
   }

##########################################################
#Create a new Cassegrain design
proc cassInit {} {
   echo You will be prompted for the following:
   echo "   " primary mirror diameter
   echo "   " primary mirror focal ratio
   echo "   " final focal ratio
#   echo "   " back focal distance as a fraction of primary mirror diameter
   echo "   " back focal distance
   echo
   while {1} {
	puts stdout "Primary mirror diameter: " nonewline
	flush stdout
	set D [gets stdin]
	if {![catch {format %f $D}]} break
	echo Try again ...
	}
   while {1} {
	puts stdout "Primary mirror focal ratio: " nonewline
	flush stdout
	set f1 [gets stdin]
	if {![catch {format %f $f1}]} break
	echo Try again ...
	}
   while {1} {
	puts stdout "Final focal ratio: " nonewline
	flush stdout
	set f2 [gets stdin]
	if {![catch {format %f $f2}]} break
	echo Try again ...
	}
   while {1} {
	puts stdout "Back-focal distance: " nonewline
	flush stdout
	set back [gets stdin]
	if {![catch {format %f $back}]} break
	echo Try again ...
	}

#Derived parameters
#Primary mirror focal length
   set f1 [expr $f1]
   set F1 [expr $f1*$D]
#Primary mirror radius of curvature
   set R1 [expr -2.*$F1]
   set curv1 [expr 1./$R1]
#Magnification
   set M [expr $f2/$f1]
#Back focal distance ratio (in terms of primary focal length)
   set beta [expr $back/$F1]
#Secondary mirror focal length
   set F2 [expr $F1*$M*(1.+$beta)/(pow($M,2)-1.)]
#Radius of curvature of secondary
   set R2 [expr -2.*$F2]
   set curv2 [expr 1./$R2]
#Delta
   set delta [expr ($beta+1.)/($M+1.)]
#Distance of secondary from focal point of primary
   set d [expr $delta*$F1]
#Mirror separation
   set sep [expr $F1-$d]
#Final focal length
   set F3 [expr $D*$f2]
#Primary conic constant
   set ck1 -1.
#Secondary conic constant
   set ck2 [expr -pow(($M + 1.) / ($M - 1.),2) \
        + pow($M,3) * (1. + $M) * ($ck1 + 1.) / ((1. + $beta) * \
	pow($M - 1.,3))]
#Focal plane curvature
   set fmed [expr (2. / ($M * $R1)) * (((pow($M,2) - 2.) * ($M - $beta) \
        + $M * ($M + 1.)) / ($M *(1. + $beta)) - $M * pow($M - $beta,2) \
	* ($ck1 + 1.) / (2. * pow(1. + $beta,2)))]
   set curv3 [expr 1.*$fmed]
#Petzval curvature
   set fpetz [expr (2. / $R1) * (($M * ($M - $beta) - ($M + 1.)) / \
        ($M * (1. + $beta)))]

#Now create internal data structures.
   set optic [opticNew]
   handleSet $optic.tel->diam $D
   handleSet $optic.tel->fr1 $f1
   handleSet $optic.tel->fr2 $f2
   handleSet $optic.tel->back $back
   handleSet $optic.tel->fl1 $F1
   handleSet $optic.tel->ck1 $ck1
   handleSet $optic.tel->rad1 $R1
   handleSet $optic.tel->finner .25
   handleSet $optic.tel->mag $M
   handleSet $optic.tel->beta $beta
   handleSet $optic.tel->fl2 $F2
   handleSet $optic.tel->ck2 $ck2
   handleSet $optic.tel->rad2 $R2
   handleSet $optic.tel->delta $delta
   handleSet $optic.tel->D $d
   handleSet $optic.tel->sep $sep
   handleSet $optic.tel->f3 $F3
   handleSet $optic.tel->fmed $fmed
   handleSet $optic.tel->fpetz $fpetz

#Now initialize surfaces.
#Size of focal plane: Use 40 arcminute diameter
   set size [expr $F3*(20./60./57.3)]
#Scale in arcsec/mm
   set scale [expr 206265./$F3]

   setFocal $optic 1 xpos 0
   setFocal $optic 1 ypos 0
   setFocal $optic 1 xrad $size
   setFocal $optic 1 yrad $size
   setFocal $optic 1 scale $scale
   setFocal $optic 1 wave .58
   setFocal $optic 1 dist 0
   setFocal $optic 1 rot 0

#Now surfaces
#Don't set optic.nsurf - setparam automatically creates new surfaces as
#needed.
   setDesign $optic "Optical Design"
   handleSet $optic.ncolor 1

   setSurf $optic 0 z -1.e14
   setIndex $optic 0 1 1

   setSurf $optic 1 z 0
   setSurf $optic 1 curv $curv1
   setSurf $optic 1 ccon $ck1
   setIndex $optic 1 1 -1
   setSurf $optic 1 stoptype 2

   setSurf $optic 2 z [expr -1.*$sep]
   setSurf $optic 2 curv $curv2
   setSurf $optic 2 ccon $ck2
   setIndex $optic 2 1 1

   setSurf $optic 3 z [expr $back]
   setSurf $optic 3 a2 [expr .5*$curv3]
   setIndex $optic 3 1 1

   rayPattern $optic 12 1
   colorcount $optic

#Compute focal length and exit pupil
   opticInfo $optic 1
   stopComp $optic
   opticinc $optic 1
   return $optic
   }

###############################################################
#Add a new wavelength.  Copy surfaces from a trace of the specified filter.

proc waveAdd {optic wave incolor} {
   set ncolor [exprGet $optic.ncolor]
   set icolor [expr $ncolor+1]

#Copy everything except ... drat, don't know glass types.
#First, focal plane stuff
   foreach elem "xoff yoff xsize ysize scale dist rot weight map fl exit" {
	setFocal $optic $icolor $elem [showFocal $optic $incolor $elem]]
	}
   setFocal $optic $icolor wave $wave

#I might want to run opticInfo at the end, but won't

   setIndex $optic 0 $icolor 1

#Clear out any old data
   set nsurf [exprGet $optic.nsurf]
   for {set isurf 1} {$isurf <= $nsurf} {incr isurf} {
	set surf [surfId $optic $isurf]
	setIndex $optic $surf $icolor 0
	}

#Now copy surfaces.  If a glass type is supplied, use it to compute new
#refraction index.  Otherwise just copy previous color.
   set surfids [surfGet $optic $incolor]
   foreach surfid $surfids {
	set oldindex [showIndex $optic $surfid $incolor]
	set type [showGlass $optic $surfid]

#My treatment of mirrors needs improvement.
	if {$type != ""} {
	   setIndex $optic $surfid $icolor [glass $type $wave]
	} else {
	   setIndex $optic $surfid $icolor $oldindex
	   }
	}
   colorcount $optic
#   opticinc $optic 1
#   stopComp $optic
   return $icolor
   }

###############################################################
#Delete a specific filter; shift things left
proc waveDel {hndl ifil} {
   set ncolor [exprGet $hndl.ncolor]
   if {$ifil > $ncolor || $ifil < 1} return
   loop i $ifil $ncolor {
	waveCopy $hndl [expr $i+1] $i
	}
   handleSet $hndl.ncolor [expr $ncolor-1]
   return
   }

###############################################################
#Copy a configuration from one filter to another.
#Copy from icolor1 to icolor2.  Both should exist

proc waveCopy {optic icolor1 icolor2} {
   set ncolor [exprGet $optic.ncolor]
   if {$icolor1 < 1 || $icolor1 > $ncolor} {
	error "Bad first color $icolor1"
	}
   if {$icolor2 < 1 || $icolor2 > $ncolor} {
	error "Bad second color $icolor2"
	}
   if {$icolor1 == $icolor2} {
	echo Identical colors - no copying done.
	return
	}

#Copy everything except ... drat, don't know glass types.
#First, focal plane stuff
   foreach elem "xoff yoff xsize ysize scale dist rot weight map fl exit \
	wave" {
	setFocal $optic $icolor2 $elem [showFocal $optic $icolor1 $elem]]
	}

#I might want to run opticInfo at the end, but won't

   setIndex $optic 0 $icolor2 1

#Clear out any old data
   set nsurf [exprGet $optic.nsurf]
   for {set isurf 1} {$isurf <= $nsurf} {incr isurf} {
	set surf [surfId $optic $isurf]
	setIndex $optic $surf $icolor2 0
	}

#Now copy surfaces.  If a glass type is supplied, use it to compute new
#refraction index.  Otherwise just copy previous color.
   set surfids [surfGet $optic $icolor1]

#Copy just the refraction index
   foreach surfid $surfids {
	set oldindex [showIndex $optic $surfid $icolor1]

#My treatment of mirrors needs improvement.
	setIndex $optic $surfid $icolor2 $oldindex
	}
   return
   }

###############################################################
#Switch the wavelength of an existing color.  Just reset refraction indices
#and update focal plane info.  

proc waveSwitch {optic icolor wave} {

#Only focal plane item to switch is wavelength
   setFocal $optic $icolor wave $wave

#I might want to run opticInfo at the end, but won't

#Reset index of surfaces.  If a glass type is supplied, use it to compute new
#refraction index.  Otherwise just copy previous color.
   set surfids [surfGet $optic $icolor]
   foreach surfid $surfids {
	set oldindex [showIndex $optic $surfid $icolor]
	set type [showGlass $optic $surfid]
	if {$type != ""} {
	   setIndex $optic $surfid $icolor [glass $type $wave]
	} else {
	   setIndex $optic $surfid $icolor $oldindex
	   }
	}
   return
   }

#######################################################################
#Set refractive indices for all filters of a given surface.  I do this
#in case I am duplicating surfaces that have only certain filters set.
#I use the stored glass name.
#This routine is similar to surfGlassSwitch but sets all indices regardless
#of initial value.

proc setIndexAll {optic surfid} {
   set ncolor [exprGet $optic.ncolor]
   set glass [showGlass $optic $surfid]
   for {set i 1} {$i <= $ncolor} {incr i} {
        set wave [showFocal $optic $i wave]
        set index [glass $glass $wave]
        setIndex $optic $surfid $i $index
        }
   return
   }

#######################################################################
#Keep only specific wavelengths for a given surface.
#
proc keepIndex {optic surfid args} {
   set ncolor [exprGet $optic.ncolor]
   for {set i 1} {$i <= $ncolor} {incr i} {
	if {[lsearch $args $i] >= 0} continue
	setIndex $optic $surfid $i 0
	}
   return
   }

#######################################################################
#Thin wrapper for C routine opticinc.
proc opticInc {optic {step ""}} {
   if {$step == ""} {
	echo Using default step of 1 arcsec
	set step 1
	}
   opticinc $optic $step
   return
   }

#######################################################################
#Thin wrapper for C routine stopcomp
proc stopComp {optic} {
   stopcomp $optic
   return
   }

#######################################################################
#Thin wrapper for C routine opticinfo
proc opticInfo {optic {filters ""}} {
   if {$filters == ""} {
	set ncolor [exprGet $optic.ncolor]
	loop i 0 $ncolor {
	   lappend filters [expr $i+1]
	   }
	}
   foreach filter $filters {
	opticinfo $optic $filter
	}
   return
   }

#######################################################################
#Get active surfaces for a specified filter.  This is a routine that I
#use often enough that it helps to have a helper routine

proc surfGet {optic icolor} {
   ray $optic 0 0 0 0 $icolor 0
   set np [exprGet $optic.diagram->np]
   set surfs ""
   loop i 0 $np {
	set isurf [exprGet $optic.diagram->indx<$i>]
	set surfid [surfId $optic $isurf]
	lappend surfs $surfid
	}
   return $surfs
   }

#######################################################################
#Reset scale factors.  Run a ray that is just off-axis and see where
#it lands.  Reset scale factor accordingly.

proc scaleGuess {hndl {ifils ""}} {

   if {$ifils == ""} {
	set ifils [range 1-[exprGet $hndl.ncolor]]
	}

#Compute the average scale factor.
   set sum 0.
   set n 0
   foreach ifil $ifils {

#Default scale
	set scale [showScale $hndl $ifil]

#Angle in arcmin - use 1 arcsec.
	set list [skytofocal $hndl 0 0 $ifil]
	set xref [lindex $list 0]
	set yref [lindex $list 1]

	set ang [expr 1./60.]
	set list [skytofocal $hndl $ang 0 $ifil]
	set xoff [lindex $list 0]
	set yoff [lindex $list 1]

#xoff, yoff is the predicted location in the focal plane for a ray 1 arcsec off
#axis.  Now run a ray and see where it really lands.
	if {![ray $hndl $xref $yref 0 0 $ifil 0]} continue
	set np [exprGet $hndl.diagram->np]
	set np1 [expr $np-1]
	set xfref [exprGet $hndl.diagram->xray<$np1>]
	set yfref [exprGet $hndl.diagram->yray<$np1>]

	if {![ray $hndl $xoff $yoff 0 0 $ifil 0]} continue
	set np [exprGet $hndl.diagram->np]
	set np1 [expr $np-1]
	set xf [exprGet $hndl.diagram->xray<$np1>]
	set yf [exprGet $hndl.diagram->yray<$np1>]
	set r [expr sqrt(pow($xf-$xfref,2) + pow($yf-$yfref,2))]
	if {$r == 0} continue

#Arcsec/mm
	set scale [expr 60.*$ang/$r]
	set sum [expr $sum + $scale]
	incr n
	}
   set scale [expr $sum/(1.*$n)]
   foreach ifil $ifils {
	setFocal $hndl $ifil scale $scale
	}
   return $scale
   }

#######################################################################
#Reset scale factors.  Set scale so it predicts the correct edge of the
#field.

proc scaleEdge {hndl {ifils ""}} {

   if {$ifils == ""} {
	set ifils [range 1-[exprGet $hndl.ncolor]]
	}

   if {$ifils == ""} return

#Compute the average scale factor.
   set sum 0.
   set n 0
   foreach ifil $ifils {

#Default scale
	set scale [showScale $hndl $ifil]

#I will do x and y directions independently.  This will handle circular
#and rectangular focal planes at the same time.
	set xrad [showFocal $hndl $ifil xrad]
	set yrad [showFocal $hndl $ifil yrad]

	if {![ray $hndl $xrad 0 0 0 $ifil 0]} continue
	set np [exprGet $hndl.diagram->np]
	set np1 [expr $np-1]
	set xf1 [exprGet $hndl.diagram->xray<$np1>]
	set yf1 [exprGet $hndl.diagram->yray<$np1>]

	if {![ray $hndl [expr -1.*$xrad] 0 0 0 $ifil 0]} continue
	set np [exprGet $hndl.diagram->np]
	set np1 [expr $np-1]
	set xf2 [exprGet $hndl.diagram->xray<$np1>]
	set yf2 [exprGet $hndl.diagram->yray<$np1>]

#Linear separation on focal plane.  In order to test for correct sign of
#scale, I will ignore the differences in y.
	set dx [expr $xf1-$xf2]

#What the difference should have been
	set dxpred [expr 2.*$xrad]
	
#Correction factor
	set factx [expr $dxpred/$dx]

#Repeat in y direction
	if {![ray $hndl 0 $yrad 0 0 $ifil 0]} continue
	set np [exprGet $hndl.diagram->np]
	set np1 [expr $np-1]
	set xf1 [exprGet $hndl.diagram->xray<$np1>]
	set yf1 [exprGet $hndl.diagram->yray<$np1>]

	if {![ray $hndl 0 [expr -1.*$yrad] 0 0 $ifil 0]} continue
	set np [exprGet $hndl.diagram->np]
	set np1 [expr $np-1]
	set xf2 [exprGet $hndl.diagram->xray<$np1>]
	set yf2 [exprGet $hndl.diagram->yray<$np1>]

#Linear separation on focal plane.  In order to test for correct sign of
#scale, I will ignore the differences in x.
	set dy [expr $yf1-$yf2]

#What the difference should have been
	set dypred [expr 2.*$yrad]
	
#Correction factor
	set facty [expr $dypred/$dy]
	set scale [expr $scale*($factx + $facty)/2.]
	set sum [expr $sum + $scale]
	incr n
	}
   if {$n == 0} {
	echo n=0, failure!
	return
	}
   set scale [expr $sum/(1.*$n)]
   foreach ifil $ifils {
	setFocal $hndl $ifil scale $scale
	}
   return $scale
   }

#######################################################################
#Reset fiels size.

proc fieldSize {hndl rad {ifils ""}} {

   if {$ifils == ""} {
	set ifils [range 1-[exprGet $hndl.ncolor]]
	}

   foreach ifil $ifils {
	set xrad [showFocal $hndl $ifil xrad]
	set sign 1.
	if {$xrad < 0} {set sign -1.}
	setFocal $hndl $ifil xrad [expr $rad*$sign]
	set yrad [showFocal $hndl $ifil yrad]
	set sign 1.
	if {$yrad < 0} {set sign -1.}
	setFocal $hndl $ifil yrad [expr $rad*$sign]
	}
   }

###############################################################
#Get the maximum field angle for a particular filter based on the mapping
#model, not actual raytraces.  Don't assume we are on-axis.

proc fieldAngMax {hndl ifil} {
   set xrad [showFocal $hndl $ifil xrad]
   set yrad [showFocal $hndl $ifil yrad]
   set maxAng 0.
   if {$xrad < 0 && $yrad < 0} {

#Rectangular field
	set pairs [list [list -1. -1.] [list -1. 1.] [list 1. -1.] \
	   [list 1. 1.] [list 0. 0.]]
   } else {
	set pairs [list [list -1. 0.] [list 1. 0.] [list 0. -1.] \
	   [list 0. 1.] [list 0. 0.]]
	}
   foreach pair $pairs {
	set xmm [expr [lindex $pair 0]*$xrad]
	set ymm [expr [lindex $pair 1]*$yrad]
	set list [focaltosky $hndl $xmm $ymm $ifil]
	set xang [lindex $list 0]
	set yang [lindex $list 1]

#Return angle in degrees
	set mag [expr sqrt(pow($xang,2) + pow($yang,2))/60.]
	set maxAng [expr max($maxAng,$mag)]
	}
   return [format %.3f $maxAng]
   }

###############################################################
#Input is a list of lists (many variations possible).
#Expand any range specifications.  Return a single list.

proc range {args} {
   set outlist ""
   foreach arg $args {
        foreach field $arg {
           if {[regexp {([0-9]+)-([0-9]+)} $field all f1 f2]} {
                loop i $f1 [expr $f2+1] {
                   lappend outlist $i
                   }
           } else {
                lappend outlist $field
                }
           }
        }
   return $outlist
   }

##########################################################################
#Convert a stand-alone design to a mosaic.
#
#In specifying a mosaic, there are two choices for ordering the surfaces ids.
#Either
#   2.01
#   3.01
#   4.01
#   2.02
#   3.02
#   4.02
# etc, or
#   2.01
#   2.02
#   3.01
#   3.02
#   4.01
#   4.02
# etc.
#
#For kent003, I used the former convention.  The latter convention works
#better if I not all surfaces are consecutive -- e.g., 2.01, 3, 4.01, etc.
#
#I believe the code does NOT care about numerical ordering of surface id,
#only that all surfaces for a given color are ordered properly in
#the optic structure.

#Each individual "color" can itself be a list

proc mosaic {hndl surfs colors} {

#First, do some checks before we start.
   foreach surf $surfs {
	if {![ctype digit $surf]} {
	   error "Surface IDs must be integers; $surf is not"
	   }
	}

#Make sure colors exist
   set ncolor [exprGet $hndl.ncolor]
   if {[llength $colors] <= 0} {
	error "Please specified some colors"
	}
   foreach combo $colors {
	set list [range $combo]
	foreach icolor $list {
	   if {$icolor > $ncolor} {
		error "Color $icolor is not in optic structure"
		}
	   }
	}

#For each surface, we create new surfaces with fract part from 1 to
#number of color combos.  Each combo can be 1 or more individual colors
#(e.g., if we specify low and high wavelengths for a filter).

   set nfract [llength $colors]
   foreach surf $surfs {
	set surfId(0) $surf.01
	set isurf [surfIndex $hndl $surf]
	setSurfId $hndl $isurf $surfId(0)
	set newid $surfId(0)
	loop i 1 $nfract {
	   set newid [opticInsert $hndl $newid 2]
	   set surfId($i) $newid
	   opticSurfCopy $hndl $surfId(0) $hndl $newid
	   }

#Zero out refractive index for all but the combo colors
	loop i 0 $nfract {
	   set combo [lindex $colors $i]
	   set list [range $combo]
	   loop j 0 $ncolor {
		if {[lsearch $list $j] < 0} {
		   setIndex $hndl $surfId($i) $j 0
		   }
		}
	   }
	}
   return
   }

