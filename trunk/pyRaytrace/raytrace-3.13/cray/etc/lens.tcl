#Some lens manipulation parameters

#Insert a lens into a design.  Specify preceding surface, z position of first
#surface, thickness, and glass type.
#I do NOT copy x, y positions or tilts.

proc lensInsert {optic fsurf zpos zthick glass {name ""}} {

   if {$name == ""} {set name LENS}

#Get surface index
   set fsurf1 [opticInsert $optic $fsurf]
   set fsurf2 [opticInsert $optic $fsurf1]

   setSurf $optic $fsurf1 z $zpos
   setSurf $optic $fsurf2 z [expr $zpos + $zthick]

#Initial curvatures, etc. are 0
#Index of refraction:  Loop through colors, check index of preceding
#surface.  If it is 0, don't set index for new lens  Otherwise, set index
#of surface 1 based on the wavelength; set the index of surface 2 to be
#the same as the preceding surface.  Be sure to get signs right.

   set ncolor [exprGet $optic.ncolor]
   set oldglass [showGlass $optic $fsurf]
   for {set i 1} {$i <= $ncolor} {incr i} {
	set oldindex [showIndex $optic $fsurf $i]
	if {$oldindex == 0} continue
	set oldglass [showGlass $optic $fsurf]

#If index is -1, need to find index of first non-mirror surface before this
#one.  Use brute strength and ignorance.
	if {$oldindex == -1} {
	   set surfids [surfIdsGet $optic $i]
	   set np [lindex $surfids $fsurf]
	   loop j 0 $np {
		set surfid [lindex $surfids $j]
		set indx [showIndex $optic $surfid $i]
		if {$indx != -1} {
		   set oldindex $indx
		   set oldglass [showGlass $optic $surfid]
		   }
		}
	   }
	set wave [showFocal $optic $i wave]
	set index1 [glass $glass $wave]
	setIndex $optic $fsurf1 $i $index1
	setIndex $optic $fsurf2 $i $oldindex
	}

   setGlass $optic $fsurf1 $glass
   setGlass $optic $fsurf2 $oldglass
   setName $optic $fsurf1 $name
   setName $optic $fsurf2 $name

#Rerun stopcomp, opticInc
   stopcomp $optic
   opticInc $optic 1
   return
   }

############################################################################
#Insert a surface into a design.  Specify preceding surface, z position of
#surface, and glass type.
#I do NOT copy x, y positions or tilts.
#This does something more useful than just calling opticInsert

proc surfInsert {optic fsurf zpos glass {mode 0} {name ""}} {

   if {$name == ""} {set name SURFACE}

#Get surface index
   set isurf [surfIndex $optic $fsurf]
   opticInsert $optic $fsurf $mode

#Get new surface name
   set isurf1 [expr $isurf + 1]
   set fsurf1 [surfId $optic $isurf1]

   setSurf $optic $fsurf1 z $zpos
   setGlass $optic $fsurf1 $glass
   setName $optic $fsurf1 $name

#Initial curvatures, etc. are 0
#Index of refraction:  Loop through colors, check index of preceding
#surface.  If it is 0, don't set index for new lens  Otherwise, set index
#of surface 1 based on the wavelength; set the index of surface 2 to be
#the same as the preceding surface.  Be sure to get signs right.

   set ncolor [exprGet $optic.ncolor]
   for {set i 1} {$i <= $ncolor} {incr i} {
	set index [showIndex $optic $fsurf $i]
	if {$index == 0} continue
	set wave [showFocal $optic $i wave]
	set index1 [glass $glass $wave]
	setIndex $optic $fsurf1 $i $index1
	}

#Rerun stopcomp, opticInc
   stopcomp $optic
   opticInc $optic 1
   return $fsurf1
   }

	

########################################################################
#Rotate a set of surfaces.
#For now, I will restrict rotation to be about the y axis at a point
#on the z axis.
#Generalization - let's add xpivot as well.  Useful for dealing with
#spectrograph designs where gratings and prisms introduce multiple
#coordinate breaks.
#
#I need to modify two sets of entries:
#   a) The location of the vertex
#   b) The orientation angles
#I assume each has some previous value already.
#rot is in degrees.

proc surfRotate {hndl fsurfs xpivot zpivot rot} {

   foreach fsurf $fsurfs {

#First, orientation of angles
	foreach param {x y z theta phi} {
	   set $param [showSurf $hndl $fsurf $param]
	   }

#theta, phi are already in radians
	set rtheta $theta
	set rphi $phi
	set PI 3.14159265358979
	set RADIAN [expr $PI/180.]
	set rrot [expr $rot*$RADIAN]
	set mux [expr sin($rtheta)*cos($rphi)]
	set muy [expr sin($rtheta)*sin($rphi)]
	set muz [expr cos($rtheta)]
	set mxp [expr $mux*cos($rrot) + $muz*sin($rrot)]
	set myp [expr $muy]
	set mzp [expr $muz*cos($rrot) - $mux*sin($rrot)]
	if {$mzp > 1.} {set mzp 1.}
	if {$mzp < -1.} {set mzp -1.}
	set theta [expr acos($mzp)]
	set phi [expr atan2($myp, $mxp)]
	if {$theta == 0.} {set phi 0.}

#Bring phi back to -pi/2 to pi/2
	if {$phi > $PI/2.} {
	   set phi [expr $phi-$PI]
	   set theta [expr -1.*$theta]
	   }

	if {$phi < -1.*$PI/2.} {
	   set phi [expr $phi+$PI]
	   set theta [expr -1.*$theta]
	   }

#Now translate the vertex
	set mux [expr $x-$xpivot]
	set muy $y
	set muz [expr $z-$zpivot]
	set mxp [expr $mux*cos($rrot) + $muz*sin($rrot)]
	set myp [expr $muy]
	set mzp [expr $muz*cos($rrot) - $mux*sin($rrot)]
	set x [expr $mxp + $xpivot]
	set y $myp
	set z [expr $mzp + $zpivot]
	foreach param {x y z theta phi} {
	   setSurf $hndl $fsurf $param [set $param]
	   }
	}
   return
   }

####################################################################

#Rotate a lens.  I will assume zpivot is the average z of all the lenses.
#Otherwise, I just call surfRotate
#Rotation about y-axis only.

proc lensRotate {hndl fsurfs rot} {
   if {$fsurfs == ""} return
   set xavg 0.
   set zavg 0.
   foreach fsurf $fsurfs {
	set xavg [expr $xavg + [showSurf $hndl $fsurf x]]
	set zavg [expr $zavg + [showSurf $hndl $fsurf z]]
	}
   set xavg [expr 1.*$xavg/[llength $fsurfs]]
   set zavg [expr 1.*$zavg/[llength $fsurfs]]
   surfRotate $hndl $fsurfs $xavg $zavg $rot
   return
   }

###################################################################

#Rotate all surfaces after a given one.  This is analogous to a Zemax
#coordinate break.  I specify the first surface id.  I use its vertex
#as the pivot point.  I need to input a filter id in order to get a 
#list of surfaces to rotate.

#rot is in degrees.

#If I want to rotate the pivot surface itself, I should call "lensRotate"
#as well.  Easier than trying to add another flag to surfBreak.

proc surfBreak {hndl surfid0 ifil rot} {
   set surfids [surfIdsGet $hndl $ifil]
   set ipivot [lsearch $surfids $surfid0]
   if {$ipivot < 0} return
   set xpivot [showSurf $hndl $surfid0 x]
   set zpivot [showSurf $hndl $surfid0 z]
   set surflist [lrange $surfids [expr $ipivot+1] end]
   surfRotate $hndl $surflist $xpivot $zpivot $rot
   return
   }
