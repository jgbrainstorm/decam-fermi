#Ghost computation functions.
#Ghosting is between surf1 and surf2.  If surf2 is blank, it is assumed
#to be the focal plane.
#CAUTION!!! I get the indices correct for filter "icolor" but not for any
#other filter; however, the other filters are not "zeroed out".

proc ghostDup {optic surf1 surf2 icolor} {

#surf1 is the interior surface, surf2 is normally the focal plane.

#Duplicate trailing surfaces to get full ghost pattern.

   set surflist [surfIdsGet $optic $icolor]
   set np [llength $surflist]

   set indx1 [lsearch $surflist $surf1]
   if {$indx1 < 0} {
	error "Bad surface name $surf1 specified for color $icolor"
	}
#echo indx is $indx

   set goptic [opticNew]
   opticCopy $optic $goptic

   set focalid [lindex $surflist end]

#If we specified surf2, use it instead of focal plane as first reflection.
   if {$surf2 != ""} {
	set focalid $surf2
	set indx2 [lsearch $surflist $focalid]
	if {$indx2 < 0} {error "Reflecting surface $surf2 is not in design"}
	set np [expr $indx2+1]
	}

#Edit surflist to contain inserted surfaces.
   set surflist [lrange $surflist 0 $indx2]

#echo ifocal is $ifocal, focalid is $focalid
#Number of surfaces to insert.  First one is the desired surface itself;
#others are intermediate between desired surface and focal plane; final one
#is focal plane again.

   set nbetween [expr $indx2-$indx1-1]
   if {$nbetween < 0} {
	error "Inserting negative number of surfaces!"
	}
   set ninsert [expr 1 + 2*$nbetween + 1]

#echo nbetween is $nbetween, ninsert is $ninsert
   set id $focalid
   loop i 0 $ninsert {
	set id [opticInsert $goptic $id]
	lappend surflist $id
	}

#Now copy surfaces for reverse traversal
   loop i 0 $nbetween {

#id1 is original surface
#id2 is same surface on our 1st traversal after reflection
#idi1 is the surface preceding the original surface (which provides the
#index on the 1st traversal)

	set id1 [lindex $surflist [expr $indx2 - $i - 1]]
	set idi1 [lindex $surflist [expr $indx2 - $i - 2]]
	set id2 [lindex $surflist [expr $indx2 + $i + 1]]

	opticSurfCopy $goptic $id1 $goptic $id2

#echo Copying $i1 to $i2

#If surface is a mirror, I don't change the index.
	set oldglass [showGlass $goptic $id1]
	if {$oldglass != "mirror"} {
	   setIndex $goptic $id2 $icolor [showIndex $goptic $idi1 $icolor]
	   setGlass $goptic $id2 [showGlass $goptic $idi1]
	   }

#Set stop type to 1 if it is 0 or 2 (aperture stop).  We want to stop
#rays that go way out of bounds.  stoptypes greater than 2 may needs
#additional care here as well.
	set stoptype [showSurf $goptic $id1 stoptype]
	if {$stoptype == 0 || $stoptype == 2} {
	   setSurf $goptic $id2 stoptype 1
	   }
	}

#Now duplicate reflecting surface
   set id2 [lindex $surflist [expr $indx2 + $nbetween + 1]]

#Cache the surface ID.
   opticSurfCopy $goptic $surf1 $goptic $id2

#Make it a mirror
   surfGlassSwitch $goptic $id2 mirror

#Set stop type to 1.  Other stoptypes may require extra thinking here
#(e.g., grating).
   setSurf $goptic $id2 stoptype 1

#Now copy surfaces for retraversal
   loop i 0 $nbetween {
	set id1 [lindex $surflist [expr $indx1 + $i + 1]]
	set id2 [lindex $surflist [expr $indx2 + $i + $nbetween + 2]]

#echo Copying $i1 to $i2
	opticSurfCopy $goptic $id1 $goptic $id2

#Set stop type to 1 if it is 0 or 2
	set stoptype [showSurf $goptic $id1 stoptype]
	if {$stoptype == 0 || $stoptype == 2} {
	   setSurf $goptic $id2 stoptype 1
	   }
	}

#And the focal plane
   set id2 [lindex $surflist [expr $indx2 + $ninsert]]

#echo Focal plane: Copying $ifocal to $i2
   opticSurfCopy $goptic $focalid $goptic $id2

#Reverse the index of the original focal plane.
   surfGlassSwitch $goptic $focalid mirror

   handleSet $goptic.nsurf [expr [exprGet $optic.nsurf] + $ninsert]
   opticInfo $goptic $icolor
   return $goptic
   }

#########################################################################
#Top level proc to create FITS files of ghost images.
#Allow multiple surfaces; create separate file for each surface.
#Use forkEach to speed things up.
#
#This is cpu intensive, npix=64 is recommended.  rayPattern = 3 gives
#coarse sampling of entrance pupil; increase to do better.

proc opticGhost {optic0 surfids {surf2 ""} icolor npix {rayPattern 3}} {
   global FORK
   if {[info exists FORK]} {set proc forkEach} else {set proc foreach}

#Copy optics structure, since I will fiddle with it.
   set optic [opticNew]
   opticCopy $optic0 $optic

   $proc surfid $surfids {
	echo Surface $surfid

#Break out code to create ghost optic structure into a separate proc
	set goptic [ghostDup $optic $surfid $surf2 $icolor]

#Run the ghosting code.  Need to have the same ray pattern for both original
#and expanded optic structures.
	rayPattern $optic $rayPattern 2
	rayPattern $goptic $rayPattern 2
	set reg [ghost $optic $goptic $icolor $npix]
	opticDel $goptic
	hdrInsWithAscii $reg.hdr DESIGN [showDesign $optic]
	hdrInsWithAscii $reg.hdr SURFID $surfid
	hdrInsWithInt $reg.hdr FILTER $icolor
	set file ghost-$surfid-$icolor.fit
	regWriteToFits $reg $file
	regDel $reg
	echo $file
	}
   opticDel $optic
   }

#########################################################################
#Form all surface combinations and compute approximate intensity of
#exit pupil ghost.
#
#usereflect:  If set, use reflection coefficients in design.  If not set,
#assume reflectivity of 1.

proc pupilGhost {optic0 icolor {usereflect ""}} {

   set surfids [surfIdsGet $optic0 $icolor]

#Copy optics structure, since I will fiddle with it.
   set optic [opticNew]
   opticCopy $optic0 $optic

   set n  0

#Field size
   opticInfo $optic $icolor
   set xrad [showFocal $optic $icolor xrad]
   set yrad [showFocal $optic $icolor yrad]

#Focal plane size.
   if {$xrad < 0 && $yrad < 0} {
	set area [expr 4.*$xrad*$yrad]
   } else {
	set area [expr abs(3.1416*$xrad*$yrad)]
	}

   set telDiam [telDiam $optic $icolor]
   set fl [showFocal $optic $icolor fl]

#Incident beamsize (angular diameter)
   set ang [expr 2.*sqrt(abs($xrad*$yrad))/$fl]

   loop i 0 [llength $surfids] {
	set surfid [lindex $surfids $i]
	loop j [expr $i+1] [llength $surfids] {
	   set surf2 [lindex $surfids $j]
	   set goptic [ghostDup $optic $surfid $surf2 $icolor]
	   opticInfo $goptic $icolor
	   set exit [showFocal $goptic $icolor exit]
	   set emag [showFocal $goptic $icolor emag]
	   set xmag [showFocal $goptic $icolor xmag]

#Reflectivity
	   set reflect1 [showSurf $goptic $surfid reflect]
	   set reflect2 [showSurf $goptic $surf2 reflect]
	   if {$xmag == 0.} {set xmag 1.}
	   set exit [showFocal $goptic $icolor exit]

#First, just compute the exit pupil intensity
#m is demagnification ( > 1 for significant intensity)
	   set m [expr $emag/$xmag]
	   if {$usereflect == ""} {
		set scale 1.
	   } else {
		set scale [expr $reflect1*$reflect2]
		}
	   set flux0 [expr $area*pow(2.*$m/$telDiam,2)*$scale/3.14]

#Blurring due to defocus - this is NOT accurate for large defocus
	   set flux [expr $flux0 * 1./(1. + \
		pow($m*$m*$exit*$ang/$telDiam,2))]
	   set index $surfid,$surf2
	   set ghost($n,flux) $flux
	   set ghost($n,index) $index
	   incr n
	   opticDel $goptic
	   }
	}
   opticDel $optic

#Sort ghosts according to inverse intensity and print out
   while {1} {
	set flag 0
	if {$n <= 1} break
	loop i1 1 $n {
	   set i0 [expr $i1-1]
	   set flux0 $ghost($i0,flux)
	   set flux1 $ghost($i1,flux)
	   if {$flux0 <= $flux1} continue
	   set flag 1
	   set index0 $ghost($i0,index)
	   set index1 $ghost($i1,index)
	   set ghost($i0,flux) $flux1
	   set ghost($i0,index) $index1
	   set ghost($i1,flux) $flux0
	   set ghost($i1,index) $index0
	   }
	if {$flag == 0} break
	}
   set format "%5s %5s     %s"
   if {$usereflect == ""} {
	echo Using default reflectivity of unity
   } else {
	echo Using supplied reflectivities
	}
   echo [format $format Surf1 Surf2 "Relative Ghost Intensity"]
   loop i 0 $n {
	set index $ghost($i,index)
	set list [split $index ,]
	set surf1 [lindex $list 0]
	set surf2 [lindex $list 1]
	set flux $ghost($i,flux)
 	if {$flux == 0.} continue
	echo [format $format $surf1 $surf2 [format %.2e $flux]]
	}
   return
   }

#########################################################################
#Form all surface combinations and compute approximate intensity of
#point source ghost.  Actually, I give distance of focus from focal plane.
#
#
#

proc pointGhost {optic0 icolor {usereflect ""}} {

   set surfids [surfIdsGet $optic0 $icolor]

#Copy optics structure, since I will fiddle with it.
   set optic [opticNew]
   opticCopy $optic0 $optic

   set n 0
   loop i 0 [llength $surfids] {
	set surfid [lindex $surfids $i]
	loop j [expr $i+1] [llength $surfids] {
	   set surf2 [lindex $surfids $j]
	   set goptic [ghostDup $optic $surfid $surf2 $icolor]
	   opticInfo $goptic $icolor
	   set focus [showFocal $goptic $icolor focus]
	   set telDiam [telDiam $goptic $icolor]
	   set fl [showFocal $goptic $icolor fl]

#Reflectivity
	   set reflect1 [showSurf $goptic $surfid reflect]
	   set reflect2 [showSurf $goptic $surf2 reflect]
	   if {$usereflect == ""} {
		set scale 1.
	   } else {
		set scale [expr $reflect1*$reflect2]
		}

#Approximate intensity of ghost  - inverse of beam area.
	   set inten [expr $scale/pow($focus*($telDiam/$fl),2)]
	   set index $surfid,$surf2
	   set ghost($n,flux) $inten
	   set ghost($n,index) $index
	   incr n
	   opticDel $goptic
	   }
	}
   opticDel $optic

#Sort ghosts according to inverse intensity and print out
   while {1} {
	set flag 0
	if {$n <= 1} break
	loop i1 1 $n {
	   set i0 [expr $i1-1]
	   set flux0 $ghost($i0,flux)
	   set flux1 $ghost($i1,flux)
	   if {$flux0 >= $flux1} continue
	   set flag 1
	   set index0 $ghost($i0,index)
	   set index1 $ghost($i1,index)
	   set ghost($i0,flux) $flux1
	   set ghost($i0,index) $index1
	   set ghost($i1,flux) $flux0
	   set ghost($i1,index) $index0
	   }
	if {$flag == 0} break
	}
   set format "%5s %5s     %s"
   if {$usereflect == ""} {
	echo Using default reflectivity of unity
   } else {
	echo Using supplied reflectivities
	}
   echo [format $format Surf1 Surf2 "Relative Ghost Intensity"]
   loop i 0 $n {
	set index $ghost($i,index)
	set list [split $index ,]
	set surf1 [lindex $list 0]
	set surf2 [lindex $list 1]
	set flux $ghost($i,flux)
 	if {$flux == 0.} continue
	echo [format $format $surf1 $surf2 [format %.2e $flux]]
	}
   return
   }





