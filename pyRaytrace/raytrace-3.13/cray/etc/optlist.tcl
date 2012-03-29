#Commands:
#	stats
#	surfId <hndl> <index>
#	surfIndex <hndl <surfId>
#	surfTrans <surfId> <paramName>
#	surfTrans <focalParam> <color>
#	setFocal <hndl> <color> <focalName> <value>
#	setFocalFlag <hndl> <color> <focalName> <value>
#	showFocal <hndl> <color> <focalName>
#	showFocalFlag <hndl> <color> <focalName>
#	setSurf <hndl> <surf> <paramName> <value>
#	setSurfFlag <hndl> <surf> <paramName> <value>
#	setColorFlag <hndl> <color> <value>
#	showSurf <hndl> <surf> <paramName>
#	showSurfFlag <hndl> <surf> <paramName>
#	setIndex <hndl> <surf> <color> <value>
#	showIndex <hndl> <surf> <color>
#	showWave <hndl> <icolor>
#	showScale <hndl> <icolor>
#   showparam
#   setflag
#   cacheFlags
#   resetFlags (from cache)
#   clearFlags (clears cache only)
#   resetOptic (from cache)
#   showflag
#   fslist
#   optdiff
#######################################################################
proc stats {} {
   global _optic

#Run the kentools "status" command if it exists
   if {[info command status] != ""} status
   if {[info exists _optic(opticache)]} {
	set cache $_optic(opticache)
   } else {
	set cache ""
	}
   puts stdout "*** Optical Designs ***"
   puts stdout "Hndl Name"
   set handles [handleList]
   foreach handle $handles {
	set hndl [lindex $handle 0]
	set type [lindex $handle 1]
	if {$type != "OPTIC"} continue
	if {$hndl == "$cache"} {
	   set info (cached)
	} else {
	   set info ""
	   }
	echo [format "%-4s %s" $hndl [lindex [exprGet $hndl.name] 0]] $info
	}
   return
   }

#######################################################################
#
#4/3/04: Encode single digit fractions as, e.g., 14.06 instead of 14.6.
#Turns out that in C code it is impossible to distinguish the two cases
#in a floating point number.
#Hope this doesn't break any code!  Watch out.


proc surfName {optic isurf} {
   return [surfId $optic $isurf]
   }

#######################################################################
#
#4/3/04: Encode single digit fractions as, e.g., 14.06 instead of 14.6.
#Turns out that in C code it is impossible to distinguish the two cases
#in a floating point number.
#Hope this doesn't break any code!  Watch out.
#
#Convert surf index to surfId
#My paradigm will be
#   isurf - index in optic array
#   surfId - idint.idfract id of a surface (e.g, 14.06)
#   surfName - text string that a user assigns (e.g., C1)
#I am working on converting all "surfName" variables to "surfId"

proc surfId {optic isurf} {
   set int [exprGet $optic.optic<$isurf>->idint]
   set fract [exprGet $optic.optic<$isurf>->idfract]
   if {$fract == 0} {
	set name $int
   } else {
	set name $int.[format %02d $fract]
	}
   return $name
   }

#######################################################################
#Safer version of surfindex.  Convert id to integer index in optic array.
#If format is, e.g., 11.6, translate to 11.06, which is what is stored
#internally.

proc surfIndex {optic id} {
   set list [split $id .]
   set int [lindex $list 0]
   set fract [fractFormat $id]
   if {[string length $fract] == 0} {
	return [surfindex $optic $int]
	}
   set id $int$fract
   return [surfindex $optic $id]
   }

#######################################################################
#Set a surfid to a new value.  I can only make restricted changes - e.g.,
#change "9" to "9.01".  I will apply tests.

proc surfIdSet {hndl oldid newid} {
   set index [surfIndex $hndl $oldid]
   set index1 [expr $index - 1]
   set index2 [expr $index + 1]
   set nsurf [exprGet $hndl.nsurf]
   if {$index1 >= 0} {
	set id1 [surfId $hndl $index1]
	if {$id1 >= $newid} {
	   echo Previous id $id1 is bigger than new id $newid
	   return
	   }
	}
   if {$index2 <= $nsurf} {
	set id2 [surfId $hndl $index2]
	if {$id2 <= $newid} {
	   echo Next id $id2 is bigger than new id $newid
	   return
	   }
	}
   set list [split $newid .]
   set int [lindex $list 0]
   set fract [lindex $list 1]
   handleSet $hndl.optic<$index>->idint $int
   if {$fract != ""} {
	set fract [string range $fract 0 1]
	regsub ^0 $fract "" fract
   } else {
	set fract 0
	}
   handleSet $hndl.optic<$index>->idfract $fract
   return
   }

#######################################################################
#Rename a surface id.  I specify the surface by its index.

proc setSurfId {optic oldid id} {

   set index [surfIndex $optic $oldid]

#Split the id into integer and fractional part
   set list [split $id .]
   set int [lindex $list 0]
   set fract [lindex $list 1]
   if {$fract == ""} {set fract 0}
   set nsurf [exprGet $optic.nsurf]
   if {$index > $nsurf} {
	error "Index $index greater than nsurf $nsurf"
	}
   handleSet $optic.optic<$index>->idint $int
   handleSet $optic.optic<$index>->idfract $fract
   return
   }

#######################################################################
#Translate the name of a parameter
#For surface (surf > 0) input is
#	"surfid parameter"
#For a focal plane param, input is
#	"param colorindex"

proc surfTrans {surf index} {
   if {![catch {format %f $surf}]} {
	if {![catch {format %d $index}]} {
	   return [list $surf $index]
	   }

#Don't translate index if surf < 0
	if {$surf < 0} {
	   return "$surf $index"
	   }
	if {$index == "curv"} {
		return [list $surf 1]
	} elseif {$index == "ccon"} {
		return [list $surf 2]
	} elseif {$index == "x" || $index == "xoff"} {
		return [list $surf 3]
	} elseif {$index == "y" || $index == "yoff"} {
		return [list $surf 4]
	} elseif {$index == "z" || $index == "zoff"} {
		return [list $surf 5]
	} elseif {$index == "phi"} {
		return [list $surf 6]
	} elseif {$index == "theta"} {
		return [list $surf 7]
	} elseif {$index == "a2" || $index == "A2"} {
		return [list $surf 8]
	} elseif {$index == "a4" || $index == "A4"} {
		return [list $surf 9]
	} elseif {$index == "a6" || $index == "A6"} {
		return [list $surf 10]
	} elseif {$index == "a8" || $index == "A8"} {
		return [list $surf 11]
	} elseif {$index == "a10" || $index == "A10"} {
		return [list $surf 12]
	} elseif {$index == "astig" || $index == "ASTIG"} {
		return [list $surf 15]
	} elseif {$index == "aphi" || $index == "APHI"} {
		return [list $surf 16]
	} elseif {$index == "a1" || $index == "A1"} {
		return [list $surf 21]
	} elseif {$index == "a3" || $index == "A3"} {
		return [list $surf 22]
	} elseif {$index == "a5" || $index == "A5"} {
		return [list $surf 23]
	} elseif {$index == "a7" || $index == "A7"} {
		return [list $surf 24]
	} elseif {$index == "a9" || $index == "A9"} {
		return [list $surf 25]
	} elseif {$index == "a11" || $index == "A11"} {
		return [list $surf 26]
	} elseif {$index == "a13" || $index == "A13"} {
		return [list $surf 27]
	} elseif {$index == "instop"} {
		return [list $surf 401]
	} elseif {$index == "outstop"} {
		return [list $surf 402]
	} elseif {$index == "stoptype"} {
		return [list $surf 403]
	} elseif {$index == "x1"} {
		return [list $surf 404]
	} elseif {$index == "y1"} {
		return [list $surf 405]
	} elseif {$index == "x2"} {
		return [list $surf 406]
	} elseif {$index == "y2"} {
		return [list $surf 407]
	} elseif {$index == "x3"} {
		return [list $surf 408]
	} elseif {$index == "y3"} {
		return [list $surf 409]
	} elseif {$index == "x4"} {
		return [list $surf 410]
	} elseif {$index == "y4"} {
		return [list $surf 411]
	} elseif {$index == "reflect"} {
		return [list $surf 412]
	} elseif {$index == "lines"} {
		return [list $surf 413]
	} elseif {$index == "order"} {
		return [list $surf 414]
	} elseif {$index == "blaze"} {
		return [list $surf 415]
	} else {
		error "Unknown index: $index"
		}
	}
   if {$surf == "xpos" || $surf == "xoff"} {
	return [list -2 $index]
   } elseif {$surf == "ypos" || $surf == "yoff"} {
	return [list -1 $index]
   } elseif {$surf == "xrad" || $surf == "xsize"} {
	return [list -2 [expr $index+100]]
   } elseif {$surf == "yrad" || $surf == "ysize"} {
	return [list -1 [expr $index+100]]
   } elseif {$surf == "scale"} {
	return [list -3 $index]
   } elseif {$surf == "wave"} {
	return [list -3 [expr $index+100]]
   } elseif {$surf == "dist"} {
	return [list -4 $index]
   } elseif {$surf == "map"} {
	return [list -4 [expr $index+100]]
   } elseif {$surf == "rot"} {
	return [list -5 $index]
   } elseif {$surf == "weight"} {
	return [list -5 [expr $index+100]]
   } elseif {$surf == "fl"} {
	return [list -6 $index]
   } elseif {$surf == "focus"} {
	return [list -6 [expr $index+100]]
   } elseif {$surf == "exit"} {
	return [list -7 $index]
   } elseif {$surf == "xmag"} {
	return [list -7 [expr $index+100]]
   } elseif {$surf == "entrance"} {
	return [list -8 $index]
   } elseif {$surf == "emag"} {
	return [list -8 [expr $index+100]]
   } else {
	error "Unknown surface type $surf"
	}
   return "$surf $index"
   }

#######################################################################
#Easier version of setparam for focal plane
proc setFocal {hndl color param value} {
   set list [surfTrans $param $color]
   set param [lindex $list 0]
   set color [lindex $list 1]
   setparam $hndl $param $color $value
   return
   }

#######################################################################
#Easier version of setinc for focal plane
proc setFocalInc {hndl color param value} {
   set list [surfTrans $param $color]
   set param [lindex $list 0]
   set color [lindex $list 1]
   setinc $hndl $param $color $value
   return
   }

#######################################################################
#Easier version of setflag for focal plane
proc setFocalFlag {hndl color param value} {
   set list [surfTrans $param $color]
   set param [lindex $list 0]
   set color [lindex $list 1]
   setflag $hndl $param $color $value
   return
   }

#######################################################################
#Easier version of setflag for focal plane
#Link multiple colors.  I think we use same increment for a given parameter,
#so no need to reset.
proc linkFocalFlag {hndl colors param value} {
   if {$value == 1} {
	echo Warning: Flag of 1 not useful for linking colors!
	}
   foreach color $colors {
	setFocalFlag $hndl $color $param $value
	}
   return
   }

#######################################################################
#Easier version of showparam for focal plane
proc showFocal {hndl color param} {
   set list [surfTrans $param $color]
   set param [lindex $list 0]
   set color [lindex $list 1]
   return [showparam $hndl $param $color]
   }

#######################################################################
#Easier version of showinc for focal plane
proc showFocalInc {hndl color param} {
   set list [surfTrans $param $color]
   set param [lindex $list 0]
   set color [lindex $list 1]
   return [showinc $hndl $param $color]
   }

#######################################################################
#Easier version of showflag for focal plane
proc showFocalFlag {hndl color param} {
   set list [surfTrans $param $color]
   set param [lindex $list 0]
   set color [lindex $list 1]
   return [showflag $hndl $param $color]
   }

#######################################################################
#Easier version of setparam for surfaces.  surf is surface ID.
proc setSurf {hndl surf param value} {
   set list [surfTrans $surf $param]
   set param [lindex $list 1]
   setparam $hndl $surf $param $value
   return
   }

#######################################################################
#Increment the value of a particular parameter
#Allow multiple surfaces (e.g., if we are shifting a lens)
proc incrSurf {hndl surfs param value} {
   foreach surfid $surfs {
	set list [surfTrans $surfid $param]
	set param [lindex $list 1]
	set oldval [showSurf $hndl $surfid $param]
	set newval [expr $oldval + $value]
	setparam $hndl $surfid $param $newval
	}
   return
   }

#######################################################################
#Easier version of setflag for surfaces.  surf is surface ID.
#If flag is > 1, we are linking surfaces.

proc setSurfFlag {hndl surfs param value} {
   foreach surfid $surfs {
	set list [surfTrans $surfid $param]
	set param [lindex $list 1]
	setflag $hndl $surfid $param $value
	return
	}
   }

#######################################################################
#Easier version of setflag for surfaces.  surf is surface ID.
#Also link mulitiple surfaces.  Make sure increment is same for all links
#surfaces

proc linkSurfFlag {hndl surfs param value} {

   if {$value == 1} {
	echo Warning: Flag value of 1 not useful for linking params!
	}

#Use small absolute value of inc.
   set inc [showSurfInc $hndl [lindex $surfs 0] $param]
   loop i 1 [llength $surfs] {
	set newinc [showSurfInc $hndl [lindex $surfs $i] $param]
	if {abs($newinc) < abs($inc)} {
	   set inc $newinc
	   }
	}
   foreach surf $surfs {
	setSurfInc $hndl $surf $param $inc
	}
   foreach surf $surfs {
	setSurfFlag $hndl $surf $param $value
	}
   return
   }

#######################################################################
#Invert the increment for a surface.  This is useful if I want parameters
#"antilinked"

proc invertSurfInc {hndl surf param} {
   setSurfInc $hndl $surf $param [expr -1.*[showSurfInc $hndl $surf $param]]
   return
   }

#######################################################################
#Easier version of setparam for surfaces.
proc showSurf {hndl surf param} {
   set list [surfTrans $surf $param]
   set param [lindex $list 1]
   return [showparam $hndl $surf $param]
   }

#######################################################################
#Easier version of showflag for surfaces.
proc showSurfFlag {hndl surf param} {
   set list [surfTrans $surf $param]
   set param [lindex $list 1]
   return [showflag $hndl $surf $param]
   }

#######################################################################
#Easier version of setinc for surfaces.
proc setSurfInc {hndl surf param value} {
   set list [surfTrans $surf $param]
   set param [lindex $list 1]
   setinc $hndl $surf $param $value
   }

#######################################################################
#Easier version of setinc for indices.
proc setIndexInc {hndl surf param value} {
   set param [expr $param+100]
   setinc $hndl $surf $param $value
   }

#######################################################################
#Easier version of showinc for surfaces.
proc showSurfInc {hndl surf param} {
   set list [surfTrans $surf $param]
   set param [lindex $list 1]
   return [showinc $hndl $surf $param]
   }

#######################################################################
#Easier version of showinc for indexes
proc showIndexInc {hndl surf param} {
   set index [expr $param+100]
   return [showinc $hndl $surf $index]
   }

#######################################################################
#Easier version of setparam for surfaces.  surf is surface ID.
proc setSurf {hndl surf param value} {
   set list [surfTrans $surf $param]
   set param [lindex $list 1]
   setparam $hndl $surf $param $value
   return
   }

#######################################################################
#Easier version of setparam for refraction indexes.  surf is surface ID.
proc setIndex {hndl surf color value} {
   set param [expr 100+$color]
   if {[catch {format %f $value}]} {
	set wave [showWave $hndl $color]
	set value [glass $value $wave]
	}
   setparam $hndl $surf $param $value
   return
   }

#######################################################################
#Easier version of setparam for refraction index gradients. surf is surface ID.
proc setIndexGrad {hndl surf color value} {
   set param [expr 500+$color]
   setparam $hndl $surf $param $value
   return
   }

#######################################################################
#Easier version of setflag for surfaces.  surf is surface ID.
proc setIndexFlag {hndl surf icolor value} {
   set param [expr 100+$icolor]
   setflag $hndl $surf $param $value
   return
   }

#######################################################################
#Easier version of setflag for surfaces.  surf is surface ID.
#Also link mulitiple surfaces.  Make sure increment is same for all links
#surfaces

proc linkIndexFlag {hndl surf colors value} {
   if {$value == 1} {
	echo Warning: flag of 1 not useful for linking indices!
	}

   set inc [showIndexInc $hndl $surf [lindex $colors 0]]
   set wave1 [showWave $hndl [lindex $colors 0]]
   loop i 1 [llength $colors] {
	set newinc [showIndexInc $hndl $surf [lindex $colors $i]]
	if {abs($newinc) < abs($inc)} {
	   set inc $newinc
	   }
	}

#Scale increment by wavelength.  This approximates a line parallel to the
#main glass ridge-line in the Abbe plot

   foreach color $colors {
	set wave [showWave $hndl $color]
	setIndexInc $hndl $surf $color [expr $inc*pow($wave/$wave1,-0.18)]
	setIndexFlag $hndl $surf $color $value
	}
   return
   }

#######################################################################
#Easier version of setparam for refraction indexes.  surf is surface ID.
proc showIndex {hndl surf color} {
   set param [expr 100+$color]
   return [showparam $hndl $surf $param]
   return
   }

#######################################################################
#Easier version of setparam for refraction indexes.  surf is surface ID.
proc showIndexGrad {hndl surf color} {
   set param [expr 500+$color]
   return [showparam $hndl $surf $param]
   return
   }

#######################################################################
#After running indexSet for a single filter, the proper refractive index
#(with sign adjusted) is store in memory.  Some procedures (e.g.,
#lensParam) rely on this.  The following proc fetches the internally
#stored index.

proc showIntIndex {hndl surf} {
   set isurf [surfIndex $hndl $surf]
   return [exprGet $hndl.optic<$isurf>->index]
   return
   }

#######################################################################
#Easier version of showflag for surfaces.
proc showIndexFlag {hndl surf color} {
   set param [expr 100+$color]
   return [showflag $hndl $surf $param]
   }

#######################################################################
#Easier version of setflag for specifying colors to include in fit.
proc setColorFlag {hndl color value} {
   setflag $hndl 0 $color $value
   return
   }

#######################################################################
#Easier version of setflag for specifying colors to include in fit.
proc linkColorFlag {hndl colors value} {
   foreach color $colors {
	setColorFlag $hndl $color $value
	}
   return
   }

#######################################################################
#Easier version of setparam for refraction indexes.  surf is surface ID.
proc showColorFlag {hndl color} {
   set param [expr 100+$color]
   return [showflag $hndl 0 $color]
   }

#######################################################################
#Get wavelength.  Index is modulo 100
proc showWave {hndl index} {
   set list [surfTrans wave [expr int(fmod($index,100))]]
   return [showparam $hndl [lindex $list 0] [lindex $list 1]]
   }

#######################################################################
#Get scale factor.  Index is modulo 100
proc showScale {hndl index} {
   set list [surfTrans scale [expr int(fmod($index,100))]]
   return [showparam $hndl [lindex $list 0] [lindex $list 1]]
   }

#######################################################################
#Set parameters for a lens.  This is like setSurf for all parameters except
#theta, where we need to translate axes as the lens rotates.
#Rotations are not commutative - so I need to set phi before setting theta.
#Need to remember that - no enforcement here.

proc setLensSurf {hndl surf1 surf2 param val} {
   set param [string tolower $param]
   if {$param != "theta"} {
	setSurf $hndl $surf1 $param $val
	setSurf $hndl $surf2 $param $val
	return
	}

#Tweak offsets
   set x1 [showSurf $hndl $surf1 x]
   set x2 [showSurf $hndl $surf2 x]
   set y1 [showSurf $hndl $surf1 y]
   set y2 [showSurf $hndl $surf2 y]
   set z1 [showSurf $hndl $surf1 z]
   set z2 [showSurf $hndl $surf2 z]
   set xavg [expr ($x1+$x2)/2.]
   set yavg [expr ($y1+$y2)/2.]
   set zavg [expr ($z1+$z2)/2.]
   set phi1 [showSurf $hndl $surf1 phi]
   set phi2 [showSurf $hndl $surf2 phi]
   set phi [expr ($phi1+$phi2)/2.]
   set theta1 [showSurf $hndl $surf1 theta]
   set theta2 [showSurf $hndl $surf2 theta]
   set theta [expr ($theta1+$theta2)/2.]

#Vector connecting 2 vertexes.
   set dx [expr $x2-$x1]
   set dy [expr $y2-$y1]
   set dz [expr $z2-$z1]

#Unit vector perpendicular to surface
   set u1 [expr sin($theta)*cos($phi)]
   set u2 [expr sin($theta)*sin($phi)]
   set u3 [expr cos($theta)]

#Thickness, hopefully with sign right.  Assume unit vector is || to vertex vec
   set thick [expr $u1*$dx + $u2*$dy + $u3*$dz]

#Now compute new vertex vector.
   set theta $val
   set dx [expr $thick/2.*sin($theta)*cos($phi)]
   set dy [expr $thick/2.*sin($theta)*sin($phi)]
   set dz [expr $thick/2.*cos($theta)]
   setSurf $hndl $surf1 x [expr $xavg - $dx]
   setSurf $hndl $surf2 x [expr $xavg + $dx]
   setSurf $hndl $surf1 y [expr $yavg - $dy]
   setSurf $hndl $surf2 y [expr $yavg + $dy]
   setSurf $hndl $surf1 z [expr $zavg - $dz]
   setSurf $hndl $surf2 z [expr $zavg + $dz]
   setSurf $hndl $surf1 $param $val
   setSurf $hndl $surf2 $param $val
   return
   }

#######################################################################
# Counterpart of setparam
proc showparam {hndl surfid index} {
   set surf [surfIndex $hndl $surfid]
   if {$surf < 0} then {
	set elem fplane
	set surf [expr -1*$surf]
   } else {
	set elem optic
	}
   if {$index >= 0 && $index < 100} then {
	set array param
	set j [expr $index]
	}
   if {$index >= 100 && $index < 200} then {
	set array n
	set j [expr $index-100]
	}
   if {$index >= 200 && $index < 300} then {
	set array pinc
	set j [expr $index-200]
	}
   if {$index >= 300 && $index < 400} then {
	set array ninc
	set j [expr $index-300]
	}
   if {$index >= 500 && $index < 600} then {
	set array ngrad
	set j [expr $index-500]
	}
   if {$index == 401} {
	return [handleShow $hndl.optic<$surf>->instop]
	}

   if {$index == 402} {
	return [handleShow $hndl.optic<$surf>->outstop]
	}

   if {$index == 403} {
	return [handleShow $hndl.optic<$surf>->stoptype]
	}

   if {$index == 404} {
	return [handleShow $hndl.optic<$surf>->x1]
	}

   if {$index == 405} {
	return [handleShow $hndl.optic<$surf>->y1]
	}

   if {$index == 406} {
	return [handleShow $hndl.optic<$surf>->x2]
	}

   if {$index == 407} {
	return [handleShow $hndl.optic<$surf>->y2]
	}

   if {$index == 408} {
	return [handleShow $hndl.optic<$surf>->x3]
	}

   if {$index == 409} {
	return [handleShow $hndl.optic<$surf>->y3]
	}

   if {$index == 410} {
	return [handleShow $hndl.optic<$surf>->x4]
	}

   if {$index == 411} {
	return [handleShow $hndl.optic<$surf>->y4]
	}
   if {$index == 412} {
	return [handleShow $hndl.optic<$surf>->reflect]
	}
   if {$index == 413} {
	return [handleShow $hndl.optic<$surf>->lines]
	}
   if {$index == 414} {
	return [handleShow $hndl.optic<$surf>->order]
	}
   if {$index == 415} {
	return [handleShow $hndl.optic<$surf>->blaze]
	}
   if {![info exists array]} {error "Unknown index $index"}
   return [handleShow $hndl.$elem<$surf>->$array<$j>]
   }

#######################################################################
#Tcl version of setflag
#Allow multiple surfs and indices
#NEW BEHAVIOR!
#Just cache the flags; don't set bits in the optic structure.
#This way, I can reset the flags and start over if I goof.

proc setflag {args} {
   global _optic

#New behavior:  Only 3 args needed.
   if {[llength $args] == 3} {
	set surfids [lindex $args 0]
	set indexes [lindex $args 1]
	set flag [lindex $args 2]

#For backward compatibility, allow handle of optic structure to be input
   } elseif {[llength $args] == 4} {
	set hndl [lindex $args 0]
	set surfids [lindex $args 1]
	set indexes [lindex $args 2]
	set flag [lindex $args 3]
   } else {
	error "setflag <surfids> <indexes> <flag>"
	}

   if {![info exists _optic(flags)]} {set _optic(flags) ""}
   foreach surfid $surfids {
       foreach iindex $indexes {

#Cache flags
	   lappend _optic(flags) [list $surfid $iindex $flag]
	   }
	}
   return
   }

#########################################################################
#Set all flags in "flags" cache.  This is now an internal routine.

proc flagSet {hndl} {
   global _optic

#First, reset all flags
   indxclr $hndl
   foreach combo $_optic(flags) {
	set surfid [lindex $combo 0]
	set iindex [lindex $combo 1]
	set flag [lindex $combo 2]
	if {[verbose]} {echo Setting flags $combo}

#Translate symbolic names to ints
	set isurf [surfIndex $hndl $surfid]
	set list [surfTrans $isurf $iindex]
	set surf [lindex $list 0]
	set index [lindex $list 1]
	if {$surf < 0} then {
	   set elem fplane
	   set surf [expr -1*$surf]
	} else {
	   set elem optic
	   }
	if {$index >= 0 && $index < 200} then {
	   set array iflag
	   set j [expr $index]
	   handleSet $hndl.$elem<$surf>->$array<$j> $flag
	   }
	}
   return
   }

#######################################################################
#Cache flags.
proc cacheFlags {} {
   global _optic
   if {![info exists _optic(flags)]} {
	echo No flags to cache
	return
	}
   set _optic(flagcache) $_optic(flags)
   unset _optic(flags)
   return
   }

#######################################################################
#Reset flags to their values cached the last time that I ran lstsq .
proc resetFlags {hndl} {
   global _optic
   if {![info exists _optic(flagcache)]} {
	echo No cached flags
	return
	}
   if {[info exists _optic(flags)]} {unset _optic(flags)}
#   set _optic(flags) $_optic(flagcache)
   foreach combo $_optic(flagcache) {
	   set surfid [lindex $combo 0]
	   set index [lindex $combo 1]
	   set flag [lindex $combo 2]
	   if {[verbose]} {echo Setting flags $combo}
	   setflag $hndl $surfid $index $flag
	   }
   return
   }

#######################################################################
#Clear cached flags.  However, I do not actually reset the flag structures.

proc clearFlags {} {
   global _optic
   if {[info exists _optic(flagcache)]} {
	unset _optic(flagcache)
	}
   if {[info exists _optic(flags)]} {
	unset _optic(flags)
	}
   return
   }

#######################################################################
#While we are at it, proc to reset optic structure
proc resetOptic {hndl} {
   global _optic
   if {![info exists _optic(opticache)]} {
	echo No optic cache
	return
	}
   opticCopy $_optic(opticache) $hndl
   return
   }
   
#######################################################################
#Counterpart of setflag
proc showflag {hndl surfid index} {
   set surf [surfIndex $hndl $surfid]
   if {$surf < 0} then {
	set elem fplane
	set surf [expr -$surf]
   } else {
	set elem optic
	}
   if {$index >= 0 && $index < 200} then {
	set array iflag
	set j [expr $index]
	}
   return [handleShow $hndl.$elem<$surf>->$array<$j>]
   }

#######################################################################
#Set parameter increment
proc setinc {hndl surfid index value} {
   set surf [surfIndex $hndl $surfid]
   if {$surf < 0} then {
	set elem fplane
	set surf [expr -$surf]
   } else {
	set elem optic
	}
   if {$index >= 0 && $index < 100} then {
	set array pinc
	set j [expr $index]
	}
   if {$index >= 100} then {
	set array ninc
	set j [expr $index-100]
	}
   handleSet $hndl.$elem<$surf>->$array<$j> $value]
   }

#######################################################################
#Show parameter increment
proc showinc {hndl surfid index} {
   set surf [surfIndex $hndl $surfid]
   if {$surf < 0} then {
	set elem fplane
	set surf [expr -$surf]
   } else {
	set elem optic
	}
   if {$index >= 0 && $index < 100} then {
	set array pinc
	set j [expr $index]
	}
   if {$index > 100} then {
	set array ninc
	set j [expr $index-100]
	}
   return [handleShow $hndl.$elem<$surf>->$array<$j>]
   }

########################################################################
#Little proc to format fractional part

proc fractFormat {id} {
   set list [split $id .]
   set int [lindex $list 0]
   set fract [lindex $list 1]

#Drop anything after 1st 2 digits
   set fract [string range $fract 0 1]
   if {[string length $fract] == 0} {
	return
	}
   regsub {^0([1-9]+)} $fract {\1} fract
   set id .[format %02d $fract]
   return $id
   }

########################################################################
#List all surfaces for a specific filter
proc fslist {hndl filter} {
   global stdout
   if {![info exists stdout]} then {set stdout stdout}
   set findx [expr $filter]
   opticInfo $hndl $filter
   puts $stdout \
"Wave.     X Pos     Y pos    X rad    Y rad   Scale Factor Rot. Dist. Map \
Wgt"
puts $stdout \
[format "%-10.2f%-10.2f%-10.2f%-9.2f%-9.2f%-10.4f%-6.4f %-6.2f%2.0f%5.1f" \
	[showFocal $hndl $findx wave] \
	[showFocal $hndl $findx xoff] \
	[showFocal $hndl $findx yoff] \
	[showFocal $hndl $findx xrad] \
	[showFocal $hndl $findx yrad] \
	[showFocal $hndl $findx scale] \
	[showFocal $hndl $findx rot] \
	[showFocal $hndl $findx dist] \
	[showFocal $hndl $findx map] \
	[showFocal $hndl $findx weight] ]
   puts stdout [format \
	   "Focal length: %.2f Focal ratio %.2f   Exit pupil: %.2f" \
	   [showFocal $hndl $findx fl] \
	   [expr [showFocal $hndl $findx fl] / \
	      [telDiam $hndl]] \
	   [showFocal $hndl $findx exit] ]
   puts $stdout \
"Surf   Curvature    Con Con  Z pos     XOFF    YOFF   PHI     THETA   FILT/INDX"
#Surface 0
   set id [surfId $hndl 0]
   puts $stdout [format \
"%2s%-3s                                                                   %7.4f" \
		[expr int($id)] \
		[fractFormat $id] \
		[handleShow $hndl.optic<0>->n<$findx>]]
   set nsurf [handleShow $hndl.nsurf]
   for {set i 1} {$i <= $nsurf} {incr i} {
	set id [surfId $hndl $i]
	if {[showIndex $hndl $id $findx] == 0} continue
	   puts $stdout [format "%2s%-3s %11.4g %7.3f %9.2f %7.1f %7.1f\
%7.3f %7.3f %3d %7.4f" \
		[expr int($id)] \
		[fractFormat $id] \
		[showSurf $hndl $id curv] \
		[showSurf $hndl $id ccon] \
		[showSurf $hndl $id z] \
		[showSurf $hndl $id x] \
		[showSurf $hndl $id y] \
		[showSurf $hndl $id phi] \
		[showSurf $hndl $id theta] \
		$filter \
		[showIndex $hndl $id $findx] ]
	}

   puts $stdout ""
   puts $stdout \
"Surf        A2         A4         A6         A8         A10   InStop OutStop T"

   for {set i 1} {$i <= $nsurf} {incr i} {
	set id [surfId $hndl $i]
	if {[showIndex $hndl $id $findx] == 0} continue
	puts $stdout [format "%2s%-3s%11.4g%11.4g%11.4g%11.4g%11.4g%8.1f%8.1f%2d" \
		[expr int($id)] \
		[fractFormat $id] \
		[showSurf $hndl $id a2] \
		[showSurf $hndl $id a4] \
		[showSurf $hndl $id a6] \
		[showSurf $hndl $id a8] \
		[showSurf $hndl $id a10] \
		[showSurf $hndl $id instop] \
		[showSurf $hndl $id outstop] \
		[showSurf $hndl $id stoptype]]
	}

   }

######################################################################
#List focal plane parameters

proc focalList {hndl} {
   set format1 "%3s  %5s %5s %5s %5s %5s %6s %6s %6s %3s %4s"
   set format2 "%3d  %5.2f %5.1f %5.1f %5.1f %5.1f %6.2f %6.4f %6.4f %3d %4.2f"
   set ncolor [exprGet $hndl.ncolor]
   puts stdout [format $format1 Fil Wave xoff yoff xrad yrad scale rot dist \
	map wgt]
   for {set i 1} {$i <= $ncolor} {incr i} {
	puts stdout [format $format2 $i [showFocal $hndl $i wave] \
	   [showFocal $hndl $i xoff] \
	   [showFocal $hndl $i yoff] \
	   [showFocal $hndl $i xrad] \
	   [showFocal $hndl $i yrad] \
	   [showFocal $hndl $i scale] \
	   [showFocal $hndl $i rot] \
	   [showFocal $hndl $i dist] \
	   [showFocal $hndl $i map] \
	   [showFocal $hndl $i weight]]
	}
   return
   }

######################################################################
#Difference input structure with initial one from setup

proc optdiff {std hndl} {

echo Difference are of the sense (first - second)
echo Looping through surfaces
   set diff [opticNew]
   opticCopy $std $diff
   for {set i 1} {$i <= [handleShow $hndl.nsurf]} {incr i} {
	loop j 1 13 {
	   handleSet $diff.optic<$i>->param<$j> \
		[expr { [handleShow $hndl.optic<$i>->param<$j>] - \
		[handleShow $std.optic<$i>->param<$j>]}]
	   }
	}
#Now loop through colors
echo Looping through colors
   for {set j 0} {$j < [handleShow $hndl.ncolor]} {incr j} {
	loop i 1 7 {
	   handleSet $diff.fplane<$i>->param<$j> \
		[expr { [handleShow $hndl.fplane<$i>->param<$j>] - \
		[handleShow $std.fplane<$i>->param<$j>]}]
	   }
	}
   return $diff
   }

######################################################################
#Partial debugging routine - list all surfaces by index, id, and name
#If I supply a surface ID, list all filters that utilize this surface

proc surfList {hndl {surfid ""}} {
   if {$surfid == ""} {
	set nsurf [exprGet $hndl.nsurf]
	set format "%5s %7s    %-s"
	puts stdout [format $format isurf surfId Name]
	for {set i 1} {$i <= $nsurf} {incr i} {
	   set id [surfId $hndl $i]
	   set name [showName $hndl $id]
	   puts stdout [format $format $i $id $name]
	   }
	return
	}

#List info for a single surface
   set glass [showGlass $hndl $surfid]
   puts stdout "Surface $surfid ($glass):"
   set ncolor [exprGet $hndl.ncolor]
   for {set ifil 1} {$ifil <= $ncolor} {incr ifil} {
	set index [showIndex $hndl $surfid $ifil]
	if {$index == 0} continue
	puts stdout [format "   Filter %2d  Index %f" $ifil $index]
	}
   return
   }

######################################################################
proc showFlags {} {
   global _optic
   if {![info exists _optic(flags)]} return
   return $_optic(flags)
   }

###################################################################
#I used to use setFlags to transfer flags from internal cache to the OPTIC
#structure.  Now I use it to initialize the internal flag cache.

proc setFlags {flags} {
   global _optic
   set _optic(flags) $flags
   return
   }

