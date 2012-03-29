#Compute incidence angle info.

#########################################################
#
#Angle of incidence deviations for spectroscopy.
#Compute and plot the incidence angles for a set of spots in one filter
#Hmmm, in DECAM the filter is a different surface in the optical beam than
#the focal plane, so I would like the incidence angle on an arbitrary surface.

proc incidence {hndl filter {surfid ""}} {
   upvar 1 xtracelist xtracelist
   upvar 1 xofflist xofflist
   upvar 1 yofflist yofflist
#

#If I specify a surface, copy the optic structure and delete all surfaces
#after this one.  This is a hack that lets me use existing infrastructure
#to get the incidence angle without a lot of extra work.

   if {$surfid != ""} {

#I input surface name.  I need to get name of surface after this one.
	set nsurf [exprGet $hndl.nsurf]
	set isurf [surfIndex $hndl $surfid]
	set hndl1 [opticNew]
	opticCopy $hndl $hndl1
	set hndl $hndl1
	if {$isurf < $nsurf} {
	   incr isurf
	   for {set i $isurf} {$i <= $nsurf} {incr i} {
		set id [surfId $hndl $isurf]
		opticRemove $hndl $id
		}
	   }
	}
   set xtracelist ""
   set xofflist ""
   set yofflist ""
   set xmax [showFocal $hndl $filter xrad]
   set ymax [showFocal $hndl $filter yrad]
   if {$xmax < 0 && $ymax < 0} then {

#Rectangular field
	set rmax [sqrt [expr $xmax*$xmax+$ymax*$ymax]]
   } else {

#Circular (or at least elliptical) field
	set rmax [max $xmax $ymax]
	}
   if {[pgQinf STATE] == "CLOSED"} then {pgBegin; pgAsk 0}
   pgEnv 0 $rmax 0 20. 0 0
   pgLabel "Radius (mm)" "Incidence angle (deg)" [showDesign $hndl]

#Alternative computation should get same answer.
   loop i 0 25 {
        set rad [expr $i*$rmax/24.]
        rtrace $hndl $rad 0 $filter 0

#Angle of incidence (tilt between chief ray and normal to focal plane) in
#radians.
	set incidence [handleShow $hndl.diagram->incidence]
	set ang [expr $incidence*57.3]
        pgPoint $rad $ang 3
	lappend xtracelist $rad
	lappend xofflist $incidence
	lappend yofflist 0
	lappend xtracelist -$rad
	lappend xofflist [expr -1*$incidence]
	lappend yofflist 0
	}
   if {$surfid != ""} {opticDel $hndl}
   return
   }

#########################################################
#
#Angle of incidence deviations for a single ray and surface.
#For filters, I want the full range of incidence angles.

proc rayIncidence {hndl xmm ymm xfract yfract filter surfid} {

#If I specify a surface, copy the optic structure and delete all surfaces
#after this one.  This is a hack that lets me use existing infrastructure
#to get the incidence angle without a lot of extra work.

   if {$surfid != ""} {

#I input surface id.  I need to get id of surface after this one.
	set nsurf [exprGet $hndl.nsurf]
	set isurf [surfIndex $hndl $surfid]
	set hndl1 [opticNew]
	opticCopy $hndl $hndl1
	set hndl $hndl1
	if {$isurf < $nsurf} {
	   incr isurf
	   for {set i $isurf} {$i <= $nsurf} {incr i} {
		set id [surfId $hndl $isurf]
		opticRemove $hndl $id
		}
	   }
	}

   ray $hndl $xmm $ymm $xfract $yfract $filter 0

#diagram.slope is the output slope at the last surface.
#vn is the normal to the lens.  Both are in world coordinates.
#Grr, I need the input ray slope.  I will recalculate from the values
#cached in xray and yray.
   set np [exprGet $hndl.diagram->np]

#Index in xray, yray of last surface
   set i [expr $np-1]

#One less
   set i1 [expr $i-1]
   foreach  var "x y z" {
	set $var [exprGet $hndl.diagram->${var}ray<$i>]
	set ${var}1 [exprGet $hndl.diagram->${var}ray<$i1>]
	}

#Compute vector slope of incoming ray
   set v1 [expr $x-$x1]
   set v2 [expr $y-$y1]
   set v3 [expr $z-$z1]
   set amp [expr sqrt($v1*$v1 + $v2*$v2 + $v3*$v3)]
   set v1 [expr $v1/$amp]
   set v2 [expr $v2/$amp]
   set v3 [expr $v3/$amp]

#Fetch normal
   set n1 [exprGet $hndl.diagram->vn<0>]
   set n2 [exprGet $hndl.diagram->vn<1>]
   set n3 [exprGet $hndl.diagram->vn<2>]

#Compute amplitude of cross-product
   set xprod [expr sqrt(pow($v2*$n3 - $v3*$n2,2) + pow($v3*$n1 - $v1*$n3,2) \
	+ pow($v1*$n2 - $v2*$n1,2))]
   set ang [expr asin($xprod)*57.3]
   if {$surfid != ""} {opticDel $hndl}
   return [format %.3f $ang]
   }

###############################################################

#Compute fractional wavelength shift for an input angle in degrees.
proc tiltShift {ang} {

#Effective index of refraction.  This matches Barr designed filters.
   set n 1.68

#The equation given by Aldering looks bogus.  I will use my own guess at
#the right answer.

#   set shift [expr 1./sqrt(1. - pow(sin($ang/57.3)/$n,2))]

#Ooops, completely blew it.  It goes in the opposite sense.  Aldering was
#right.
   set shift [expr sqrt(1. - pow(sin($ang/57.3)/$n,2))]
   return $shift
   }

#########################################################################
#Compute average change in filter wavelength for a given spot position.
#Average over entrance pupil.

proc filterShift {hndl xmm ymm filter surfid} {
   if {$surfid != ""} {

#I input surface id.  I need to get id of surface after this one.
	set nsurf [exprGet $hndl.nsurf]
	set isurf [surfIndex $hndl $surfid]
	set hndl1 [opticNew]
	opticCopy $hndl $hndl1
	set hndl $hndl1
	if {$isurf < $nsurf} {
	   incr isurf
	   for {set i $isurf} {$i <= $nsurf} {incr i} {
		set id [surfId $hndl $isurf]
		opticRemove $hndl $id
		}
	   }
	}

   set avg 0.
   set n 0
   set nfract [exprGet $hndl.diagram->nfract]

   loop i 0 $nfract {
	set xfract [exprGet $hndl.diagram->xfract<$i>]
	set yfract [exprGet $hndl.diagram->yfract<$i>]
	if {![ray $hndl $xmm $ymm $xfract $yfract $filter 1]} continue

	set ang [rayIncidence $hndl $xmm $ymm $xfract $yfract $filter $surfid]
	incr n
	set z [tiltShift $ang]
	set dz [expr $z-1.]
	set avg [expr $avg + $dz]
	}
   if {$n > 0} {
	set avg [expr $avg/$n]
	}
   if {$surfid != ""} {
	opticDel $hndl
	}
   return [format %.4f $avg]
   }
