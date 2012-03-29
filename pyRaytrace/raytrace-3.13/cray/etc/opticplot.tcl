#TCL version of FORTRAN routine "opticplot".  Draw side view of
#optics and some same ray traces.  Enhancement: Specify which filter.

########################################################################
#First, a procedure to convert from lens-o-centric to global coords.
#For now, I convert just x,z pairs.

proc lens2Global {optic xs zs surfid} {

#Don't use phi for now, which is wrong.
# set phi [showSurf $optic $name phi]
   set theta [showSurf $optic $surfid theta]
   set xvert [showSurf $optic $surfid xoff]
#   set yvert [showSurf $optic $surfid yoff]
   set zvert [showSurf $optic $surfid zoff]
   if {$theta != 0} {
	set x1 [expr $xs*cos($theta) + $zs*sin($theta)]
	set z1 [expr -$xs*sin($theta) + $zs*cos($theta)]
	set xs $x1
	set zs $z1
	}
   set x [expr $xs + $xvert]
   set z [expr $zs + $zvert]
   return [list $x $z]
   }

#########################################################
#A handy procedure to get a list of surfindex for a given filter

proc surfIndexGet {optic filter} {
   set isurfs ""
   set nsurf [exprGet $optic.nsurf]
   for {set isurf 0} {$isurf <= $nsurf} {incr isurf} {
	set id [surfId $optic $isurf]
	if {[showIndex $optic $id $filter] == 0} continue
	lappend isurfs $isurf
	}
   return $isurfs
   }

#########################################################
#A handy procedure to get a list of surfids for a given filter

proc surfIdsGet {optic filter} {
   set surfids ""
   set nsurf [exprGet $optic.nsurf]
   for {set isurf 0} {$isurf <= $nsurf} {incr isurf} {
	set id [surfId $optic $isurf]
	if {[showIndex $optic $id $filter] == 0} continue
	lappend surfids $id
	}
   return $surfids
   }

###########################################################################
#Make a plot of the lens design and run some rays through it.
#I use the routine "rayPupil" which locates center of entrance pupil if
#it is severely displaced from the paraxial location.  This routine is
#not needed if I have run "pupilInfo", but I may not have done so
#when opticPlot is first called.

proc opticPlot {optic nmin nmax filter {noray ""}} {
   global XMIN XMAX YMIN YMAX
   set nsurf [exprGet $optic.nsurf]
   set ncolor [exprGet $optic.ncolor]
#   set nmin [expr max($nmin,1)]
#   set nmax [expr min($nmax,$nsurf)]
   set rmax 0.
   set zmin 1.e14
   set zmax -1.e14

   set opticName [lindex [exprGet $optic.name] 0]

#Aperture stop surface
   set apsurf ""

   set surfids [surfIdsGet $optic $filter]

#Specify limits based only on surfaces actually being plotted.
   foreach surfId $surfids {
	if {$surfId < $nmin || $surfId > $nmax} continue
	set xvert [showSurf $optic $surfId xoff]
	set yvert [showSurf $optic $surfId yoff]
	set rvert [expr sqrt(pow($xvert,2) + pow($yvert,2))]
	set zvert [showSurf $optic $surfId zoff]
	set outstop [showSurf $optic $surfId outstop]
	set zmin [expr min($zmin,$zvert)]
	set zmax [expr max($zmax,$zvert)]

#If aperture center is off-center, make rmax bigger.
	set rmax [expr max($rmax,[expr abs($outstop) + $rvert])]
	}

#Get exit pupil.  This is location of exit pupil relative to the focal plane.
#Exit pupil calculation is now broken, alas.

   set exit [showFocal $optic $filter exit]
#   set zmin [expr min($zmin,$exit)]
#   set zmax [expr max($zmax,$exit)]
   set zrange [expr $zmax - $zmin]
   set zmid [expr ($zmin + $zmax)/2.]

#Pad limits.
   set rmax [expr $rmax*1.2]
   set rmin [expr -1.*$rmax]
   set zmin [expr $zmid - 0.6 * $zrange]
   set zmax [expr $zmid + 0.6 * $zrange]
   if {$rmax-$rmin < $zmax-$zmin} {
	set rmax [expr ($zmax-$zmin)/2.]
	set rmin [expr -1.*$rmax]
	}

   if {[info exists XMIN]} {set rmin $XMIN}
   if {[info exists XMAX]} {set rmax $XMAX}
   if {[info exists YMIN]} {set zmin $YMIN}
   if {[info exists YMAX]} {set zmax $YMAX}

#Begin plot
   plotInit b
   set color [pgQci]
   pgSci 1
   pgEnv $rmin $rmax $zmax $zmin 1 0
   pgLabel Radius Z [lindex $opticName 0]

#Aperture stop surface
   set apsurf ""

   foreach var "xleft zleft xright zright" {
	global $var
	if {[info exists $var]} {unset $var}
	}

   foreach surfId $surfids {
	if {$surfId < $nmin || $surfId > $nmax} continue
	set isurf [surfIndex $optic $surfId]
	set instop [showSurf $optic $surfId instop]
	set outstop [showSurf $optic $surfId outstop]
	set stoptype [showSurf $optic $surfId stoptype]
	if {$stoptype == 0 || $stoptype == 5 || $stoptype == 6} {
	   pgSci 1
	} elseif {$stoptype == 1} {
	   pgSci 3
	   set apsurf $surfId
	} else {
	   pgSci 2
	   }
	if {$stoptype == 2} {
	   set apsurf $surfId
	   }
	set xxl ""
	set zzl ""
	set xxr ""
	set zzr ""
	for {set j -20} {$j <= 20} {incr j} {
	   set r2 [expr pow($outstop * $j/20.,2)]
	   set xs [expr sqrt($r2)]
	   if {$xs < $instop} continue
	   set zs [zsurf $optic $isurf $xs 0.]
	   if {$j < 0} {set xs [expr -1.*$xs]}
	   if {$j == -20} {
		set xleft($surfId) $xs
		set zleft($surfId) $zs
		}
	   if {$j == 20} {
		set xright($surfId) $xs
		set zright($surfId) $zs
		}
	   set list [lens2Global $optic $xs $zs $surfId]
	   set x [lindex $list 0]
	   set z [lindex $list 1]
	   if {$j <= 0} {
		lappend xxl $x
		lappend zzl $z
		}
	   if {$j >= 0} {
		lappend xxr $x
		lappend zzr $z
		}
	   }
	pgLine $xxl $zzl
	pgLine $xxr $zzr
	}

#Get fancy - fill in edges of lenses.
#Each side of lens has 3 lines - the side itself and the top, bottom
#lines connecting the side to the outer stop of the lens.  For the smaller
#of the two surfaces, this line forms a "shelf" extending the surface out
#to the side itself.
   pgSci 1
   loop k 0 [llength $surfids] {
	if {$k == 0 || $k == 1} continue
	set id1 [lindex $surfids [expr $k-1]]
	set id [lindex $surfids $k]
	if {$id < $nmin || $id > $nmax} continue
	if {$id1 < $nmin || $id1 > $nmax} continue
	set outstop1 [showSurf $optic $id1 outstop]
	set outstop [showSurf $optic $id outstop]
	set index1 [showIndex $optic $id1 $filter]
	set index [showIndex $optic $id $filter]

#To identify lenses, I could also match up surfName, not just indexes.
	if {$index1*$index <= 0} continue
	if {abs($index1) <= 1} continue
	set minleft [expr min($xleft($id1),$xleft($id))]
	set maxright [expr max($xright($id1),$xright($id))]

#Still need to convert from lens-o-centric to global
#	pgLine "$minleft $xleft($name1)" "$zleft($name1) $zleft($name1)"
	set list [lens2Global $optic $minleft $zleft($id1) $id1]
	set x1 [lindex $list 0]
	set z1 [lindex $list 1]
	set list [lens2Global $optic $xleft($id1) $zleft($id1) $id1]
	set x2 [lindex $list 0]
	set z2 [lindex $list 1]
	pgLine "$x1 $x2" "$z1 $z2"

#	pgLine "$maxright $xright($name1)" "$zright($name1) $zright($name1)"
	set list [lens2Global $optic $maxright $zright($id1) $id1]
	set x1 [lindex $list 0]
	set z1 [lindex $list 1]
	set list [lens2Global $optic $xright($id1) $zright($id1) $id1]
	set x2 [lindex $list 0]
	set z2 [lindex $list 1]
	pgLine "$x1 $x2" "$z1 $z2"

#	pgLine "$minleft $xleft($name)" "$zleft($name) $zleft($name)"
	set list [lens2Global $optic $minleft $zleft($id) $id]
	set x1 [lindex $list 0]
	set z1 [lindex $list 1]
	set list [lens2Global $optic $xleft($id) $zleft($id) $id]
	set x2 [lindex $list 0]
	set z2 [lindex $list 1]
	pgLine "$x1 $x2" "$z1 $z2"

#	pgLine "$maxright $xright($name)" "$zright($name) $zright($name)"
	set list [lens2Global $optic $maxright $zright($id) $id]
	set x1 [lindex $list 0]
	set z1 [lindex $list 1]
	set list [lens2Global $optic $xright($id) $zright($id) $id]
	set x2 [lindex $list 0]
	set z2 [lindex $list 1]
	pgLine "$x1 $x2" "$z1 $z2"

#	pgLine "$minleft $minleft" "$zleft($name1) $zleft($name)"
	set list [lens2Global $optic $minleft $zleft($id1) $id1]
	set x1 [lindex $list 0]
	set z1 [lindex $list 1]
	set list [lens2Global $optic $minleft $zleft($id) $id]
	set x2 [lindex $list 0]
	set z2 [lindex $list 1]
	pgLine "$x1 $x2" "$z1 $z2"

#	pgLine "$maxright $maxright" "$zright($name1) $zright($name)"
	set list [lens2Global $optic $maxright $zright($id1) $id1]
	set x1 [lindex $list 0]
	set z1 [lindex $list 1]
	set list [lens2Global $optic $maxright $zright($id) $id]
	set x2 [lindex $list 0]
	set z2 [lindex $list 1]
	pgLine "$x1 $x2" "$z1 $z2"
	}

#Plot exit pupil
   pgSci 7
   set name [surfId $optic $nmax]
   set r1 [expr $rmax*3./4.]
   set list [lens2Global $optic -$rmax $exit $id]
   set x1 [lindex $list 0]
   set z1 [lindex $list 1]
   set list [lens2Global $optic -$r1 $exit $id]
   set x2 [lindex $list 0]
   set z2 [lindex $list 1]
   pgLine "$x1 $x2" "$z1 $z2"

   set list [lens2Global $optic $rmax $exit $id]
   set x1 [lindex $list 0]
   set z1 [lindex $list 1]
   set list [lens2Global $optic $r1 $exit $id]
   set x2 [lindex $list 0]
   set z2 [lindex $list 1]
   pgLine "$x1 $x2" "$z1 $z2"


#Plot a few rays
   if {$noray != ""} return
   set PI 3.141593

#Size of focal plane
   set xmm [showFocal $optic $filter xsize]
   set ymm [showFocal $optic $filter ysize]
   if {$xmm < 0 || $ymm < 0} {
	set rad [expr sqrt(pow($xmm,2) + pow($ymm,2))]
   } else {
	set rad $xmm
	}

#Finner is a useful surrogate for knowing the inner stop of the telescope.
#Normally I don't set this in the aperture stop surface, since it is set
#by the secondary mirror.

   set finner [exprGet (*$optic.tel).finner]

#Place spots at center of field and at max. x position.  Keep y=0.
   foreach ixspot "0 1" {
	if {$ixspot == 0} {
	   pgSci 5
	} else {
	   pgSci 6
	   }
	set x [expr $ixspot*$rad]
	set y 0

#Trace a couple of rays for each position
	foreach frad "-1 -$finner $finner 1" {
	   rayPupil $optic $x $y $frad 0 $filter 0
	   set np [exprGet (*$optic.diagram).np]
	   set xx ""
	   set zz ""
	   loop i 0 $np {
		lappend xx [exprGet (*$optic.diagram).xray<$i>]
		lappend zz [exprGet (*$optic.diagram).zray<$i>]
		}
	   pgLine $xx $zz
	   }
	}
   pgSci $color
   return
   }

#############################################################
#Plot outline of baffle

proc bafflePlot {chain} {
   set color [pgQci]
   pgSci 3
   chainEach hndl $chain {
	set type [exprGet -enum $hndl.type]
	if {$type == "ANNULUS"} {
	   foreach var "r1 r2 z" {
		set $var [exprGet $hndl.$var]
		}
	   pgLine "$r1 $r2" "$z $z"
	   set r1 [expr -$r1]
	   set r2 [expr -$r2]
	   pgLine "$r1 $r2" "$z $z"
	} else {
	   foreach var "m z0 z1 z2" {
		set $var [exprGet $hndl.$var]
		}
	   set r1 [expr $m*($z1-$z0)]
	   set r2 [expr $m*($z2-$z0)]
	   pgLine "$r1 $r2" "$z1 $z2"
	   set r1 [expr -$r1]
	   set r2 [expr -$r2]
	   pgLine "$r1 $r2" "$z1 $z2"
	   }
	}
   pgSci $color
   return
   }

###########################################################################
#Analyze physical parameters of lenses
#margin is fractional increase in radius of the surface to provide for margin
#of error.  The surface height is adjusted accordingly.
#rim is the increase in radius (mm) to proivde room for the support system.
#The height is kept fixed.

#According to Bigelow, a rim of 1 cm (radius) is fine for an element 1 meter
#in diameter.

proc lensParam {optic filter {margin 0} {rim 0}} {


   set nsurf [exprGet $optic.nsurf]
   set ncolor [exprGet $optic.ncolor]

#Trace chief ray and get list of surfaces
   rayPupil $optic 0 0 0 0 $filter 0
   set surfids [surfIdsGet $optic $filter]
   puts stdout "Using margin of [format %.4f $margin]"
   puts stdout "Using rim width of $rim (mm)"
   set format1 "%7s %7s %7s %12s %8s %8s %7s %7s %6s"
   puts stdout [format $format1 Surfs  Diam Height Thin/Thick \
	"D+rim" "Vol(cc)" Wgt(lb) Blnk(kg) Glass]

#Link surfaces to construct lenses.
   indexSet $optic $filter
   loop k 0 [llength $surfids] {
	if {$k == 0} continue
	set id1 [lindex $surfids [expr $k-1]]
	set id [lindex $surfids $k]

#Fetch internal index
	set index1 [showIntIndex $optic $id1]
	set index [showIntIndex $optic $id]
	set glass1 [showGlass $optic $id1]
	set glass [showGlass $optic $id]

#If we reverse direction, don't call it a lens.
	if {$index1*$index <= 0} continue

#If first index is not air or vacuum, continue
	if {abs($index1) <= 1} continue

#Store shape of surface in lens-o-centric coords.
	set xx1 ""
	set zz1 ""
	set xx ""
	set zz ""
	set diffz ""

	set outstop1 [abs [showSurf $optic $id1 outstop]]
	set outstop1 [expr $outstop1*(1.+$margin)]
	set outstop [abs [showSurf $optic $id outstop]]
	set outstop [expr $outstop*(1.+$margin)]
	set outstopmax [expr max($outstop1,$outstop)]

#Get separation of lenses.  Use vertices
	set xvert1 [showSurf $optic $id1 x]
	set yvert1 [showSurf $optic $id1 y]
	set zvert1 [showSurf $optic $id1 z]

	set xvert [showSurf $optic $id x]
	set yvert [showSurf $optic $id y]
	set zvert [showSurf $optic $id z]

	set sep [expr sqrt(pow($xvert-$xvert1,2) + pow($yvert-$yvert1,2) \
	   + pow($zvert-$zvert1,2))]

	set isurf1 [surfIndex $optic $id1]
	set isurf [surfIndex $optic $id]
	for {set j 0} {$j <= 20} {incr j} {
	   set xs [expr $outstopmax * $j/20.]
	   set xs1 $xs

	   if {$xs1 > $outstop1} {set xs1 $outstop1}
	   set zs1 [zsurf $optic $isurf1 $xs1 0.]
	   lappend xx1 $xs
	   lappend zz1 $zs1

	   set xs2 $xs
	   if {$xs2 > $outstop} {set xs2 $outstop}
	   set zs [zsurf $optic $isurf $xs2 0.]
	   lappend xx $xs
	   lappend zz $zs
	   lappend diffz [expr $zs - $zs1]
	   }

#Approximate volume of finished piece.  Works even with tilts.
	set vol 0.
	for {set j 1} {$j <= 20} {incr j} {
	   set j1 [expr $j-1]
	   set x1 [lindex $xx $j1]
	   set x [lindex $xx $j]
	   set z11 [lindex $zz1 $j1]
	   set z12 [lindex $zz1 $j]
	   set z21 [lindex $zz $j1]
	   set z22 [lindex $zz $j]

	   if {$index1 > 0} {
		set dz1 [expr $sep + $z21 - $z11]
		set dz [expr $sep + $z22 - $z12]
	   } else {
		set dz1 [expr $sep - $z21 + $z11]
		set dz [expr $sep - $z22 + $z12]
		}
	   set vol [expr $vol + 3.14*($x1+$x) * ($x*$dz1 - $x1*$dz) + \
		.667*3.14*($x*$x + $x*$x1 + $x1*$x1)*($dz - $dz1)]
	   }
	set zmin1 [eval min $zz1]
	set zmax1 [eval max $zz1]

	set zmin [eval min $zz]
	set zmax [eval max $zz]

	set diffmin [eval min $diffz]
	set diffmax [eval max $diffz]

#Compute the peak-peak thickness of the lens.
#Which is the leading and which is the trailing surface?
	if {$index1 > 0} {

#"height" is the lens height, which gives the thickness of the blank before
#carving.
	   set height [expr $zmax + $sep - $zmin1]

#"thin" is the thinnest cross-section of the lens in the z direction
	   set thin [expr $diffmin + $sep]

#"thick is the thickest cross-section of the lens in the z direction
	   set thick [expr $diffmax + $sep]
	} else {
	   set height [expr $zmax1 + $sep - $zmin]
	   set thin [expr $sep - $diffmax]
	   set thick [expr $sep - $diffmin]
	   }
	set diam [expr 2.*$outstopmax]

#Outer diameter including rim
	set diamrim [expr $diam + 2.*$rim]

#Volume with rim
	set x1 $x
	set x [expr $x1 + $rim]
	set vol [expr $vol + 3.14*$dz*($x*$x - $x1*$x1)]

#Express volume in cc.  Weight in lbs.
	set vol [expr $vol/1.e3]

#Weight in lbs:
	set density [density $glass1]
	set weight [expr 2.2*$vol*$density/1.e3]

#Weight of blank in kg (need to convert mm^3 to cm^3 and g to kg):
	set blank [expr (3.1416*$diamrim*$diamrim/4.)*$height*$density/1.e6]

	set format2 "%3.0f/%-3.0f %7.1f %7.1f %6.1f/%-6.1f %7.1f \
	   %7.0f %6.0f %7.0f %8s"
	puts stdout [format $format2 $id1 $id $diam $height \
	   $thin $thick $diamrim $vol $weight $blank $glass1]
	}
   return
   }

###########################################################################
#Quick routine to get "sag" of a surface.  Ignores tilts.

proc surfSag {optic surfId} {
   set xrad [showSurf $optic $surfId outstop]
   set zcen [zsurf $optic $surfId 0. 0.]
   set zrad [zsurf $optic $surfId $xrad 0.]

#Sag to nearest micron
   set sag [format %.3f [expr abs($zrad -$zcen)]]
   return $sag
   }

###########################################################################
#Compute transmission of lenses.  There is no unique answer here.
#Run a set of rays for a location at the edge of the field.

proc lensTrans {optic filter} {

#Trace chief ray and get list of surfaces
   set xmm [showFocal $optic $filter xrad]
   set ymm 0
   rayPupil $optic 0 0 0 0 $filter 0
   set surfids [surfIdsGet $optic $filter]

   set transtot 1.
   set format1 "  Surfaces  %-10s %10s %10s"
   puts stdout [format $format1 Glass Thick(mm) Transmission]
   set totpath 0

#Link surfaces to construct lenses.
   indexSet $optic $filter
   loop k 0 [llength $surfids] {
	if {$k == 0} continue
	set k1 [expr $k-1]
	set name1 [lindex $surfids $k1]
	set name [lindex $surfids $k]

	set index1 [showIntIndex $optic $name1]
	set index [showIntIndex $optic $name]

#If we reverse direction, don't call it a lens.
	if {$index1*$index <= 0} continue

#If first index is not air or vacuum, continue
	if {abs($index1) <= 1} continue

#Store shape of surface in lens-o-centric coords.
	set xx1 ""
	set zz1 ""
	set xx ""
	set zz ""
	set diffz ""
	set path 0.
	set n 0
	foreach pair {{-1 0} {1 0} {0 -1} {0 0} {0 1}} {
	   set i [lindex $pair 0]
	   set j [lindex $pair 1]
		rayPupil $optic $xmm 0 $j $i $filter 0
		set x1 [exprGet $optic.diagram->xray<$k>]
		set y1 [exprGet $optic.diagram->yray<$k>]
		set z1 [exprGet $optic.diagram->zray<$k>]
#??? I don't want this.
#		set k1 [expr $k+1]
		set x [exprGet $optic.diagram->xray<$k1>]
		set y [exprGet $optic.diagram->yray<$k1>]
		set z [exprGet $optic.diagram->zray<$k1>]
		set path [expr $path + sqrt(pow($x-$x1,2) + pow($y-$y1,2) + \
		   pow($z-$z1,2))]
		incr n
	   }
	set glass [showGlass $optic $name1]
	set path [expr $path/$n]
	set totpath [expr $totpath + $path]
	set wave [showFocal $optic $filter wave]

#Transmission - 25 mm reference thickness
	set trans [transmit $glass $wave]

#Correct for thickness
	set trans [expr pow($trans,$path/25.)]
	set format2 "%5.0f/%-5.0f %-10s %10.1f %10.2f"
	puts stdout [format $format2 $name1 $name $glass $path $trans]
	if {$trans > 0.} {set transtot [expr $transtot*$trans]}
	}
   puts stdout [format "Total transmission %.2f" $transtot]
   puts stdout [format "Total path length  %.1f mm" $totpath]
   return
   }

######################################################################
#For an aspheric surface, compute maximum deviation from a sphere
proc asphere {optic surf} {
   set max -1.e10
   set min 1.e10
   set isurf [surfIndex $optic $surf]
   set outstop [showSurf $optic $surf outstop]
   set c [showSurf $optic $surf curv]
   if {$c != 0.} {
	set rtrue [expr 1./$c]
   } else {
	set rtrue 0.
	}
   set ksum 0.
   set weight 0
   set r2list ""
   set zlist ""
   set csum 0
   set cweight 0
   for {set j -50} {$j <= 50} {incr j} {
	set xmm [expr $outstop*$j/50.]
	set ymm 0
	set z [zsurf $optic $isurf $xmm $ymm]
	set r2 [expr $xmm*$xmm + $ymm*$ymm]
	set r $xmm
	set factor [expr 1.-($c)*($c)*$r2]
	if {$factor < 0.} {
	   set zc [expr 1./$c]
	} else {
	   set zc [expr $c*$r2/(1.+sqrt($factor))]
	   }
	set max [max $max [expr $z-$zc]]
	set min [min $min [expr $z-$zc]]

#Estimate of conic constant
	if {$z != 0 && $c != 0} {
	   set k [expr (-$r2*$c + 2.*$z - $c*$z*$z)/($c*$z*$z)]
	   set ksum [expr $ksum + $r2*$k]
	   set weight [expr $weight + $r2]
	   }

#Estimate best-fitting curvature
	if {$j == 0} {
	   set z0 $z
	} elseif {$j > 0} {
	   set dz [expr $z - $z0]
	   set cbest [expr 2.*$dz/($r2+pow($dz,2))]
	   set csum [expr $csum + $cbest]
	   set cweight [expr $cweight + 1.]
	   }
	lappend rlist $r
	lappend r2list $r2
	lappend zlist $z
	}
   echo Peak deviations from a sphere: max [format %.2f $max] min \
	[format %.2f $min]

   if {$weight > 0.} {
	set k [expr $ksum/$weight]
	echo
	echo Approximate conic constant: [format %.5f $k]
	set max -1.e10
	set min 1.e10
	loop j 0 [llength $r2list] {
	   set z [lindex $zlist $j]
	   set r2 [lindex $r2list $j]
	   set factor [expr 1.-($c)*($c)*$r2*(1.+$k)]
	   set zc [expr $c*$r2/(1.+sqrt($factor))]
	   set max [max $max [expr $z-$zc]]
	   set min [min $min [expr $z-$zc]]
	   }
	echo Peak deviations from asphere: max [format %.2f $max] min \
	   [format %.2f $min]
	}

   if {$cweight > 0.} {
	set cbest [expr $csum/$cweight]
	if {$cbest != 0.} {
	   set rbest [expr 1./$cbest]
	} else {
	   set rbest 0.
	   }
	set c $cbest
	echo Initial rbest [format %.1f $rbest]

#Get better estimate of rbest.  Do least squares fit!
	loop iter 0 2 {
	   set nparam 2
	   set l 0
	   loop j 0 $nparam {
		set rhs($j) 0.
		for {set k 0} {$k <= $j} {incr k} {
		   set cmat($l) 0.
		   incr l
		   }
		}

	   set nobs [llength $r2list]
	   loop i 0 $nobs {
		set r2 [lindex $r2list $i]
		set z [lindex $zlist $i]

#Predicted z
		set zcomp [expr $c*$r2/(1. + sqrt(1.-$c*$c*$r2))]
		set res [expr $z - $zcomp]
		set deriv(0) 1.
		set deriv(1) [expr $r2/(1. + sqrt(1.-$c*$c*$r2)) * (1. + \
		   $c*$zcomp/sqrt(1.-$c*$c*$r2))]
		set l 0
		loop j 0 $nparam {
		   set rhs($j) [expr $rhs($j) + $res*$deriv($j)]
		   for {set k 0} {$k <= $j} {incr k} {
			set cmat($l) [expr $cmat($l) + $deriv($j) * $deriv($k)]
			incr l
			}
		   }
		}
	   syminv cmat rhs
	   set offset $rhs(0)
	   set c [expr $c + $rhs(1)]
	   if {$c != 0.} {
		set rbest [expr 1./$c]
	   } else {
		set rbest 0.
		}
	   }
	echo Actual radius of curvature: [format %.1f $rtrue]
#	echo Best fitting radius of curvature: [format %.1f $rbest]  Offset \
#	   [format %.1f $offset]
	set c [expr 1./$rbest]

#Forget all the above - just find rbest so deviation at edge is 0.
	set z [lindex $zlist end]
	set r2 [lindex $r2list end]
	set c [expr 2.*$z/($r2 + $z*$z)]
	set rbest [expr 1./$c]
	echo Best fitting radius of curvature: [format %.1f $rbest]
	set max -1.e10
	set min 1.e10
	set slopemax 0.
	set curvmax 0.
	set totcurvmax 0.
	set r [lindex $rlist 0]
	loop j 0 [llength $r2list] {
	   set rold $r
	   set z [lindex $zlist $j]
	   set r2 [lindex $r2list $j]
	   set r [lindex $rlist $j]
	   set k 0.
	   set factor [expr 1.-($c)*($c)*$r2*(1.+$k)]
	   set zc [expr $c*$r2/(1.+sqrt($factor))]

#Cache height differentials so we can calculate slopes and curvatures.
	   set zdiff($j) [expr $z-$zc]
	   set zcache($j) $z
	   set max [max $max $zdiff($j)]
	   set min [min $min $zdiff($j)]
	   if {$j == 0} continue
	   set slope($j) [expr ($zdiff($j) - $zdiff([expr $j-1])) / ($r-$rold)]

#More accurate - take cross-product.  However, this is really not necessary.
	   set slopemax [expr max($slopemax, abs($slope($j)))]
	   if {$j <= 1} continue

#Curvature of height differentials.  Assume we have uniform spacing in radius.
	   set r0 [lindex $rlist [expr $j-2]]
	   set r1 [lindex $rlist [expr $j-1]]
	   set r2 [lindex $rlist $j]
	   set z0 $zdiff([expr $j-2])
	   set z1 $zdiff([expr $j-1])
 	   set z2 $zdiff($j)
	   set curv [curv3p $r0 $z0 $r1 $z1 $r2 $z2]
	   set curvmax [expr max($curvmax,abs($curv))]

#For fun, compute total curvature as well.  This may be more important?
	   set r0 [lindex $rlist [expr $j-2]]
	   set r1 [lindex $rlist [expr $j-1]]
	   set r2 [lindex $rlist $j]
	   set z0 [lindex $zlist [expr $j-2]]
	   set z1 [lindex $zlist [expr $j-1]]
 	   set z2 [lindex $zlist $j]
	   set totcurv [curv3p $r0 $z0 $r1 $z1 $r2 $z2]
#echo totcurv $totcurv z $zcache($j)
	   set totcurvmax [expr max($totcurvmax,abs($totcurv))]
	   }
	echo Peak deviations from best-fit sphere: max [format %.2f $max] min \
	   [format %.2f $min]
	set slopemax [expr $slopemax*1.e3]
	if {$curvmax != 0.} {
	   set radmax [expr 1./$curvmax]
	} else {
	   set radmax 0.
	   }
	if {$totcurvmax != 0.} {
	   set totradmax [expr 1./$totcurvmax]
	} else {
	   set totradmax 0.
	   }
	echo Peak slope w.r.t. sphere is [format %.2f $slopemax] microns/mm.
	echo Small radius of curvature w.r.t. sphere is \
	   [format %.1f $radmax] mm.
	echo Small total radius of curvature is [format %.1f $totradmax] mm.

#Make a plot
	plotInit a
	pgEnv 0 [lindex $rlist end] $min $max 0 0
	pgLabel "Radius" "Sag w.r.t. best fitting sphere" "Surface $surf"
	loop j 0 [llength $rlist] {
	   set r [lindex $rlist $j]
	   pgPoint $r $zdiff($j) 3
	   }
	plotInit b
	pgEnv 0 [lindex $rlist end] -50. 50. 0 0
	pgLabel "Radius" "Slope (microns/mm)" "Surface $surf"
	loop j 1 [llength $rlist] {
	   set r [lindex $rlist $j]
	   set s [expr $slope($j)*1.e3]
	   pgPoint $r $s 3
	   }
	}

   return
   }

######################################################################
#Plot chief rays at many field points.

proc chiefPlot {optic filter} {
   set PI 3.141593

#Size of focal plane
   set xmm [showFocal $optic $filter xsize]
   set ymm [showFocal $optic $filter ysize]
   if {$xmm < 0 || $ymm < 0} {
	set rad [expr sqrt(pow($xmm,2) + pow($ymm,2))]
   } else {
	set rad $xmm
	}

#Finner is a useful surrogate for knowing the inner stop of the telescope.
#Normally I don't set this in the aperture stop surface, since it is set
#by the secondary mirror.

   set finner [exprGet (*$optic.tel).finner]

#Place spots at center of field and at max. x position.  Keep y=0.
   set pgSci 1
   foreach ixspot "0 .2 .4 .6 .8 1." {
	incr pgSci
	if {$pgSci == 8} {set pgSci 2}
	pgSci $pgSci
	set x [expr $ixspot*$rad]
	set y 0
	foreach frad "0" {
	   rayPupil $optic $x $y $frad 0 $filter 0
	   set np [exprGet (*$optic.diagram).np]
	   set xx ""
	   set zz ""
	   loop i 0 $np {
		lappend xx [exprGet (*$optic.diagram).xray<$i>]
		lappend zz [exprGet (*$optic.diagram).zray<$i>]
		}
	   pgLine $xx $zz
	   }
	}
   pgSci 1
   return
   }

#####################################################################
#Top level routine to plot chief rays for ghosting.
#I input first and last surfaces to plot and ghosting surface.
#Ghosting is between surfaces ghost1 and ghost2.

proc ghostPlot {hndl surf1 surf2 ghost1 ghost2 {ifil 1}} {
   set hndl2 [ghostDup $hndl $ghost1 $ghost2 $ifil]
   opticPlot $hndl $surf1 $surf2 $ifil no
   chiefPlot $hndl2 $ifil
   opticDel $hndl2
   return
   }

#####################################################################
#Helper function - compute curvature, given 3 points in a plane
#Horizontal or vertical lines cause crash - beware!

proc curv3p {x0 y0 x1 y1 x2 y2} {

#Compute intersection of bisectors of line segments to get center of
#circle, then compute radius straightforwardly.

   set xc1 [expr ($x1 + $x0)/2.]
   set yc1 [expr ($y1 + $y0)/2.]
   set m1 [expr ($y1-$y0)/($x1-$x0)]
   set a1 [expr -1./$m1]
   set b1 [expr $yc1 + $xc1/$m1]

   set xc2 [expr ($x2 + $x1)/2.]
   set yc2 [expr ($y2 + $y1)/2.]
   set m2 [expr ($y2-$y1)/($x2-$x1)]
   set a2 [expr -1./$m2]
   set b2 [expr $yc2 + $xc2/$m2]

   set xc [expr ($b2-$b1)/($a1-$a2)]
   set yc [expr $a1*$xc + $b1]
   set rad [expr sqrt(pow($x0-$xc,2) + pow($y0-$yc,2))]
   set curv [expr 1./$rad]
   return $curv
   }

########################################################################
#
#Plot the shape of one surface.  Useful when modeling surfaces as zerikes,
#e.g.  Scale is NOT maintained

proc zshape {hndl fsurf} {
   set isurf [surfindex $hndl $fsurf]
   set rad [expr abs([showSurf $hndl $fsurf outstop])]
   set xlist ""
   set zlist ""
   set zmin 1.e10
   set zmax -1.e10
   loop i -20 21 {
	set x [expr $rad*($i/20.)]
	set z [zsurf $hndl $isurf $x 0.]
	lappend xlist $x
	lappend zlist $z
	set zmin [min $zmin $z]
	set zmax [max $zmax $z]
	}
   set zdiff [expr $zmax-$zmin]
   if {$zdiff == 0} return
   set sag [format %.3f [expr abs($zdiff)]]
   plotInit b
   set zmax [expr $zmax + .1*$zdiff]
   set zmin [expr $zmin - .1*$zdiff]
   pgEnv -$rad $rad $zmin $zmax 0 0
   pgLine $xlist $zlist
   pgLabel "Radius" "Z Height" "Surface shape for surface $fsurf"
   echo Sag is $sag mm
   return
   }

