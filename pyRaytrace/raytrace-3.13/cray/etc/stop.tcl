#Code to add stops to a design

#First, add an arbitrary stop after surface isurf; specify z position,
#inner radius (mm), outer radius (mm)

proc stopAdd {optic fsurf z inner outer} {
   set surfId [opticInsert $optic $fsurf]
   loop i 0 [exprGet $optic.ncolor] {
	set index [showIndex $optic $fsurf [expr $i+1]]
	setIndex $optic $surfId [expr $i+1] $index
	}
   set glass [showGlass $optic $fsurf]
   setGlass $optic $surfId $glass
   setSurf $optic $surfId z $z
   setSurf $optic $surfId instop $inner
   setSurf $optic $surfId outstop $outer
   setSurf $optic $surfId stoptype 1
   setName $optic $surfId STOP
   return
   }

#################################################################
#Rotate some points counterclockwise
proc stopRotate {x y ang} {
   set rang [expr $ang*3.141593/180.]
   set x1 [expr $x*cos($rang) - $y*sin($rang)]
   set y1 [expr $y*cos($rang) + $x*sin($rang)]
   return [list $x1 $y1]
   }

#################################################################
#Next, add a strut type stop.  Specify previous surface, z position,
#position angle, and width (in mm).  Assume that it goes the entire radius
#of the aperture.

proc strutAdd {optic fsurf z width ang} {
   set surfId [opticInsert $optic $fsurf]
   set halfwidth [expr $width/2.]
   loop i 0 [exprGet $optic.ncolor] {
	set index [showIndex $optic $fsurf [expr $i+1]]
	setIndex $optic $surfId [expr $i+1] $index
	}
   set glass [showGlass $optic $fsurf]
   setGlass $optic $surfId $glass
   setSurf $optic $surfId z $z

#Get radius of aperture
   set radius [expr [telDiam $optic]/2.]
   set x 0
   set y $halfwidth
   set list [stopRotate $x $y $ang]
   setSurf $optic $surfId x1 [lindex $list 0]
   setSurf $optic $surfId y1 [lindex $list 1]
   set x $radius
   set y $halfwidth
   set list [stopRotate $x $y $ang]
   setSurf $optic $surfId x2 [lindex $list 0]
   setSurf $optic $surfId y2 [lindex $list 1]
   set x $radius
   set y [expr -1.*$halfwidth]
   set list [stopRotate $x $y $ang]
   setSurf $optic $surfId x3 [lindex $list 0]
   setSurf $optic $surfId y3 [lindex $list 1]
   set x 0
   set y [expr -1.*$halfwidth]
   set list [stopRotate $x $y $ang]
   setSurf $optic $surfId x4 [lindex $list 0]
   setSurf $optic $surfId y4 [lindex $list 1]
   setSurf $optic $surfId stoptype 3
   setName $optic $surfId STRUT
   return
   }

#######################################################################
#Make a face-on plot of stops.

proc stopPlot {optic} {
   plotInit b
   pgSci 1
   set xmax 0
   set nsurf [exprGet $optic.nsurf]
   loop i 0 $nsurf {
	set id [surfId $optic $i]
	set type [showSurf $optic $id stoptype]
	set outstop [showSurf $optic $id outstop]
	if {$type == 2} {set xmax $outstop}
	}
   if {$xmax == 0} {
	set xmax [expr [telDiam $optics]/2.]
	}
   if {$xmax == 0} return
   pgEnv -$xmax $xmax -$xmax $xmax 0 0
   pgSfs 2
   loop i 0 $nsurf {
	set id [surfId $optic $i]
	set type [showSurf $optic $id stoptype]
	set outstop [showSurf $optic $id outstop]
	if {$type == 0} continue
	if {$type == 1} {
	   set instop [showSurf $optic $id instop]
	   set outstop [showSurf $optic $id outstop]
	   set x [showSurf $optic $id xoff]
	   set y [showSurf $optic $id yoff]
	   pgSci 7
	   pgSfs 1
	   pgCirc $x $y $instop
	   pgSfs 2
	   pgCirc $x $y $outstop
	   }
	if {$type == 2} {
	   set instop [showSurf $optic $id instop]
	   set outstop [showSurf $optic $id outstop]
	   set x [showSurf $optic $id xoff]
	   set y [showSurf $optic $id yoff]
	   pgSfs 2
	   pgCirc $x $y $instop
	   pgCirc $x $y $outstop
	   }
	if {$type == 3} {
	   set x1 [showSurf $optic $id x1]
	   set y1 [showSurf $optic $id y1]
	   set x2 [showSurf $optic $id x2]
	   set y2 [showSurf $optic $id y2]
	   set x3 [showSurf $optic $id x3]
	   set y3 [showSurf $optic $id y3]
	   set x4 [showSurf $optic $id x4]
	   set y4 [showSurf $optic $id y4]
	   pgSfs 1
	   pgSci 6
	   pgPoly "$x1 $x2 $x3 $x4 $x1" "$y1 $y2 $y3 $y4 $y1"
	   }
	}
   return
   }

##########################################################################
#Get the mirror diameter.  Fetch from aperture stop.  This is done to
#help get rid of tel structure.

#Note: I should really input a filter number, fetch all surfaces for this
#filter, then get aperture stop just for that filter.  This is necessary
#in case the stop is different for different filters.

#I now do this.  Default filter is 1.

proc telDiam {optic {ifil 1}} {
   set diam 0
   set surfids [surfIdsGet $optic $ifil]

   foreach surfid $surfids {
       set type [showSurf $optic $surfid stoptype]
        if {$type == 2} {
           set rad [showSurf $optic $surfid outstop]
           set diam [expr abs($rad)*2.]
           break
           }
        }
   if {$diam == 0} {error "No nonzero aperture stop found"}
   return $diam
   }

########################################################################
#For surfaces that are "linked", the outer stops that are computed by
#stopcomp should be linked as well, but are not.  If I supply no arguments,
#I will look for all surfaces that have same integer ID and set outstop to
#the biggest.  If I supply a list of surfaces, I will only process those
#surfaces.
#
#In some designs, linked surfaces differ in the integer identifier.  Those
#will have to be dealt with separately - can't do everything here.
#Also, the only design where I actually did that is the SDSS design.

proc stopLink {hndl {surfs ""}} {
   set nsurf [exprGet $hndl.nsurf]
   if {$surfs == ""} {
	for {set i 1} {$i <= $nsurf} {incr i} {
	   set surfId [surfId $hndl $i]
	   set idint [expr int($surfId)]
	   if {[lsearch $surfs $idint] < 0} {
		lappend surfs $idint
		}
	   }
	}
   foreach surf $surfs {
	set id($surf) ""
	set outstop($surf) ""
	}
   for {set i 1} {$i <= $nsurf} {incr i} {
	set surfId [surfId $hndl $i]
	set idint [expr int($surfId)]
	if {[lsearch $surfs $idint] < 0} continue
	lappend id($idint) $surfId
	lappend outstop($idint) [expr abs([showSurf $hndl $surfId outstop])]
	}
   foreach surf $surfs {
	if {[llength $id($surf)] <= 1} continue
	set stop [eval max $outstop($surf)]
	foreach surfId $id($surf) {
	   setSurf $hndl $surfId outstop $stop
	   }
	}
   return
   }
