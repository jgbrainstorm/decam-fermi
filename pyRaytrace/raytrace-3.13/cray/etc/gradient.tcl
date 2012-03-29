#Support for index gradients in glass.
#Gradients are intended for approximate analysis only [eg., impact of
#temp. gradients].

#Approach: In C code, add support for a radial gradient in the refractive
#index.  In TCL, add support to take a single lens and split it into
#multiple lenses with the same overall shape.  This will allow me to
#implement, crudely, a 2-dimensional gradient pattern.

#Split a lens into n pieces.
proc lensSplit {hndl surfid1 ifil n} {

#Get fsurf2 - the surface after fsurf1
   set ids [surfIdsGet $hndl $ifil]
   set indx [lsearch $ids $surfid1]
   if {$indx < 0} {
	error "Bad first surface $surfid1; not in design"
	}
   if {$indx == [llength $ids]-1} {
	error "First surface $surfid1 is last surface"
	}
   set surfid2 [lindex $ids [expr $indx+1]]
   echo Surface are $surfid1 $surfid2
   set x1 [showSurf $hndl $surfid1 x]
   set x2 [showSurf $hndl $surfid2 x]
   set y1 [showSurf $hndl $surfid1 y]
   set y2 [showSurf $hndl $surfid2 y]
   set z1 [showSurf $hndl $surfid1 z]
   set z2 [showSurf $hndl $surfid2 z]
   set theta1 [showSurf $hndl $surfid1 theta]
   set theta2 [showSurf $hndl $surfid2 theta]
   set curv1 [showSurf $hndl $surfid1 curv]
   set curv2 [showSurf $hndl $surfid2 curv]
   set glass [showGlass $hndl $surfid1]
   set name [showName $hndl $surfid1]
   set newid $surfid1
   loop i 1 $n {

#I could interpolate lots of other parameters here as well - e.g., asphere
#terms.
	set x [expr $x1 + ($x2-$x1)*($i/(1.*$n))]
	set y [expr $y1 + ($y2-$y1)*($i/(1.*$n))]
	set z [expr $z1 + ($z2-$z1)*($i/(1.*$n))]
	set theta [expr $theta1 + ($theta2-$theta1)*($i/(1.*$n))]
	set newid [surfInsert $hndl $newid $z $glass]
	set curv [expr $curv1 + ($curv2-$curv1)*($i/(1.*$n))]
	setSurf $hndl $newid x $x
	setSurf $hndl $newid y $y
	setSurf $hndl $newid theta $theta
	setSurf $hndl $newid curv $curv
	setName $hndl $newid $name
	}
   stopComp $hndl
   return
   }

######################################################################
#Increment a refractive index by a factor that is supposed to represent
#a change in temperature.

proc indexInc {hndl surfid ifil fact} {
   set index [showIndex $hndl $surfid $ifil]
   set index [expr $index + $fact]
   setIndex $hndl $surfid $ifil $index
   return
   }
