#Let's modify spectro2u for different focal lengths
#Coll is collimator focal length
#cam is camera focal length

proc modify {h coll cam} {

#Hardwire focal ratio
   set fr 2.8

#old scale
   set oldscale [showFocal $h 2 scale]
   set ncolor [exprGet $h.ncolor]
   set fsurf [exprGet $h.nsurf]

#Front end
   setSurf $h 0 z $coll
   setSurf $h 1 outstop [expr $coll/(2.*$fr)]
   setSurf $h 2 z [expr 2.*$coll]
   setSurf $h 2 curv [expr -1./(2.*$coll)]

#Back end - scale z, curv
   set scale [expr -206265./$cam]
   set fact [expr $scale/$oldscale]
echo Fact = $fact
   foreach surf [range 6-$fsurf] {
	set z [showSurf $h $surf z]
	set curv [showSurf $h $surf curv]
	set z [expr $z/$fact]
	set curv [expr $curv*$fact]
	setSurf $h $surf z $z
	setSurf $h $surf curv $curv
	}
   foreach ifil [range 1-$ncolor] {
	setFocal $h $ifil scale $scale
	}
   return
   }
