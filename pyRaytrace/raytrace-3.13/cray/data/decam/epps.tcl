#Create an optical design from harland Epps corrector design.

proc epps {} {
   set fid [open epps.txt]
   set line [gets $fid]
   set optic [opticNew]

#Filter wavelengths
   set waveList(1) .39
   set waveList(2) .49
   set waveList(3) .51
   set waveList(4) .59
   set waveList(5) .59
   set waveList(6) .81
   set waveList(7) .76
   set waveList(8) 1.

#Focal plane
   foreach ifil [array names waveList] {
	setFocal $optic $ifil xpos 0
	setFocal $optic $ifil ypos 0
	setFocal $optic $ifil xrad 154.3
	setFocal $optic $ifil yrad 154.3
	setFocal $optic $ifil wave $waveList($ifil)
	setIndex $optic 0 $ifil 1
	}
   handleSet $optic.name "Epps Design"
   handleSet $optic.ncolor [array size waveList]
   setSurf $optic 0 z -1.e14

   set z 0

#Insert stop
   stopAdd $optic 0 -8442. 700 4000
   set isurf 1
   setName $optic $isurf STOP
   while {1} {
	set line [gets $fid]
echo $line
	if {[eof $fid]} break
	set surf [lindex $line 0]
	set rad [lindex $line 1]
	set thick [lindex $line 2]
	set glass [lindex $line 3]

#If surface number is repeated, this line has aspheric info
	if {[info exists lastsurf] && $surf == $lastsurf} {
	   set param $rad
	   setSurf $optic $isurf $param $thick
	   continue
	   }

#I've let surface 1 be the stop, so increment numbers here.
	incr isurf
	set lastsurf $surf
	if {$rad != 0.} {set curv [expr 1./$rad]} else {set curv 0}

#Primary - I will create 4 surfaces, one for each filter, and allow refocus
#for each filter.
	if {$isurf == 2} {
	   foreach filter "1 2 3 4" {
		set surfid $isurf.[format %02d $filter]
		setSurf $optic $surfid z $z
		setSurf $optic $surfid curv $curv
#Parameters not in file
		setSurf $optic $surfid stoptype 2
		setSurf $optic $surfid outstop 2000.
#Paraboloid on primary
		setSurf $optic $surfid ccon -1.0
		set ifil1 [expr $filter*2-1]
		set ifil2 [expr $ifil1+1]
		foreach ifil "$ifil1 $ifil2" {
		   set index [glass $glass $waveList($ifil)]
		   set index [expr -1.*$index]
		   setIndex $optic $surfid $ifil $index
		   }
		setGlass $optic $surfid $glass
		setName $optic $surfid PRIMARY
		}
	} else {
	   setSurf $optic $isurf z $z
	   setGlass $optic $isurf $glass
	   setName $optic $isurf $isurf
	   setSurf $optic $isurf curv $curv
	   foreach ifil [array names waveList] {
		set index [glass $glass $waveList($ifil)]
		set index [expr -1.*$index]
		setIndex $optic $isurf $ifil $index
		}
	   }
	set z [expr $z + $thick]
	}
   close $fid
   incr isurf
   setSurf $optic $isurf z $z
   foreach ifil [array names waveList] {
	setIndex $optic $isurf $ifil -1
	}

   setName $optic $isurf FOCAL
   rayPattern $optic 6 1
   colorcount $optic

#Scale factor
   set scale 0
   set n 0
   foreach ifil [array names waveList] {
	opticInfo $optic $ifil
	set fl [showFocal $optic $ifil fl]
	set scale [expr $scale + 206265./$fl]
	incr n
	}
   set scale [expr -1.*$scale/$n]

#Hardwire - paraxial calculation is not very good.
   set scale 23.59
   foreach ifil [array names waveList] {
	setFocal $optic $ifil scale $scale
	}
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }
