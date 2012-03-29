##########################################################################
#Read an atmos file.
#First line is a header
#Next line gives nsurf, nwave, other stuff
#Next line has wavelengths (in microns)
#Each surface is a group of 6 lines.

#This file probably does not work for mirrors anymore - I changed my
#conventions for storing refractive indices to not require negative
#values anymore, but atmos has always used my original convention!
#Not yet fixed.

proc atmRead {file} {
   set ext [file extension $file]

#Default extension is .atm
   if {$ext == ""} {set file $file.atm}
   if {![file exists $file]} {
	error "No file found $file"
	}
   set fid [open $file]

#First line has titles.
   set line [gets $fid]
   set optic [opticNew]

#Next line has countinfo
   set line [split [gets $fid] ,]
   set nsurf [lindex $line 0]
   set nwave [lindex $line 1]

#Next line has wavelengths
   set line [split [gets $fid] ,]
   loop i 0 $nwave {
	set wave [lindex $line $i]
	set indx [expr $i+1]
	set waveList($indx) $wave
	}

#Focal plane - fill in fixed quantities.
   foreach ifil [array names waveList] {
	setFocal $optic $ifil xpos 0
	setFocal $optic $ifil ypos 0
	setFocal $optic $ifil wave $waveList($ifil)
	setFocal $optic $ifil weight 1
	setIndex $optic 0 $ifil 1
	}
   handleSet $optic.name "$file"
   handleSet $optic.ncolor [array size waveList]
   setSurf $optic 0 z -1.e14

   set z 0
   set surf 0
   set sign 1.
   set astop 0
   loop i 0 $nsurf {
	incr surf
	set line [split [gets $fid] ,]
	set stop [lindex $line 0]
	set rad [lindex $line 1]
	set thick [lindex $line 2]

#Meaning of stop depends on surface type, which we don't have yet.
	set line [split [gets $fid] ,]

#gtype is actually the glass catalog.
#It is conjectured that 0 = air, 1 = air after a mirror, 2 = Schott,
#3 = Ohara, 4 = Corning, 5 = Hoya, 6 = Direct
	set gtype [lindex $line 0]

#gtype 3 is OHara.  For some reason, glass name is 3rd item in list
	set gindex 1
	if {$gtype == 3} {set gindex 2}
	set glass [lindex [lindex $line $gindex] 0]

	set line [split [gets $fid] ,]
	loop j 0 $nwave {
	   set jndx [expr $j+1]
	   set index($jndx) [lindex $line $j]
	   }

#Conic is on next line.
	set line [split [gets $fid] ,]
	set surftype [lindex $line 0]
	set param [lindex $line 1]

#Next line has asphere stuff
	set line [split [gets $fid] ,]
	set a2 [lindex $line 0]
	set a4 [lindex $line 1]
	set a6 [lindex $line 2]
	set a8 [lindex $line 3]

#Last  line is of unknown purpose
	set line [split [gets $fid] ,]

#gtype actually gives the catalog source or specified mirror or air,
#so I could redo the following logic:
	if {$gtype != 2 && $gtype != 3} {
	   if {$index(1) == 1} {
		set glass air
	   } elseif {$index(1) == -1} {
		set glass mirror
	   } else {
		set glass unknown
		}
	   }
	if {$rad != 0.} {
	   set curv [expr 1./$rad]
	} else {
	   set curv 0.
	   }

#Revert back some of my conventions
	if {$glass == "Fused-Silica"} {set glass SIO2}
	if {$glass == "reflect"} {set glass mirror}

	set name S$surf
	setSurf $optic $surf z $z
	setSurf $optic $surf curv $curv
	foreach ifil [array names waveList] {

#Sign is correct for old scheme - not for new (?!)
	   set index($ifil) [expr $index($ifil)*$sign]
	   setIndex $optic $surf $ifil $index($ifil)
	   }
	if {$glass == "mirror"} {
	   set sign [expr -1.*$sign]
	   }
	setGlass $optic $surf $glass
	setName $optic $surf $name
	set z [expr $z + $thick]
	if {$surftype == 4} {
	   setSurf $optic $surf stoptype 1
	   setSurf $optic $surf instop $stop
	} else {
	   setSurf $optic $surf stoptype 0
	   setSurf $optic $surf outstop $stop
	   }

#Surftype means what?  If 2, conic, if 3, conic + asphere?
	if {$surftype == 2 || $surftype == 3} {
	   setSurf $optic $surf ccon $param
	   }

	foreach p "a2 a4 a6 a8" {
	   setSurf $optic $surf $p [set $p]
	   }

#Check if we should use this as the aperture stop
#Here I am really guessing.
	if {$surftype == 1 && $astop == 0} {
	   setSurf $optic $surf stoptype 2
	   set astop 1
#I will set z of my aperture stop to 0 later.
	   set zapp [showSurf $optic $surf z]
	   }
	}

#Last stop value is istop
   set istop $stop
   if {![info exist istop]} {set istop 10.}

#Reposition z
   for {set i 1} {$i <= [exprGet $optic.nsurf]} {incr i} {
	set surfid [surfId $optic $i]
	setSurf $optic $surfid z [expr [showSurf $optic $surfid z] - $zapp]
	}

#Remaining focal plane parameters
   foreach ifil [array names waveList] {
	setFocal $optic $ifil xrad $istop
	setFocal $optic $ifil yrad $istop
	}

   rayPattern $optic 6 1
   colorcount $optic

#Scale factor.  Use a single value for all filters.
   set sum 0.
   set n 0
   foreach ifil [array names waveList] {
	opticInfo $optic $ifil
	set fl [showFocal $optic $ifil fl]
	set scale [expr 206265./$fl]
	set sum [expr $sum + $scale]
	incr n
	}

#Don't multiply by sign, since it is not always correct.
   set scale [format %.3f [expr $sum/$n]]

   foreach ifil [array names waveList] {
	setFocal $optic $ifil scale $scale
	}
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }
