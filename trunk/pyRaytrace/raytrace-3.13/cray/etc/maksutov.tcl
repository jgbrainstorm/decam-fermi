#Enter information for a Maksutov design and convert to cray format.
#Input format is inspired by http://www.atmsite.org/contrib/Fejes/
#	Optimization1/Optimization1.html

proc maksutov {} {
   puts stdout "Enter parameters for a Maksutov Cassegrain design."
   while {1} {
	puts stdout "Primary mirror diameter: " nonewline; flush stdout
	set D [gets stdin]
	if {![catch {format %f $D}]} break
	puts stdout "Illegal number - try again"
	}
   while {1} {
	puts stdout "Final focal ratio: " nonewline; flush stdout
	set f2 [gets stdin]
	if {![catch {format %f $f2}]} break
	puts stdout "Illegal number - try again"
	}
   while {1} {
	puts stdout "Number of wavelengths: " nonewline; flush stdout
	set nwave [gets stdin]
	if {![catch {format %f $nwave}]} break
	puts stdout "Illegal number - try again"
	}
   set waves ""
   loop i 0 $nwave {
	while {1} {
	   puts stdout "Design wavelength [expr $i+1]: " nonewline
	   flush stdout
	   set wave [gets stdin]
	   if {![catch {format %f $wave}]} break
	   puts stdout "Illegal number - try again"
	   }
	lappend waves $wave
	}
   while {1} {
	puts stdout "Front of corrector, Radius, Thickness: " nonewline
	flush stdout
	set list [gets stdin]
	set R1 [lindex $list 0]
	set dz1 [lindex $list 1]
	if {![catch {format %f $R1}] && ![catch {format %f $dz1}]} break
	puts stdout "Illegal number - try again"
	}
   while {1} {
	puts stdout "Back of corrector, Radius, Thickness: " nonewline
	flush stdout
	set list [gets stdin]
	set R2 [lindex $list 0]
	set dz2 [lindex $list 1]
	if {![catch {format %f $R2}] && ![catch {format %f $dz2}]} break
	puts stdout "Illegal number - try again"
	}
   while {1} {
	puts stdout "Primary mirror, Radius, Thickness: " nonewline
	flush stdout
	set list [gets stdin]
	set R3 [lindex $list 0]
	set dz3 [lindex $list 1]
	if {![catch {format %f $R3}] && ![catch {format %f $dz3}]} break
	puts stdout "Illegal number - try again"
	}
   while {1} {
	puts stdout "Secondary mirror, Radius, Thickness: " nonewline
	flush stdout
	set list [gets stdin]
	set R4 [lindex $list 0]
	set dz4 [lindex $list 1]
	if {![catch {format %f $R4}] && ![catch {format %f $dz4}]} break
	puts stdout "Illegal number - try again"
	}

   while {1} {
	puts stdout "Focal plane, Radius: " nonewline
	flush stdout
	set list [gets stdin]
	set R5 [lindex $list 0]
	if {![catch {format %f $R5}]} break
	puts stdout "Illegal number - try again"
	}

#Fill in only essential telescope params.

#F1 is primary mirror focal length
   set F1 [expr $R3 / 2.]
   set f1 [expr $F1 / $D]
   set back [expr $dz3 + $dz4]

#F3 is overall focal length
   set F3 [expr $f2*$D]
   set beta [expr $back/$F3]
   set M [expr $F3/$F1]

#Now create internal data structures.
   set optic [opticNew]
   handleSet $optic.tel->diam $D
   handleSet $optic.tel->fr1 $f1
   handleSet $optic.tel->fr2 $f2
   handleSet $optic.tel->back $back
   handleSet $optic.tel->fl1 $F1
   handleSet $optic.tel->finner .25
   handleSet $optic.tel->mag $M
   handleSet $optic.tel->beta $beta
   handleSet $optic.tel->sep $dz3
   handleSet $optic.tel->f3 $F3

#Now initialize surfaces.
#Size of focal plane: Use 60 arcminute diameter
   set size [expr abs($F3*(30./60./57.3))]

#Scale in arcsec/mm
   set scale [expr 206265./$F3]

   loop i 0 $nwave {
	set icolor [expr $i+1]
	setFocal $optic $icolor xpos 0
	setFocal $optic $icolor ypos 0
	setFocal $optic $icolor xrad $size
	setFocal $optic $icolor yrad $size
	setFocal $optic $icolor scale $scale
	setFocal $optic $icolor wave [lindex $waves $i]
	setFocal $optic $icolor dist 0
	setFocal $optic $icolor rot 0
	}

#Now surfaces
#Don't set optic.nsurf - setparam automatically creates new surfaces as
#needed.
   handleSet $optic.name "Maksutov-Cassegrain"
   handleSet $optic.ncolor 1

   setSurf $optic 0 z -1.e14
   loop i 0 $nwave {
	setIndex $optic 0 [expr $i+1] 1
	}

#Front of corrector
   set z 0
   setSurf $optic 1 z 0
   setSurf $optic 1 curv [expr 1./$R1]
   loop i 0 $nwave {
	setIndex $optic 1 [expr $i+1] BK7
	}
   set z [expr $z + $dz1]

#Back of corrector
   setSurf $optic 2 z $z
   setSurf $optic 2 curv [expr 1./$R2]
   loop i 0 $nwave {
	setIndex $optic 2 [expr $i+1] 1
	}
   set z [expr $z + $dz2]

#Primary
   setSurf $optic 3 z $z
   setSurf $optic 3 curv [expr 1./$R3]
   loop i 0 $nwave {
	setIndex $optic 3 [expr $i+1] -1
	}
   set z [expr $z + $dz3]

#Scondary
   setSurf $optic 4 z $z
   setSurf $optic 4 curv [expr 1./$R4]
   loop i 0 $nwave {
	setIndex $optic 4 [expr $i+1] 1
	}
   set z [expr $z + $dz4]

#Focal plane
   setSurf $optic 5 z $z
   setSurf $optic 5 curv [expr 1./$R5]
   loop i 0 $nwave {
	setIndex $optic 5 [expr $i+1] 1
	}
   colorcount $optic
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }

######################################################################
#Enter information for a lens and convert to cray format.
#Input format is inspired by http://www.atmsite.org/contrib/Fejes/
#	Optimization1/Optimization1.html

proc lens {} {
   puts stdout "Enter parameters for a lens design."
   while {1} {
	puts stdout "Primary diameter: " nonewline; flush stdout
	set D [gets stdin]
	if {![catch {format %f $D}]} break
	puts stdout "Illegal number - try again"
	}
   while {1} {
	puts stdout "Final focal ratio: " nonewline; flush stdout
	set fratio [gets stdin]
	if {![catch {format %f $fratio}]} break
	puts stdout "Illegal number - try again"
	}

   set flen [expr $D*$fratio]

#Initialize surfaces.
   set optic [opticNew]

   handleSet $optic.name "Lens"
   handleSet $optic.ncolor 1
   handleSet $optic.tel->diam $D
   handleSet $optic.tel->fr1 $fratio
   handleSet $optic.tel->fr2 $fratio
   handleSet $optic.tel->back 0
   handleSet $optic.tel->fl1 $flen
   handleSet $optic.tel->finner 0
   handleSet $optic.tel->mag 1
   handleSet $optic.tel->beta 0
   handleSet $optic.tel->sep 0
   handleSet $optic.tel->f3 $flen

#Size of focal plane: Use 60 arcminute diameter
   set size [expr abs($flen*(30./60./57.3))]

#Scale in arcsec/mm
   set scale [expr 206265./$flen]

   while {1} {
	puts stdout "Number of wavelengths: " nonewline; flush stdout
	set nwave [gets stdin]
	if {![catch {format %f $nwave}]} break
	puts stdout "Illegal number - try again"
	}
   set waves ""
   loop i 0 $nwave {
	while {1} {
	   puts stdout "Design wavelength [expr $i+1]: " nonewline
	   flush stdout
	   set wave [gets stdin]
	   if {![catch {format %f $wave}]} break
	   puts stdout "Illegal number - try again"
	   }
	lappend waves $wave
	}

   loop i 0 $nwave {
	set icolor [expr $i+1]
	setFocal $optic $icolor xpos 0
	setFocal $optic $icolor ypos 0
	setFocal $optic $icolor xrad $size
	setFocal $optic $icolor yrad $size
	setFocal $optic $icolor scale $scale
	setFocal $optic $icolor wave [lindex $waves $i]
	setFocal $optic $icolor dist 0
	setFocal $optic $icolor rot 0
	}

#Set wavelength in use flags.
   setSurf $optic 0 z -1.e14
   loop i 0 $nwave {
	setIndex $optic 0 [expr $i+1] 1
	}

#Input compound lens design
   while {1} {
	puts stdout "Number of elements: " nonewline; flush stdout
	set nelem [gets stdin]
	if {![catch {format %f $nelem}]} break
	puts stdout "Illegal number - try again"
	}

#Now input the surfaces parameters
   set z 0
   set nsurf 1
   loop i 0 $nelem {
	while {1} {
	   puts stdout "Surface $nsurf radius, thickness, glass: " nonewline
	   flush stdout
	   set list [gets stdin]
	   set R [lindex $list 0]
	   set dz [lindex $list 1]
	   set glass [lindex $list 2]
	   if {![catch {format %f $R}] && ![catch {format %f $dz}] && \
	       $glass != ""} break
	   puts stdout "Illegal number - try again"
	   }
	setSurf $optic $nsurf z $z
	setSurf $optic $nsurf curv [expr 1./$R]
	loop i 0 $nwave {
	   setIndex $optic $nsurf [expr $i+1] $glass
	   }
	set z [expr $z + $dz]
	incr nsurf
	}

   while {1} {
	puts stdout "Focal plane, Radius: " nonewline
	flush stdout
	set list [gets stdin]
	set R [lindex $list 0]
	if {![catch {format %f $R}]} break
	puts stdout "Illegal number - try again"
	}

   setSurf $optic $nsurf z $z
   if {$R != 0} {setSurf $optic $nsurf curv [expr 1./$R]}
   loop i 0 $nwave {
	setIndex $optic $nsurf [expr $i+1] 1
	}

   colorcount $optic
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }



