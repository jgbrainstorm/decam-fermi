#TCL replacement for readin.

proc opticRead {file} {
   global env
   global _optic
   set file [dataFileFind $file dat]
   set hndl [opticNew]
   set fid [open $file]
   handleSet $hndl.name [file tail $file]

#Leading comment is not stored.
   set comment [string range [gets $fid] 1 end]
   echo $comment
   set line [string trim [gets $fid]]

#Primary diameter
   handleSet $hndl.tel->diam [lindex $line 0]
   set line [string trim [gets $fid]]

#Primary mirror focal ration
   handleSet $hndl.tel->fr1 [lindex $line 0]
   set line [string trim [gets $fid]]

#Final focal ration
   handleSet $hndl.tel->fr2 [lindex $line 0]
   set line [string trim [gets $fid]]

#Back focal distance
   handleSet $hndl.tel->back [lindex $line 0]
   set line [string trim [gets $fid]]

#Primary mirror inner fraction
   handleSet $hndl.tel->finner [lindex $line 0]

#Need to call telcomp.
#I need fl1 in opticin
   handleSet $hndl.tel->fl1 [expr [exprGet $hndl.tel->diam] * [exprGet \
	$hndl.tel->fr1]]

#Final focal length - not really needed.
   handleSet $hndl.tel->f3 [expr [exprGet $hndl.tel->diam] * [exprGet \
	$hndl.tel->fr2]]

   opticzero $hndl
   handleSet $hndl.nsurf 0
   set comment ""
   set lastsurf 0
   set newsurf 0
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	if {[string index $line 0] == "#"} {
	   set comment [string range $line 1 end]
	   continue
	   }
	set isurf [lindex $line 0]
	set param [lindex $line 1]
	set val [lindex $line 2]

#If index of refraction is a glass name, convert to a numerical index.
#Must have input the wavelength previously, but I don't check for that.
	if {$isurf >= 0 && $param > 100 && $param <= 200} {
	   if {[catch {format %f $val}]} {
		set wave [getWave $hndl $param]
		set val [glass $val $wave]
		}
	   if {$isurf != $lastsurf} {
		set lastsurf $isurf
		set newsurf 1
		}
	   }
	setparam $hndl $isurf $param $val
	if {$newsurf == 1} {
	   set nsurf [exprGet $hndl.nsurf]
	   handleSet $hndl.optic<$nsurf>->comment $comment
	   set newsurf 0
	   set comment ""
	   }
	}
   close $fid

#Look for aperture stop - if none supplied, use surface 1.
   set iapp 0
   for {set isurf 1} {$isurf <= [exprGet $hndl.nsurf]} {incr isurf} {
	if {[exprGet $hndl.optic<$isurf>->stoptype] == 2} {
	   set iapp $isurf
	   break
	   }
	}
   if {$iapp == 0} {
	set iapp 1
	handleSet $hndl.optic<$iapp>->stoptype 2
	}

#Make sure aperture stop is non-zero, positive
   if {[exprGet $hndl.optic<$iapp>->outstop] <= 0} {
	handleSet $hndl.optic<$iapp>->outstop [expr [exprGet $hndl.tel->diam] \
	   /2.]
	}
   colorcount $hndl

#Convert to new index convention
   indexConvert $hndl

   opticinc $hndl 1

#Define ray pattern
   rayPattern $hndl 6 1

#Display first color
   if {[info command .r] != ""} {
	opticDisplay $hndl 1
	set nsurf [exprGet $hndl.nsurf]
	opticPlot $hndl 1 $nsurf 1
	}

   return $hndl
   }
