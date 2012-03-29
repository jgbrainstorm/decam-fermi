############################################################################
#Convert oslo format to my "zlist" format.
#DANGER WILL ROBINSON!!!  CRAY and OSLO  both use the .len extension.
#Since I am just reading, I will make the bold assumption that I need
#not care.

proc oslo2zlist {file {outfile ""}} {
   set ext [file extension $file]
   if {$ext == ""} {
	set file $file.len
	}
   if {![file exists $file]} {
	error "File $file not found!"
	}
   set root [file root [file tail $file]]
##   set outfile $root.txt

#Read and store oslo prescription first.
   set fid [open $file]
   set design() ""
   set field ""

#For sanity, check for "magic number" on 1st line
   set line [gets $fid]
   if {![regexp {//[ ]+OSLO} $line]} {
	close $fid
	error "Does not appear to be an OSLO file, based on 1st line"
	}

#Skip all ancillary info, including wavelength for now.
   set surf 0

#It seems that default stop surface is 1.
   set astopsurf 1
   set type STANDARD
   set scale 1.
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set command [lindex $line 0]

#A new lens is initialized with the "LEN" command, and surface 0 is
#created immediately.
	if {$command == "LEN" || $command == "NXT"} {
	   if {$command == "LEN"} {
		set surf 0
	   } else {
		incr surf
		}
	   lappend design() $surf
	   set design($surf,curv) 0
	   set design($surf,zthick) 0
	   set design($surf,diam) 10
	   set design($surf,type) STANDARD
	   set design($surf,glass) air
	   set design($surf,name) S$surf
	   set type STANDARD
	   continue
	   }

	if {$command == "PK"} {
	   set var [lindex $line 1]
	   set offset [lindex $line 2]
	   set surf1 [expr $surf + $offset]

#Don't know what fact1 does.  fact2 appears to be a multiplier
	   set fact1 [lindex $line 3]
	   set fact2 [lindex $line 4]
	   if {$var == "THM"} {
		set var TH
		set fact2 -1.
		}
	   if {[info exists design($surf1,$var)]} {
		set command $var
		set param $design($surf1,$var)
		if {$fact2 != ""} {
		   set param [expr $param*$fact2]
		   }
		set line [list $command $param]
		}
	   }

	set design($surf,$command) [lindex $line 1]

#Don't know if this is correct interpretation of units.

	if {$command == "UNI"} {
	   set scale [lindex $line 1]
	} elseif {$command == "AST"} {
	   set design($surf,stop) 1
	   set astopsurf $surf
	} elseif {$command == "RD"} {
	   set rd [lindex $line 1]
	   set rd [expr $rd*$scale]
	   if {$rd != 0.} {
		set curv [expr 1./$rd]
	   } else {
		set curv 0.
		}
	   set design($surf,curv) $curv
	   set design($surf,CV) $curv
	} elseif {$command == "CV"} {
	   set curv [lindex $line 1]
	   set curv [expr $curv/$scale]
	   set design($surf,curv) $curv
	   if {$curv != 0.} {
		set rd [expr 1./$curv]
	   } else {
		set rd 0.
		}
	   set design($surf,RD) $rd
	} elseif {$command == "TH"} {
	   set zthick [lindex $line 1]
	   set zthick [expr $zthick*$scale]
	   set design($surf,zthick) $zthick
	} elseif {$command == "CC"} {
	   set ccon [lindex $line 1]
	   set design($surf,ccon) $ccon
	} elseif {$command == "AP"} {

#What zemax calls "diam" is actually the radius.  OSLO also uses radius.
#I will perpetuate the weird Zemax convention.
	   set diam [lindex $line 1]
	   if {$diam == "CHK"} {
		set diam [lindex $line 2]
		}
	   set diam [expr $diam*$scale]
	   set design($surf,diam) $diam
	} elseif {$command == "GLA"} {
	   set glass [lindex $line 1]

#Translate name.  A lot more could be done here.
	   if {$glass == "F_SILICA"} {set glass Fused-Silica}
	   global glassName
	   if {![info exists glassName]} opusName
	   if {[info exists glassName($glass)]} {
		set glass $glassName($glass)
		}
	   set design($surf,glass) $glass
	} elseif {$command == "RFL"} {
	   set design($surf,glass) mirror

#Is this also a mirror?
	} elseif {$command == "RFH"} {
	   set design($surf,glass) mirror
	} elseif {$command == "AC"} {
	   set design($surf,a2) [expr [lindex $line 1] * pow($scale,4)]
	} elseif {$command == "AD"} {
	   set design($surf,a4) [expr [lindex $line 1] * pow($scale,4)]
	} elseif {$command == "AE"} {
	   set design($surf,a6) [expr [lindex $line 1] * pow($scale,6)]
	} elseif {$command == "AF"} {
	   set design($surf,a8) [expr [lindex $line 1] * pow($scale,8)]
	} elseif {$command == "AG"} {
	   set design($surf,a10) [expr [lindex $line 1] * pow($scale,10)]
	} elseif {$command == "ASP"} {
	   set type "EVENASPH"
	   set design($surf,type) $type
	} elseif {$command == "A1" && $type == "EVENASPH"} {
	   set design($surf,a2) [expr [lindex $line 1] * pow($scale,2)]
	} elseif {$command == "A2" && $type == "EVENASPH"} {
	   set design($surf,a4) [expr [lindex $line 1] * pow($scale,4)]
	} elseif {$command == "A3" && $type == "EVENASPH"} {
	   set design($surf,a6) [expr [lindex $line 1] * pow($scale,6)]
	} elseif {$command == "A4" && $type == "EVENASPH"} {
	   set design($surf,a8) [expr [lindex $line 1] * pow($scale,8)]
	} elseif {$command == "A5" && $type == "EVENASPH"} {
	   set design($surf,a10) [expr [lindex $line 1] * pow($scale,10)]
	} elseif {$command == "AS1"} {
	   set design($surf,a2) [expr [lindex $line 1] * pow($scale,2)]
	} elseif {$command == "AS2"} {
	   set design($surf,a4) [expr [lindex $line 1] * pow($scale,4)]
	} elseif {$command == "AS3"} {
	   set design($surf,a6) [expr [lindex $line 1] * pow($scale,6)]
	} elseif {$command == "AS4"} {
	   set design($surf,a8) [expr [lindex $line 1] * pow($scale,8)]
	} elseif {$command == "AS5"} {
	   set design($surf,a10) [expr [lindex $line 1] * pow($scale,10)]
	   }
	}

#In case the design did not expicitly set the aperture stop ...
   set design($astopsurf,stop) 1

#Last surface is the image surface
   set imagesurf $surf
   close $fid

#Temporary
   if {$outfile != ""} {
	set stdout [open $outfile w]
   } else {
	set stdout stdout
	}

   set isurf 0
   puts $stdout "Surface        Radius       Thickness      \
           Glass                Name"

   foreach surf $design() {
	if {$surf == 0} continue
	if {$design($surf,type) == "COORDBRK"} continue
	incr isurf
	set curv $design($surf,curv)
	if {$curv != 0} {
	   set radius [expr 1./$curv]
	} else {
	   set radius 0.
	   }
	set thick $design($surf,zthick)
	set glass $design($surf,glass)
	set name $design($surf,name)
	set line [format "%3d   %15.6f %15.3f %15s %15s " $isurf $radius \
	    $thick $glass $name]
	puts $stdout $line
	foreach var "ccon a2 a4 a6 a8 a10" {
	   if {[info exists design($surf,$var)]} {
		set val $design($surf,$var)
		if {$val != 0.} {
		   puts $stdout "       $var = [format %.7g $val]"
		   }
		}
	   }
	set type $design($surf,type)
	if {[info exists design($surf,stop)]} {
	   set radius $design($surf,diam)
	   puts $stdout "       astop = [format %.7g $radius]"
	   }
	if {$surf == $imagesurf} {
	   set radius $design($surf,diam)
	   puts $stdout "       istop = [format %.7g $radius]"
	   }
	}
   if {$outfile != ""} {close $stdout}
   return
   }

######################################################################
#Read a zemax file in one go.  Use a temporary file to store data.

proc osloRead {file} {
   set design [file root [file tail $file]]
   set tempfile /tmp/[clock seconds].txt
   if {[catch {oslo2zlist $file $tempfile} msg]} {
	catch {exec rm $tempfile}
	echo $msg
	error "Failure!"
	}
   if {[catch {set hndl [zread $tempfile]} msg]} {
	catch {exec rm $tempfile}
	echo $msg
	error "Failure!  Check for aberrant OPTIC handle"
	}
   catch {exec rm $tempfile}
   handleSet $hndl.name $design
   return $hndl
   }
