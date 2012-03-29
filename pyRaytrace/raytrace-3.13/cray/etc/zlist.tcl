#List a design in a format that can be copied into Zemax.
#If flag is 0, list thickness.  If flag is 1, list absolute z pos.

proc zlist {hndl icolor {file ""} {flag 0}} {
   set nsurf [exprGet $hndl.nsurf]
   set surfids ""
   if {$file != ""} {
	set stdout [open $file w]
   } else {
	set stdout stdout
	}
    for {set i 0} {$i <= $nsurf} {incr i} {
	set surfid [surfId $hndl $i]
	set index [showIndex $hndl $surfid $icolor]
	if {$index == 0} continue
	lappend surfids $surfid
	}
   if {$flag == 0} {
	puts $stdout "Surface        Radius       Thickness      \
	   Glass                Name"
   } else {
	puts $stdout "Surface        Radius       ZPos           \
	   Glass                Name"
	}
   set index0 1
   loop i 0 [expr [llength $surfids]] {
	set surfid [lindex $surfids $i]
	set surfid1 [lindex $surfids [expr $i+1]]
	foreach var "curv ccon a2 a4 a6 a8 a10 theta" {
	   set $var [showSurf $hndl $surfid $var]
	   }
	set z [showSurf $hndl $surfid z]
	if {$surfid1 != ""} {
	   set z1 [showSurf $hndl $surfid1 z]
	} else {
	   set z1 $z
	   }
	set thick [expr $z1 - $z]
	if {$curv != 0.} {
	   set radius [expr 1./$curv]
	} else {
	   set radius 0.
	   }
	set index [showIndex $hndl $surfid $icolor]
	if {abs($index) > 1} {
	   set glass [showGlass $hndl $surfid]
	   if {$glass == ""} {set glass UNKNOWN}
	   if {$glass == "SIO2"}  {set glass "Fused-Silica"}
	} elseif {$index0*$index < 0} {
	   set glass reflect
	} else {
	   set glass air
	   }
	set index0 $index
	set name [showName $hndl $surfid]
	if {$flag == 0} {
	   puts $stdout [format "%3d   %15.6f %15.3f %15s %15s " $i $radius \
		$thick $glass $name]
	} else {
	   puts $stdout [format "%3d   %15.6f %15.3f %15s %15s " $i $radius \
		$z $glass $name]
	   }
	foreach var "ccon a2 a4 a6 a8 a10 theta" {
	   if {[set $var] != 0.} {
		puts $stdout "       $var = [format %.7g [set $var]]"
		}
	   }
	set stoptype [showSurf $hndl $surfid stoptype]
	if {$stoptype == 2} {
	   set radius [format %.1f [showSurf $hndl $surfid outstop]]
	   puts $stdout "       astop = $radius"
	   }
	if {$stoptype == 5} {
	   set lines [showSurf $hndl $surfid lines]
	   set order [showSurf $hndl $surfid order]
	   puts $stdout "       grating = $lines $order"
	   }

#Last surface - print outer stop as an "image stop"
	if {$surfid1 == ""} {
	   set radius [format %.1f [showSurf $hndl $surfid outstop]]
	   puts $stdout "       istop = $radius"
	   }
	}

#Print out wavelengths.  OK, do I want to print everything, or just
#the input wavelength?  For now, I will print everything.
   set ncolor [exprGet $hndl.ncolor]
   for {set ifil 1} {$ifil <= $ncolor} {incr ifil} {
	set flag ""
	if {$ifil == $icolor} {
	   set flag " X"
	   }
	set wave [showFocal $hndl $ifil wave]
	puts $stdout "wave = $wave$flag"
	}
   if {$file != ""} {close $stdout}
   return
   }

##########################################################################
#Helper proc

proc waveRead {fid} {
   set waveList ""
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {[string index $line 0] == "#"} continue
	if {[llength $line] == 0} continue
	set wave [lindex $line 1]
	lappend waveList $wave
	}
   return $waveList
}

##########################################################################
#Read back zlist file.
#For now, only read back files with thickness.
#First line is a header
#Next lines are surfid, curv, thick, glass
#   or param = val
#param can be astop, istop, ccon, a4, a6, etc.

#I need to convert inches to mm.  Let's do it here.
#Scale = 25.4 to make that happen

proc zread {file {scale 1.}} {
   set ext [file extension $file]

#Default extension is .txt
   if {$ext == ""} {set file $file.txt}
   if {![file exists $file]} {
	error "No file found $file"
	}
   set fid [open $file]

#First line has titles.
   set line [gets $fid]
   set optic [opticNew]

   handleSet $optic.name "$file"
   handleSet $optic.ncolor [array size waveList]

#I am now writing out and reading back surface 0 (object plane), so the
#following should get overwritten.  I'll keep just in case ...
   setSurf $optic 0 z -1.e14

   set z 0
   set surf 0
   set waveList ""
   set surfList ""
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set var [lindex $line 0]

#If surf is a number, this is a new surface.  If surf is alphanumeric,
#it is a modifier for the current surface.
	if {[ctype digit $var]} {
	   set surf $var
	   if {[lsearch $surfList $surf] < 0} {
		lappend surfList $surf
		}
	   set rad [lindex $line 1]
	   set thick [lindex $line 2]
	   set glass [lindex $line 3]

#Apply scaling
	   set rad [expr $rad*$scale]
	   set thick [expr $thick*$scale]
	   if {$rad != 0.} {
		set curv [expr 1./$rad]
	   } else {
		set curv 0.
		}

#Revert back some of my conventions.
	   if {$glass == "Fused-Silica"} {set glass SIO2}
	   if {$glass == "Fused-Silica"} {set glass SIO2}
	   if {$glass == "F-SILICA"} {set glass SIO2}
	   if {$glass == "SILICA"} {set glass SIO2}
	   if {$glass == "reflect"} {set glass mirror}

	   set name [lindex $line 4]
	   setSurf $optic $surf z $z
	   setSurf $optic $surf curv $curv
	   setGlass $optic $surf $glass
	   setName $optic $surf $name
	   set z [expr $z + $thick]
	} else {
	   set val [lindex $line 2]
	   set val2 [lindex $line 3]
	   foreach p "ccon a2 a4 a6 a8 a10 theta" {
		if {"$var" == "$p"} {
		   if {$var == "a2"} {
			set val [expr $val/pow($scale,1)]
			}
		   if {$var == "a4"} {
			set val [expr $val/pow($scale,3)]
			}
		   if {$var == "a6"} {
			set val [expr $val/pow($scale,5)]
			}
		   if {$var == "a8"} {
			set val [expr $val/pow($scale,7)]
			}
		   if {$var == "a10"} {
			set val [expr $val/pow($scale,9)]
			}
		   setSurf $optic $surf $var $val
		   }
		}

#Aperture stop
	   if {$var == "astop"} {
		set val [expr $val*$scale]
		setSurf $optic $surf stoptype 2
		setSurf $optic $surf outstop $val

#I will set z of my aperture stop to 0 later.
		set zapp [showSurf $optic $surf z]
		}

#Grating
	   if {$var == "grating"} {
		setSurf $optic $surf stoptype 5
		setSurf $optic $surf lines $val
		setSurf $optic $surf order $val2
		}

#Stop in image plane.
	   if {$var == "istop"} {
		set val [expr $val*$scale]
		set istop $val
		setSurf $optic $surf outstop $val
		}
#Wavelength.  Can be anywhere in the file
	   if {$var == "wave"} {
		lappend waveList $val
		}
	   }
	}
   close $fid
   if {![info exists istop]} {set istop 10.}

#Focal plane - fill in fixed quantities.

#Use preset wavelengths if none supplied.
   if {$waveList == ""} {
	set waveList ".39 .50 .65 .90"
	}

   set ifil 0
   foreach wave $waveList {
	incr ifil
	setFocal $optic $ifil xpos 0
	setFocal $optic $ifil ypos 0
	setFocal $optic $ifil wave $wave
	setFocal $optic $ifil weight 1

#The following is not wanted if I specifically input surface 0.
	if {[lsearch $surfList 0] < 0} {
	   setIndex $optic 0 $ifil 1
	   }

	foreach surfid $surfList {
	   set glass [showGlass $optic $surfid]
	   set index [glass $glass $wave]
	   setIndex $optic $surfid $ifil $index
	   }

#Remaining focal plane parameters
	setFocal $optic $ifil xrad $istop
	setFocal $optic $ifil yrad $istop
	}

#Reposition z
   for {set i 0} {$i <= [exprGet $optic.nsurf]} {incr i} {
	set surfid [surfId $optic $i]
	setSurf $optic $surfid z [expr [showSurf $optic $surfid z] - $zapp]
	}

   rayPattern $optic 6 1
   colorcount $optic

#Scale factor.  Different from input scale!
#Use a single value for all filters.
   set sum 0.
   set n 0
   set ifil 0
   foreach wave $waveList {
	incr ifil
	opticInfo $optic $ifil
	set fl [showFocal $optic $ifil fl]
	set scale [expr 206265./$fl]
	set sum [expr $sum + $scale]
	incr n
	}

#If scale factor is negative, difficult to discern here.
   set scale [format %.3f [expr $sum/$n]]

   set ifil 0
   foreach wave $waveList {
	incr ifil
	setFocal $optic $ifil scale $scale
	}
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }

############################################################################
#Convert zemax format to my "zlist" format.

proc zemax2zlist {file {outfile ""}} {
   set ext [file extension $file]
   if {$ext == ""} {
	set file $file.zmx
	}
   if {![file exists $file]} {
	error "File $file not found!"
	}
   set root [file root [file tail $file]]
##   set outfile $root.txt

#Read and store zemax prescription first.
   set fid [open $file]
   set design() ""
   set field ""

#Skip all ancillary info, including wavelength for now.
   set surf 0
   set type STANDARD
   set scale 1.
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set command [lindex $line 0]
	if {$command == "SURF"} {
	   set surf [lindex $line 1]
	   lappend design() $surf
	   set design($surf,curv) 0
	   set design($surf,zthick) 0
	   set design($surf,diam) 10
	   set design($surf,type) STANDARD
	   set design($surf,glass) air
	   set design($surf,name) S$surf
	   set type STANDARD
	} elseif {$command == "UNIT"} {
	   set unit [lindex $line 1]
	   if {$unit == "MM"} {
		set scale 1.
	   } elseif {$unit == "CM"} {
		set scale 10.
	   } elseif {$unit == "IN"} {
		set scale 25.4
		}
	} elseif {$command == "TYPE"} {
	   set type [lindex $line 1]
	   set design($surf,type) $type
	} elseif {$command == "STOP"} {
	   set design($surf,stop) 1
	} elseif {$command == "CURV"} {
	   set curv [lindex $line 1]
	   set curv [expr $curv/$scale]
	   set design($surf,curv) $curv
	} elseif {$command == "DISZ"} {
	   set zthick [lindex $line 1]
	   if {$zthick == "INFINITY"} {
		set zthick 1.e14
	   } else {
		set zthick [expr $zthick*$scale]
		}
	   set design($surf,zthick) $zthick
	} elseif {$command == "CONI"} {
	   set ccon [lindex $line 1]
	   set design($surf,ccon) $ccon
	} elseif {$command == "DIAM"} {
	   set diam [lindex $line 1]
	   set diam [expr $diam*$scale]
	   set design($surf,diam) $diam
	} elseif {$command == "GLAS"} {
	   set glass [lindex $line 1]

#Translate name.  A lot more could be done here.
	   if {$glass == "MIRROR"} {set glass reflect}
	   if {$glass == "F_SILICA"} {set glass Fused-Silica}
	   global glassName
	   if {![info exists glassName]} opusName
	   if {[info exists glassName($glass)]} {
		set glass $glassName($glass)
		}
	   set design($surf,glass) $glass
	} elseif {$command == "PARM"} {
	   if {![info exists design($surf,type)]} continue
	   if {$design($surf,type) != "EVENASPH"} continue
	   set num [lindex $line 1]
	   set param a[expr 2*$num]
	   set val [lindex $line 2]
	   set val [expr $val/pow($scale,2*$num-1)]
	   if {$val != 0} {set design($surf,$param) $val}
	} elseif {$command == "COMM"} {
	   set name [lindex $line 1]
	   set design($surf,name) $name
	   }
	}

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
	foreach var "ccon a2 a4 a6 a8 a10 theta" {
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

proc zemaxRead {file} {
   set design [file root [file tail $file]]
   set tempfile /tmp/[clock seconds].txt
   if {[catch {zemax2zlist $file $tempfile} msg]} {
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
