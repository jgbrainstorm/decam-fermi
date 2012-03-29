#Create tcl bindings that emulate the Dervish pgPlot bindings.

###############################################################
#Refresh screen from within getNextChar loop of ile.tcl
#This is needed really only for xwin device
proc localupdate {} {
   global _plplot
   if {![info exists _plplot(dev)]} return

#Need to cache device type, because calling plgdev causes plflush already!
   if {$_plplot(dev) == "xwin"} {plflush}
   }

###############################################################
proc pgUpdt {} {
   plflush
   return
   }

###############################################################
proc pgBegin {{type ""} {nx 1} {ny 1}} {
   global _plplot
   global INVERT


#If plot is open, close first
   if {[info exists _plplot(state)]} {
	if {$_plplot(state) == "open"} {
	   pgEnd
	   }
	}

#Cache the command - in case I need to rerun pgBegin from pgGeomSet
   set _plplot(command) [list pgBegin $type $nx $ny]
   set _plplot(nx) $nx
   set _plplot(ny) $ny

#If invert is set, we set background to white, foreground to black.
   if {$type == ""} {
	set file ""
	set dev xwin
   } else {
	set line [split $type /]

#If file name includes slashes, need to accommodate here
	if {[llength $line] == 1} {lappend line xwin}
	set file [join [lrange $line 0 end-1] /]
	set dev [string tolower [lindex $line end]]
	}

#I comvert file device to one that plplot likes below.  Keep track of
#original format - this is useful if I need to do a convert of file types
#at end.
   set _plplot(outdev) $dev

#Don't know if I can get a listing of all devices, unlike pgPlot
   if {[string range $dev 0 1] == "xw"} {set dev xwin}

#plplot calls it "pbm" (which is 1 bit deep) but really creates "ppm".
#I will allow user to input "ppm".
   if {$dev == "ppm"} {
	set dev pbm
	if {$file == ""} {
	   error "Must specify output file name"
	   }
	}

#If user inputs "jpg", I will create ppm, then convert using cjpeg later.
#This is easier than trying to compile in support for all graphics formats
#in the main program - too many libraries to track!
   if {$dev == "jpg"} {
	set dev pbm
	if {$file == ""} {
	   error "Must specify output file name"
	   }
	set root [file root $file]
	set file $root.ppm
	}

#Likewise for gif - conversion program is ppmtogif
   if {$dev == "gif"} {
	set dev pbm
	if {$file == ""} {
	   error "Must specify output file name"
	   }
	set root [file root $file]
	set file $root.ppm
	}

#Likewise for png - conversion program is pnmtopng
   if {$dev == "png"} {
	set dev pbm
	if {$file == ""} {
	   error "Must specify output file name"
	   }
	set root [file root $file]
	set file $root.ppm
	}

#Color v. mono postscript
   if {$dev == "cps"} {
	set dev psc
	set INVERT 1
	}

   if {$dev == "vcps"} {
	set dev vpsc
	set INVERT 1
	}

#vps for portrait-mode postscript
   set orient 0
   if {$dev == "vps"} {
	set orient 1
	set dev ps
	}
   if {$dev == "vpsc"} {
	set orient 1
	set dev psc
	}

   if {$file != ""} {
	plsfnam $file
	}
   plsdev $dev
   if {$orient == 1} {
	plsetopt -portrait
	}

   set _plplot(dev) $dev
   set _plplot(file) $file

   if {[info exists INVERT]} {
	plscolbg 255 255 255
	plscol0 15 0 0 0
	}

#   plinit
#Allow multiple plots per page
   plstar $nx $ny
   set _plplot(state) open

#Default white fg
   plcol0 15
   plspause 0
   pgNext
   return
   }

###############################################################
#Advance to next subpage without clearing the screen.
#pgPage will clear the screen.
proc pgNext {} {
   pladv
   plvpor .15 .9 .15 .9
   return
   }

###############################################################
proc pgEnv {xmin xmax ymin ymax arg1 arg2} {
   plenv $xmin $xmax $ymin $ymax $arg1 $arg2
   return
   }

###############################################################
proc pgAsk {arg} {
   plspause $arg
   return
   }

###############################################################
#Convert pgplot index to plplot index
proc pgcolToPlcol {index} {
   set out 15
   if {![ctype digit $index]} {error "Input $index must be 0-7"}
   if {$index == 0} {set out 0}
   if {$index == 1} {set out 15}
   if {$index == 2} {set out 1}
   if {$index == 3} {set out 3}
   if {$index == 4} {set out 9}
   if {$index == 5} {set out 11}
   if {$index == 6} {set out 13}
   if {$index == 7} {set out 2}
   return $out
   }

###############################################################
#Set color index.  In pgplot, we have 8 color indices.  I will do a horrible
#hack and assume input is original pgPlot convention; convert to plplot
#equivalent.

proc pgSci {index} {
   global _pgplot
   set _pgplot(ci) $index
   set out [pgcolToPlcol $index]
   plcol0 $out
   return
   }

###############################################################
#Set fill style.  The only way that an area is filled is via the pgPoly
#command.  So we cache the style index, and use it to direct other
#commands to use either pgPoly or pgLine.  With a bit more work, I can
#pass the index on to plpoly as well.

proc pgSfs {index} {
   global _plplot
   if {![ctype digit $index]} {error "In pgSfs, $index must be 1 or 2"}
   set _plplot(fill) $index
   return
   }

###############################################################
proc pgQfs {} {
   global _plplot
   if {![info exists _plplot(fill)]} {set _plplot(fill) 1}
   return $_plplot(fill)
   }
   
###############################################################
#Query for information.
proc pgQinf {item} {
   global _plplot
   set item [string tolower $item]
   if {$item == "state"} {
	if {[info exists _plplot(state)]} {
	   return "OPEN"
	} else {
	   return "CLOSED"
	   }
	}
   if {$item == "file"} {
	if {[info exists _plplot(file)]} {
	   return $_plplot(file)
	} else {
	   return
	   }
	}
   if {$item == "dev" || $item == "device" || $item == "type"} {
	if {[info exists _plplot(dev)]} {
	   set dev $_plplot(dev)
	   if {$dev == "xwin"} {set dev xwindow}
	   return [string toupper $dev]
	} else {
	   return
	   }
       }
   return
   }

###############################################################
#Query color index.

proc pgQci {} {
#Gak! plplot does not have a routine to query for the current color index!
#Return 1 (white) just to be friendly
   global _pgplot
   if {[info exists _pgplot(ci)]} {
	return $_pgplot(ci)
   } else {
	return 1
	}
   }

############################################################################
#Set color representation

proc pgScr {index r g b} {
   set out [pgcolToPlcol $index]

#Rescale to 0 to 255
   foreach var "r g b" {
	set $var [expr int([set $var]*255+0.5)]
	if {[set $var] < 0} {set $var 0}
	if {[set $var] > 255} {set $var 255}
	}
   plscol0 $out $r $g $b
   return
   }

############################################################################
#Query for color representation

proc pgQcr {index} {
   set out [pgcolToPlcol $index]
   set list [plgcol0 $out]

#Convert to floating point
   set r [format %.4f [expr [lindex $list 0]/255.]]
   set g [format %.4f [expr [lindex $list 1]/255.]]
   set b [format %.4f [expr [lindex $list 2]/255.]]
   return "$r $g $b"
   }

######################################################################
#Set character scale (seems to work only work for some text commands).
proc pgSch {scale} {
   global _pgplot
   set _pgplot(ch) $scale
   plschr 0 $scale
   plssym 0 $scale
   return
   }

#####################################################################
proc pgQch {} {
   global _pgplot
   if {[info exists _pgplot(ch)]} {
	return $_pgplot(ch)
   } else {
	return 1
	}
   }

######################################################################
proc pgSlw {width} {
   global _pgplot
   set _pgplot(lw) $width
   plwid $width
   return
   }

######################################################################
proc pgSls {style} {
   global _pgplot
   set _pgplot(ls) $style
   pllsty $style
   return
   }

######################################################################
proc pgQls {} {
   global _pgplot
   if {[info exists _pgplot(ls)]} {
	return $_pgplot(ls)
   } else {
	return 1
	}
   }

##########################################################################
#Get line width
proc pgQlw {} {
   global _pgplot
   if {[info exists _pgplot(lw)]} {
	return $_pgplot(lw)
   } else {
	return 1
	}
   }

##########################################################################
#Draw a line
proc pgLine {xpts ypts} {
   if {[llength $xpts] != [llength $ypts]} {
	error "x and y lists must be same length"
	}
   if {[info command _x] != ""} {_x delete}
   if {[info command _y] != ""} {_y delete}
   matrix _x f [llength $xpts]
   matrix _y f [llength $ypts]
   loop i 0 [llength $xpts] {
	_x $i = [lindex $xpts $i]
	}
   loop i 0 [llength $ypts] {
	_y $i = [lindex $ypts $i]
	}
   plline [llength $xpts] _x _y
   _x delete
   _y delete
   return
   }

##########################################################################
#Draw a set of points
proc pgPoint {xpts ypts sym} {
   if {[llength $xpts] != [llength $ypts]} {
	error "x and y lists must be same length"
	}
   if {[info command _x] != ""} {_x delete}
   if {[info command _y] != ""} {_y delete}
   matrix _x f [llength $xpts]
   matrix _y f [llength $ypts]
   loop i 0 [llength $xpts] {
	_x $i = [lindex $xpts $i]
	}
   loop i 0 [llength $ypts] {
	_y $i = [lindex $ypts $i]
	}
   plpoin [llength $xpts] _x _y $sym
   _x delete
   _y delete
   return
   }

##########################################################################
#Draw a polygon
proc pgPoly {xpts ypts} {
   if {[llength $xpts] != [llength $ypts]} {
	error "x and y lists must be same length"
	}
   if {[info command _x] != ""} {_x delete}
   if {[info command _y] != ""} {_y delete}

#Append beginning onto end
   lappend xpts [lindex $xpts 0]
   lappend ypts [lindex $ypts 0]

   matrix _x f [llength $xpts]
   matrix _y f [llength $ypts]
   loop i 0 [llength $xpts] {
	_x $i = [lindex $xpts $i]
	}
   loop i 0 [llength $ypts] {
	_y $i = [lindex $ypts $i]
	}

#I should not fill in if fill style is 2
   if {[pgQfs] != 2} {
	plfill [llength $xpts] _x _y
   } else {
	plline [llength $xpts] _x _y
	}
   _x delete
   _y delete
   return
   }

#######################################################################
proc pgRect {x1 x2 y1 y2} {
   set index [pgQfs]
   if {$index == 1} {
	set cmd pgPoly
   } else {
	set cmd pgLine
	}
   set xlist ""
   set ylist ""
   lappend xlist $x1 $x2 $x2 $x1 $x1
   lappend ylist $y1 $y1 $y2 $y2 $y1
   $cmd $xlist $ylist
   return
   }

#######################################################################
proc pgCirc {xcen ycen rad} {
   set index [pgQfs]
   if {$index == 1} {
	set cmd pgPoly
   } else {
	set cmd pgLine
	}
   set xlist ""
   set ylist ""
   set nstep 100
   loop i 0 [expr $nstep+1] {
	set ang [expr 2.*3.14159*$i/$nstep]
	set x [expr $xcen + $rad*cos($ang)]
	set y [expr $ycen + $rad*sin($ang)]
	lappend xlist $x
	lappend ylist $y
	}
   $cmd $xlist $ylist
   return
   }

#######################################################################
proc pgEnd {} {
   global _plplot
   plend

   if {![info exists _plplot]} return

#For certain formats, need to run a conversion program at the end.
   if {$_plplot(outdev) == "jpg"} {
	set root [file root $_plplot(file)]
	set file $root.jpg
	set command [auto_execok cjpeg]
	if {$command != ""} {
	   echo Converting $_plplot(file) to $file
	   set code [catch {exec $command -sample 1x1 $_plplot(file) > $file} ]
	   if {$code == 0} {
		exec rm $_plplot(file)
		}
	} else {
	   echo Conversion program cjpeg not found\; file in ppm format
	   }
	}
   if {$_plplot(outdev) == "gif"} {
	set root [file root $_plplot(file)]
	set file $root.gif
	set command [auto_execok ppmtogif]
	if {$command != ""} {
	   echo Converting $_plplot(file) to $file
	   set code [catch {exec $command $_plplot(file) > $file}]

#For some reason the return code on my machine is always 1, even for success.
	   if {$code == 0 || $code == 1} {
		exec rm $_plplot(file)
		}
	} else {
	   echo Conversion program ppmtogif not found\; file in ppm format
	   }
	}
   if {$_plplot(outdev) == "png"} {
	set root [file root $_plplot(file)]
	set file $root.png
	set command [auto_execok pnmtopng]
	if {$command != ""} {
	   echo Converting $_plplot(file) to $file
	   set code {exec $command $_plplot(file) > $file}
	   if {$code == 0} {
		exec rm $_plplot(file)
		}
	} else {
	   echo Conversion program pnmtopng not found\; file in ppm format
	   }
	}
   if {[info exists _plplot(state)]} {unset _plplot(state)}
   unset _plplot
   return
   }

#######################################################################
proc pgVport {{xmin .15} {xmax .9} {ymin .15} {ymax .9}} {
   plvpor $xmin $xmax $ymin $ymax
   return
   }

#######################################################################
proc pgWindow {xmin xmax ymin ymax} {

#plwind seems to need a viewport, which seems to need a new page.
#I will call pladv in pgBegin or pgPage
   plwind $xmin $xmax $ymin $ymax
   return
   }

#######################################################################
#Default box for now - need to figure out pgBox options, if any.
proc pgBox {{xopt bcnst} {xtick 0} {nxsub 0} {yopt bcnstv} {ytick 0} \
      {nysub 0}} {
   plbox $xopt $xtick $nxsub $yopt $ytick $nysub
   return
   }

#####################################################################
proc pgLabel {x y t} {
   pllab $x $y $t
   return
   }

####################################################################
#Write text relative to world coordinates
#Allow optional justification - pgPlot extension!
proc pgText {x y text {just 0}} {
   plptex $x $y 1 0 $just $text
   return
   }

####################################################################
#Write text relative to viewport boundaries.
proc pgMtext {side displ pos just text} {
   plmtex $side $displ $pos $just $text
   return
   }

###################################################################
#Set geometry - pixel devices only
proc pgGeomSet {xwid ywid xoff yoff} {
   global _plplot

#Hmm, in plplot, one sets geometry before calling plinit, which makes more
#sense than Dervish way, but requires some monkeying about
   plend
   set _plplot(state) closed
   plspage $xwid $ywid $xwid $ywid $xoff $yoff

#Rerun pgBegin
   if {[info exists _plplot(command)]} {
	eval $_plplot(command)
	}
   return
   }

#######################################################################
#End of page
proc pgPage {} {
   pleop
   pladv
   plvpor .15 .9 .15 .9
   return
   }

#########################################################################
#Create a new canvas

proc canvInit {{canv ""}} {
   set top ""
   if {$canv == "" || [winfo class $canv] != "Canvas"} {
	foreach letter "a b c d e f g h i j k l m n o p q r s t u v w x y z" {
	   if {[info command .$letter] == ""} {
	      set top .$letter
	      break
	      }
	   }
	if {$top == ""} {
	   echo Cannot find next toplevel window!
	   return
	   }
	toplevel $top
	set canv $top.c
	canvas $canv -width 640 -height 480 -bg black
	pack $canv -fill both -expand 1
   } else {
	pgEnd
	}
   pgBegin $canv/tk
   return $canv
   }


###########################################################################
#Some viewport and window queries

#Get viewport coordinates in device land.  This is all I support.
#pgPlot allows querying for physical coords as well.  Yuk!

#NOTE!!! derivsh returns the info in a list with surrounding braces, which
#is rather pointless.  My pgQvp does not include the braces.
#Also, dervish pgQvp needs a unit specifier; I will just default to 0.

proc pgQvp {{unit 0}} {
   if {$unit == 0} {
	return [plgvpd]
	}
   error "pgQvp: Unit specifier $unit is invalid!"
   }

#########################################################################
#Get window coord limits

proc pgQwin {} {
   return [plgvpw]
   }

#########################################################################
#Get page size.
#In pgPlot, the command "pgQvsz 3" will return the page size in pixels.
#I am not sure what happens when one creates multiple viewports (plots)
#on a page.  For now, I will create a fictitious command to return the
#page size.

#Return params are resx, resy, xlen, ylen, xoff, yoff

proc pgQpsz {} {
   return [plgpage]
   }
