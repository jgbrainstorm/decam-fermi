##############################################################
#Plotting setup routines.  Use a Tk window for plotting.
#
##############################################################
#If I supply a subwindow, assume I am plotting to a Tk canvas.  Else,
#assume I want X window.
proc plotInit {{subwin ""} {nx 1} {ny 1}} {
   if {[pgQinf STATE] == "OPEN"} {
	set device [pgQinf device]
	set type [pgQinf type]
	if {$type == "XWINDOW"} {
	   pgAsk 0
	   pgPage
	} elseif {$type == "TK"} {

#If no subwindow supplied, continue using existing window.
	   if {$subwin != ""} {
		pgEnd
#subwin must exist
		if {[info command .r.bottom.$subwin] == ""} {
		   error "Subwindow $subwin does not exist"
		   }
		pgBegin .r.bottom.$subwin/TK $nx $ny
		pgAsk 0
		pgPage
		}
	   }
#For any other device type already open, such as postscript or gif,
#don't do anything else.
   } else {

#If no subwindow specified, use X
	if {$subwin == ""} {
	   pgBegin "" $nx $ny
	   pgGeomSet 400 400 10 25
	   pgUpdt
	   pgAsk 0
	   pgPage
	} else {
	   if {[info command .r] == ""} {
		toplevel .r
		wm title .r RAYTRACE
		frame .r.top
		pack .r.top -side top -fill x

#Create the table in a separate procedure
		set tframe .r.top.tframe
		frame $tframe
		pack $tframe -side left -fill x -expand 1
		tableCreate $tframe
#Label is intended for ppm versions of diffraction patterns, but I don't
#really use it.  We'll create and pack but not allow expand.
#		label .r.top.label
#		pack .r.top.label -side left -fill both
		frame .r.bottom
		canvas .r.bottom.a -bg black -width 400 -height 400
		canvas .r.bottom.b -bg black -width 400 -height 400
		pack .r.bottom -side top -fill both -expand 1

		pack .r.bottom.a -side left -fill both -expand 1
		pack .r.bottom.b -side left -fill both -expand 1
		}
	   if {[info command .r.bottom.$subwin] == ""} {
		error "Subwindow $subwin does not exist"
		}
	   pgBegin .r.bottom.$subwin/TK $nx $ny
	   pgAsk 0
	   pgPage
	   }
	}
   pgSch 1
   pgSci 1
   return
   }

#################################################################
proc tableCreate {top} {
   global widget

#Need two frames, one to hold table and right scrollbar, another to
#hold bottom scrollbar and spacer
   frame $top.uframe -bg grey30
   set table $top.uframe.table
   table $table -titlerows 4 -titlecols 2 \
	-variable tktable -bg white -cols 34 \
	-colwidth 16 -font "Mono 10"
   $table tag configure title -bg grey -fg black \
	-font "Helvetica 10 bold" -relief raised
   $table tag row data 1
   $table tag raise data
   $table tag configure data -font "Mono 10" -bg white

#Cache table name
   set widget(table) $table
   set rscroll $top.uframe.rscroll
   scrollbar $rscroll -command "$table yview" -troughcolor grey30 -jump 0 \
	-width 15
   $table configure -yscrollcommand "$rscroll set"
   frame $top.lframe -bg grey30
   set cscroll $top.lframe.cscroll
   scrollbar $cscroll -command "$table xview" -troughcolor grey30 -orient \
	horizontal -jump 0 -width 15
   $table configure -xscrollcommand "$cscroll set"
   pack $top.uframe -side top -fill x -expand 1
   pack $table -side left -fill x -expand 1
   pack $rscroll -side right -fill y
   pack $top.lframe -side top -fill x
   pack $cscroll -side left -fill x -expand 1
   return
   }

###################################################################
proc plotEnd {} {
   pgEnd
   if {[info command .r] != ""} {destroy .r}
   return
   }

###################################################################
#Quick way to set the plotting scale
proc plotScale {scale xoff yoff} {
   global XMIN XMAX YMIN YMAX
   if {$scale == 0} {
	if {[info exists XMIN]} {unset XMIN}
	if {[info exists XMAX]} {unset XMAX}
	if {[info exists YMIN]} {unset YMIN}
	if {[info exists YMAX]} {unset YMAX}
	return
	}
   set XMIN [expr -1.*$scale + $xoff]
   set XMAX [expr $scale + $xoff]
   set YMIN [expr -1.*$scale + $yoff]
   set YMAX [expr $scale + $yoff]
   return
   }

######################################################################
proc scaleSet {{val 0} {xoff 0} {yoff 0}} {
   plotScale $val $xoff $yoff
   return
   }

#####################################################################
#Easier to type ...
proc opticList {optic icolor} {
   opticDisplay $optic $icolor
   return
   }

#####################################################################
#Display an optical design in a tkTable.

proc opticDisplay {optic icolor} {
   global tktable widget
   if {[info exists tktable]} {unset tktable}
   opticInfo $optic $icolor

#Focal plane info
   set titles "wave xpos ypos xrad  yrad scale rot dist map weight fl focus \
	exit xmag entrance emag"
   set tktable(0,0) Filter
   set tktable(1,0) $icolor
   loop i 0 [llength $titles] {
	set itable [expr $i+1]
	set tktable(0,$itable) [string toupper [lindex $titles $i]]
	}
   loop i 0 [llength $titles] {
	set itable [expr $i+1]
	set title [lindex $titles $i]
	set $title [showFocal $optic $icolor $title]
	set tktable(1,$itable) [set $title]
	}

#Display object zpos.  I now entertain the fact that this might be at a
#finite distance
   set itable [expr [llength $titles]+1]
   set tktable(0,$itable) "Object ZPOS"
   set tktable(1,$itable) [format %.3g [showSurf $optic 0 z]]

#Now display the surface data
   set nsurf [handleShow $optic.nsurf]

#Construct list of surfaces for this icolor
   set surflist ""
   for {set isurf 1} {$isurf <= $nsurf} {incr isurf} {
	set id [surfId $optic $isurf]
	set index [showIndex $optic $id $icolor]
	if {$index == 0} continue
	lappend surflist $isurf
	}
   
   set nrow [expr [llength $surflist]+4]
   $widget(table) configure -rows $nrow

   set jtable 3
   set inext 0
   foreach isurf $surflist {
	set id [surfId $optic $isurf]
	set index [showIndex $optic $id $icolor]
	incr jtable
	set id [surfId $optic $isurf]
	set tktable($jtable,0) $id
	set name [showName $optic $id]
	set tktable($jtable,1) $name
	set titles "curv ccon z x y phi theta a2 a4 a6 a8 a10 astig aphi \
	   instop outstop \
	   stoptype x1 y1 x2 y2 x3 y3 x4 y4 reflect lines order"
	loop i 0 [llength $titles] {
	   set title [lindex $titles $i]
	   set $title [showSurf $optic $id $title]
	   }

#Add extra entries.
#Radius of curvature
	lappend titles index glass
	if {$curv == 0.} {
	   set radius 0.
	} else {
	   set radius [format %.5f [expr 1./$curv]]
	   }
	set ix [expr [lsearch $titles curv] + 1]
	set titles [linsert $titles $ix radius]

#Compute thickness.  I assume things align along z axis.
	incr inext
	if {$inext < [llength $surflist]} {
	   set z1 [showSurf $optic [surfId $optic \
		[lindex $surflist $inext]] z]
	   set thick [format %.5f [expr $z1 - $z]]
	} else {
	   set thick 0
	   }
	set ix [expr [lsearch $titles z] + 1]
	set titles [linsert $titles $ix thick]

#Glass type
	lappend titles index glass
	set glass [showGlass $optic $id]

#Set titles first time through only
	if {$jtable == 4} {
	   set tktable(3,0) "SURF ID"
	   set tktable(3,1) NAME
	   loop i 0 [llength $titles] {
		set itable [expr $i+2]
		set tktable(3,$itable) [string toupper [lindex $titles $i]]
		}
	   }

	loop i 0 [llength $titles] {
	   set itable [expr $i+2]
	   set title [lindex $titles $i]
	   set tktable($jtable,$itable) [set $title]
	   }
	}
   wm title .r "RAYTRACE: [showDesign $optic], Filter $icolor"
   return
   }

######################################################################
#Difference input structure with initial one from setup

proc optdiff {std hndl} {

echo Difference are of the sense (first - second)
echo Looping through surfaces
   set diff [opticNew]
   opticCopy $std $diff
   for {set i 1} {$i <= [handleShow $hndl.nsurf]} {incr i} {
	loop j 1 13 {
	   handleSet $diff.optic<$i>->param<$j> \
		[expr { [handleShow $hndl.optic<$i>->param<$j>] - \
		[handleShow $std.optic<$i>->param<$j>]}]
	   }
	}
#Now loop through colors
echo Looping through colors
   for {set j 1} {$j <= [handleShow $hndl.ncolor]} {incr j} {
	loop i 1 7 {
	   handleSet $diff.fplane<$i>->param<$j> \
		[expr { [handleShow $hndl.fplane<$i>->param<$j>] - \
		[handleShow $std.fplane<$i>->param<$j>]}]
	   }
	}
   return $diff
   }

####################################################################
#In a plotting window, fetch the position of the pointer, convert to
#world coordinates.  Return an empty list if pointer is NOT in the
#window.

proc pointerCoord {relx rely window} {

#Recompute relx, rely
   set list [winfo pointerxy $window]
   set absx [lindex $list 0]
   set absy [lindex $list 1]
   set rootx [winfo rootx $window]
   set rooty [winfo rooty $window]
#echo rootx $rootx root y $rooty
#echo absx $absx absy $absy
#   set relx [expr $absx-$rootx]
#   set rely [expr $absy-$rooty]

#Get the window size and the viewport coords in device land
#I want widget-based width and height, which remain intact even if window is
#resized.

#These are cached by the tk driver in the "page" size parameters.
#Note that canvas width, height are not used by the driver.
#Instead, it picks up the window width, height, which chage whenever the
#window is resized.
   set line [pgQpsz]

   set width [lindex $line 2]
   set height [lindex $line 3]

#Get viewport limits (expressed as fraction of width and height)
   set list [pgQvp]
   set dxmin [lindex $list 0]
   set dxmax [lindex $list 1]
   set dymin [lindex $list 2]
   set dymax [lindex $list 3]

#Get viewport world coordinate limits
   set list [pgQwin]
   set wxmin [lindex $list 0]
   set wxmax [lindex $list 1]
   set wymin [lindex $list 2]
   set wymax [lindex $list 3]

#Finally ... the desired world coords.
   set wx [expr $wxmin + ($relx-$width*$dxmin)/($width*($dxmax-$dxmin)) * \
	($wxmax - $wxmin)]

#Need to flip y, since top is origin in X land but bottom is origin in plot
#land.
   set wy [expr $wymin + (($height-$rely)-$height*$dymin) / \
	($height*($dymax-$dymin)) * ($wymax - $wymin)]
   return [list $wx $wy]
   }

######################################################################
proc pgCurs {} {
   global _PGCURS
   set dev [pgQinf DEVICE]
   if {$dev != "TK"} return
   set window [pgQinf FILE]
   bind $window <Button-1> {global _PGCURS; set _PGCURS [list %x %y %W]}
   vwait _PGCURS
   bind $window <Button-1> ""
   return $_PGCURS
   }

######################################################################
#Query pointer and return world coords

proc pgCoord {} {
   set list [pgCurs]
   set x [lindex $list 0]
   set y [lindex $list 1]
   set win [lindex $list 2]
   set list [pointerCoord $x $y $win]
   return $list
   }

