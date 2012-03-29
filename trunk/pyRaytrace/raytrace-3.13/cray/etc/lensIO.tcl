#Restructure the cray data files.  The old FORTRAN-oriented style is
#woefully out of date.

#What is the format?  I will adopt the Yanny parameter file format.
#Will eventually want a way to add arbitrary comments

proc lensWrite {optic file args} {
   global _history
   set ext [file extension $file]
   if {$ext == ""} {
	set file $file.lns
	set ext .lns
	}
   set fid [open $file w]

#First, some basic information
#Ray Pattern
   set nstep [exprGet $optic.diagram->nstep]
   set ntype [exprGet $optic.diagram->ntype]

#Write out only the tail - don't normally want the directory info.
   puts $fid "file [file tail $file]"
   puts $fid ""
   puts $fid "#Ray Pattern information"
   puts $fid "nstep $nstep"
   puts $fid "ntype $ntype"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   char comment<80>;"
   puts $fid "   } HISTORY;"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   char param<10>;"
   puts $fid "   double val;"
   puts $fid "   } TEL;"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   int filter;"
   puts $fid "   char param<10>;"
   puts $fid "   double val;"
   puts $fid "   } FOCAL;"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   int filter;"
   puts $fid "   int mode;"
   puts $fid "   } PUPIL;"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   double surfid;"
   puts $fid "   char param<10>;"
   puts $fid "   double val;"
   puts $fid "   } SURF;"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   double surfid;"
   puts $fid "   int filter;"
   puts $fid "   double val;"
   puts $fid "   } INDEX;"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   double surfid;"
   puts $fid "   char glass<80>;"
   puts $fid "   } GLASS;"
   puts $fid ""
   puts $fid "typedef struct {"
   puts $fid "   double surfid;"
   puts $fid "   char comment<80>;"
   puts $fid "   } COM;"
   puts $fid ""

#Write out history records
#These are input via the "args" convention
   if {[string length $args] > 0} {
	foreach arg $args {
	   historyAdd $optic $arg
	   }
	}
   if {[info exists _history($optic)]} {
	loop i 0 $_history($optic) {
	   puts $fid [list history $_history($optic,$i)]
	   }
	puts $fid ""
	}

#Write out non-zero telescope parameters
   foreach param "diam fr1 fr2 back finner fiberdiam" {
	set $param [exprGet $optic.tel->$param]
	if {$param != 0} {puts $fid [format "tel %10s %15.6g" \
	   $param [set $param]]}
	}
   puts $fid ""
   puts $fid "#Focal plane parameters"
   puts $fid [format "#     %6s %10s %15s" Filter Parameter Value]

#Focal plane parameters
   set ncolor [exprGet $optic.ncolor]
   for {set i 1} {$i <= $ncolor} {incr i} {
	puts $fid ""
	foreach param "xoff yoff xsize ysize scale wave dist rot weight map" {
	   set val [showFocal $optic $i $param]
	   if {$val != 0} {
		puts $fid [format "focal %6d %10s %15.6g" $i $param $val]
		}
	   }
	}
   puts $fid ""

#Entrance pupil parameters
   set ncolor [exprGet $optic.ncolor]
   for {set i 1} {$i <= $ncolor} {incr i} {
	set mode [exprGet $optic.pupil<$i>->mode]
	if {$mode != 0} {
	   puts $fid [format "pupil %6d %10d" $i $mode]
	   }
	}
   puts $fid ""

#Surface parameters
   puts $fid "Surface parameters"
   set nsurf [exprGet $optic.nsurf]

   puts $fid [format "#         %-6s %10s %15s" SurfId Parameter Value]
   for {set i 0} {$i <= $nsurf} {incr i} {
	puts $fid ""
	set surfid [surfId $optic $i]

#Comment (also used to cache a name for the surface - e.g., PRIMARY).
	set format "com       %-6s %s"
	set com [string trim [showComment $optic $surfid]]
	if {[string length $com] > 0} {
	   puts $fid [format $format $surfid $com]
	   }

#Main parameters
	set format "surf      %-6s %10s %15.9g"

#I need to write out at least one parameter for each surface here because of
#the fact that I now break out the refraction indices into a separate
#data structure.  I must ensure that the surface is defined.
	set ndef 0
	foreach param "curv ccon x y z theta phi a2 a4 a6 a8 a10 \
		a1 a3 a5 a7 a9 a11 a13 \
		astig aphi instop outstop stoptype order lines blaze \
		 x1 y1 x2 y2 x3 y3 x4 y4 reflect" {
	   set val [showSurf $optic $surfid $param]
	   if {$val != 0.} {
		puts $fid [format $format $surfid $param $val]
		incr ndef
		}
	   }

#Default - write out curvature 0 if no other params were written.
	if {$ndef == 0} {puts $fid [format $format $surfid curv 0]}
	set format "index     %-6s %10d %15.9g"
	for {set icolor 1} {$icolor <= $ncolor} {incr icolor} {
	   set val [showIndex $optic $surfid $icolor]
	   if {$val != 0} {
		puts $fid [format $format $surfid $icolor $val]
	 	}
	   }

#Glass type
	set format "glass     %-6s %20s"
	set glass [showGlass $optic $surfid]
	if {$glass != ""} {
	   puts $fid [format $format $surfid $glass]
	   }
	}
   close $fid
   return
   }

###################################################################
#First, a procedure to find a lens data file.  This will also be used by
#opticRead, hopefully

proc dataFileFind {file ext} {
   global env

#Strip leading . if supplied for extention
   set ext [string trimleft $ext .]

#First, check file existing as supplied veratim
   if {![file exists $file]} {
	set fext [file extension $file]
	if {$fext == ""} {set file $file.$ext}
	}
   if {![file exists $file]} {

#Search through CRAY data directories.  I will use the POSIX find command.
	set result [exec find $env(CRAY_DIR)/data -name $file]
	set result [split $result \n]
	set file [lindex $result 0]	
	}
   if {![file exists $file]} then {error "File $file does not exist"}
   return $file
   }

####################################################################
#Read and write design name

proc showDesign {optic} {

#Unnecessary action - braces are placed around a name when exprGet
#references a character array.  Thus, the lindex function.

   set design [lindex [exprGet $optic.name] 0]

#Strip off file extension
   set design [file tail [file root $design]]
   return $design
   }

###################################################################
#Set overall name of a design.
proc setDesign {optic name} {
   handleSet $optic.name $name
   return
   }

###################################################################
#Set fiber diameter - used when constraining incidence angles.
#Diameter is in mm.

proc setFiber {optic diam} {
   handleSet $optic.tel->fiberdiam $diam
   return
   }

###################################################################
#Show fiber diameter - used when constraining incidence angles.
#Diameter is in mm.

proc showFiber {optic} {
   set diam [exprGet $optic.tel->fiberdiam]
   return $diam
   }

####################################################################
#Read back a lens file.

#TCL replacement for readin.

proc lensRead {file} {
   global env
   global _optic
   global _history
   set file [dataFileFind $file lns]
   set extension [file extension $file]
   set hndl [opticNew]
   set fid [open $file]
   handleSet $hndl.name [file tail $file]

#Read in parameter array
   paramRead $file hdr table

#typedefs should be tel, focal, surf, index, and com
#Do some error checking
   foreach type "tel focal surf index com history" {
	if {[lsearch $table() $type] < 0} {
	   error "Data type $type not found in $file"
	   }
	}
   loop i 0 $table(tel,counter) {
	set param $table(tel,$i,param)
	set val $table(tel,$i,val)
	if {[catch {handleSet $hndl.tel->$param $val}]} {
	   error "Bad param $param or val $val for datatype TEL"
	   }
	}

#I need fl1 in opticin
   handleSet $hndl.tel->fl1 [expr [exprGet $hndl.tel->diam] * [exprGet \
	$hndl.tel->fr1]]

#Final focal length - not really needed.
   handleSet $hndl.tel->f3 [expr [exprGet $hndl.tel->diam] * [exprGet \
	$hndl.tel->fr2]]

   opticzero $hndl

#Focal plane
   loop i 0 $table(focal,counter) {
	set filter $table(focal,$i,filter)
	set param $table(focal,$i,param)
	set val $table(focal,$i,val)
	if {[catch {setFocal $hndl $filter $param $val}]} {
	   error "Bad param $param or val $val for datatype FOCAL"
	   }
	}

#Entrance pupil
   if {[info exists table(pupil,counter)]} {
	loop i 0 $table(pupil,counter) {
	   set filter $table(pupil,$i,filter)
	   set mode $table(pupil,$i,mode)
	   if {[catch {handleSet $hndl.pupil<$filter>->mode $mode}]} {
		error "Bad mode $mode for datatype PUPIL"
		}
	   }
	}

#Surfaces
   handleSet $hndl.nsurf 0
   set lastsurf 0
   set newsurf 0
   loop i 0 $table(surf,counter) {
	set surfid $table(surf,$i,surfid)
	set param $table(surf,$i,param)
	set val $table(surf,$i,val)
	if {$surfid != $lastsurf} {
	   set lastsurf $surfid
	   set newsurf 1
	   }
	set list [surfTrans $surfid $param]
	set iparam [lindex $list 1]
	if {[catch {setparam $hndl $surfid $iparam $val}]} {
	   error "Bad surface: surfid $surfid param $param value $val"
	   }
	if {$newsurf == 1} {
	   set nsurf [exprGet $hndl.nsurf]
	   set newsurf 0
	   }
	}

#Refraction indices
   loop i 0 $table(index,counter) {
	set surfid $table(index,$i,surfid)
	set filter $table(index,$i,filter)
	set val $table(index,$i,val)
	if {[catch {setIndex $hndl $surfid $filter $val}]} {
	   error "Bad refraction index: surfid $surfid filter $filter val $val"
	   }
	}

#Glass (not in all files)
   if {[info exists table(glass,counter)]} {
	   loop i 0 $table(glass,counter) {
		set surfid $table(glass,$i,surfid)
		set glass $table(glass,$i,glass)
		setGlass $hndl $surfid $glass
		}
	}

#Comments
   loop i 0 $table(com,counter) {
	set surfid $table(com,$i,surfid)
	set comment $table(com,$i,comment)
	set isurf [surfIndex $hndl $surfid]
	if {[catch {handleSet $hndl.optic<$isurf>->comment \
	   [string range $comment 0 78]}]} {
	   error "Bad comment: surfid $surfid comment $comment"
	   }
	}

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

#If this is .len format, convert to my new index convention.
   if {$extension == ".len"} {
	indexConvert $hndl
	}
#Define ray pattern
   set names [keylkeys hdr]
   if {[lsearch $names nstep] >= 0 && [lsearch $names ntype] >= 0} {
	set nstep [keylget hdr nstep]
	set ntype [keylget hdr ntype]
	rayPattern $hndl $nstep $ntype
   } else {
	rayPattern $hndl 6 1
	}

#We seem to need to call opticinc before calling pupilInfo.  opticinc calls
#opticInfo.
   opticinc $hndl 1

#Compute pupil info if pupil mode is 1.
   set ncolor [exprGet $hndl.ncolor]
   for {set i 1} {$i <= $ncolor} {incr i} {
	if {[exprGet $hndl.pupil<$i>->mode] > 0} {
	   pupilInfo $hndl $i
	   }
	}

#Display first color
   if {[info command .r] != ""} {
	opticDisplay $hndl 1
	}

#Cache comments.  Track by handle name.  This could lead to problems if
#I run opticCopy - will need to update it to understand _history
   historyInit $hndl
   loop i 0 $table(history,counter) {
	historyAdd $hndl $table(history,$i,comment)
	}

   return $hndl
   }

###########################################################################
#For symmetry, lensDel is like opticDel
proc lensDel {hndl} {
   opticDel $hndl
   return
   }

###########################################################################

#History commands

proc historyInit {hndl} {
   global _history
   if {[info exists _history($hndl)]} {
	loop i 0 $_history($hndl) {
	   unset _history($hndl,$i)
	   }
	unset _history($hndl)
	}
   return
   }

######################################################################
proc historyAdd {hndl msg} {
   global _history
   if {![info exists _history($hndl)]} {
	set _history($hndl) 0
	}
   set _history($hndl,$_history($hndl)) [string range $msg 0 64]
   incr _history($hndl)
   return
   }

########################################################################
proc historyList {hndl} {
   global _history
   if {![info exists _history($hndl)]} return
   loop i 0 $_history($hndl) {
	echo $i: $_history($hndl,$i)
	}
   return
   }

#####################################################################
proc historyDel {hndl n} {
   global _history
   if {![info exists _history($hndl)]} return
   if {$n >= $_history($hndl)} return
   if {$n < 0} return
   loop i $n $_history($hndl) {
	set _history($hndl,[expr $i-1]) $_history($hndl,$i)
	}
   incr _history($hndl) -1
   unset _history($hndl,$_history($hndl))
   return
   }

#########################################################################
#Convert refractive indices to new convention

proc indexConvert {hndl} {
   set ncolor [exprGet $hndl.ncolor]
   set nsurf [exprGet $hndl.nsurf]
   for {set ifil 0} {$ifil <= $ncolor} {incr ifil} {
	set flip 0
	set index1 0
	for {set isurf 0} {$isurf <= $nsurf} {incr isurf} {
	   set fsurf [surfId $hndl $isurf]
	   set index [showIndex $hndl $fsurf $ifil]
	   set glass [showGlass $hndl $fsurf]
	   if {$index == 0} continue
	   set flip [expr $index1*$index]
	   set index1 $index
	   if {$index == -1 && $glass == "air"} {
		set index 1.
		}

#If index is negative and not a mirror
	   if {$index == -1} {set glass mirror}
	   if {$index < 0 && $glass != "mirror"} {
		set index [expr abs($index)]
		}
	   if {$flip < 0 && $index == 1} {
		set index -1.
		set glass mirror
		}
	   setIndex $hndl $fsurf $ifil $index

#Let's hope I get the same glass for all filters!
	   setGlass $hndl $fsurf $glass
	   }
	}
   return
   }

#######################################################################
#Convert a design from .len to .lns format

proc lensConvert {file} {
   set ext [file extension $file]
   set root [file root $file]
   if {$ext == ""} {
	set ext .len
	}
   set hndl [lensRead $root$ext]
   lensWrite $hndl $root.lns
   opticDel $hndl
   return
   }
