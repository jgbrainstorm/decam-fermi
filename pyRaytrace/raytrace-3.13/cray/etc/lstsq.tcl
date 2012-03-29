#TCL proc for running least squares.
#Procedures:
#	lstsq <optic> [lscale]
#	fastLstsqInit (must have set flags first)
#	fastLstsq <optic>
#
#	flagList
#
#	clearSpot
#	addSpot  <x> <y> <weight> <lscale>
#	listSpot
#
#Known problem: If one of the lstsq routines bombs, we leave an LSTSQ struct
#behind.

proc lstsq {optic {lscale ""}} {
   global _optic

#Set flags in optic structure
   flagSet $optic

#Cache flags
   cacheFlags

#Cache current optic structure
   if {![info exists _optic(opticache)] || \
	   [catch {handleGet $_optic(opticache)}]} {
	set _optic(opticache) [opticNew]
	}
   opticCopy $optic $_optic(opticache)
   flush stdout

   set _optic(lscale) $lscale
   spotQuery $lscale
   if {![info exists _optic(spots)]} {
        error "Please use spotQuery to enter spot positions!"
        }

#Initialize
   set lstsq [genericNew LSTSQ]

#Replay actual spots
   spotSet $lstsq

#Backward compatibility:
#   If lscale is supplied, set all spot flags to this value.
   if {$lscale != ""} {
	loop i 0 [exprGet $lstsq.nspotdef] {
	   handleSet $lstsq.lscale<$i> $lscale
	   }
	}
   lstsq1 $optic $lstsq
   lstsq1a $lstsq
   lstsq1b $lstsq

#Echo flags and cache current parameter values.
   if {[verbose]} {flagList $optic}
   set niter 0

#Automatically run 1 iteration to set spot centers.  This will not be
#counted in the iteration count.

   lstsq2 $optic $lstsq

#The following is a bit of a hack
   handleSet $lstsq.fractmax 1.
   lstsq3 $optic $lstsq

#Number of iterations to run on this pass
   set totiter 1
   while {1} {

#Now print out residuals before running the next iteration.

#The goal is to parallelize the next statement, which is where all the CPU
#goes.

	if {[verbose]} {
	   puts stdout "Iteration $niter"
	   catch {rmsList $optic $lstsq}
	   }

#Decrement total number of iterations this pass.  If 0, ask if we want to
#continue
	incr totiter -1
	if {$totiter <= 0} {
	   while {1} {
		set yn [getLine "Iterate ? "]
		set yn [string index [string tolower $yn] 0]
		if {$yn == "" || $yn == "y" || $yn == "n"} break
		}
	   if {$yn != "y"} break
	   while {1} {
		set line [getLine "niter, fractmax = "]
		set totiter [lindex $line 0]
		set fractmax [lindex $line 1]
		if {$fractmax == ""} {set fractmax 1.}
		if {$totiter == ""} {set totiter 1}
		if {[catch {format %d $totiter}]} continue
		break
		}
	   }
	incr niter

	lstsq2 $optic $lstsq

#The following is a bit of a hack
	handleSet $lstsq.fractmax $fractmax
	lstsq3 $optic $lstsq
	if {[verbose]} {updateList $optic}

#Check for convergence
	set fract [exprGet $lstsq.fract]
	if {$fract == 0} break
	}
   lstsq4 $lstsq

#More caching
   foreach var "rmsx rmsy nxy rmscx rmscy ncxy rmsi ni fractmax" {
	set $var [exprGet $lstsq.$var]
	}
   set _optic(niter) $niter
   set _optic(rmsx) $rmsx
   set _optic(rmsy) $rmsy
   set _optic(nxy) $nxy
   set _optic(rmscx) $rmscx
   set _optic(rmscy) $rmscy
   set _optic(ncxy) $ncxy
   set _optic(rmsi) $rmsi
   set _optic(ni) $ni
   set _optic(fractmax) $fractmax
   genericDel $lstsq
   return
   }

##############################################################
#Replay spots from cache and update lstsq structure.

proc spotSet {lstsq} {
   global _optic
   if {![info exists _optic(spots)]} {
	spotQuery 0
	}
   spotInit $lstsq
   foreach quad $_optic(spots) {
	eval spotAdd $lstsq $quad
       if {[verbose]} {echo Using spot $quad}
	}
   return
   }

##############################################################
#List spots in cache

proc listSpot {} {
   spotList
   return
   }

proc spotList {} {
   global _optic
   if {![info exists _optic(spots)]} {return}
   set format "%8s %8s %8s %8s"
   puts stdout [format $format xfract yfract weight lscale]
   foreach quad $_optic(spots) {
	puts stdout [eval format {$format} $quad]
	}
   return
   }

##############################################################
#Default procedure to input and cache spot positions.  Support grid,
#radial layouts.

#If clip = yes, we limits spots to the unit circle, even if layout is
#otherwise a grid.

proc spotQuery {lscale {clip no}} {
   global _optic

   if {$lscale == ""} {set lscale 0}

#Allow flexible input for "clip"
   set clip [string tolower $clip]
   if {$clip == "clip"} {set clip yes}

   set line [getLine "XMIN, XMAX: "]
   if {$line == ""} return

#Cache all parameters
   set xmin [lindex $line 0]
   set xmax [lindex $line 1]
   set _optic(xmin) $xmin
   set _optic(xmax) $xmax

   flush stdout
   set line [getLine "YMIN, YMAX: "]
   if {$line == ""} return

   set ymin [lindex $line 0]
   set ymax [lindex $line 1]
   set _optic(ymin) $ymin
   set _optic(ymax) $ymax

   set nzero 0
   foreach var "xmin xmax ymin ymax" {
	if {[set $var] == 0} {incr nzero}
	}
   set nspot [expr ($xmax-$xmin+1)*($ymax-$ymin+1)]
   if {$nspot == 0} {set nspot 1}

#Cache spot positions
   set _optic(spots) ""

#Weighting: In old scheme, I needed to have info on field shape.  I used
#uniform weight if rectangular, increasing weight with radius if circular.
#This is OK in some cases, not others.  Since I don't have field shape
#available yet, I will use a different scheme - if two or more out of xmin,
#xmax, #ymin, ymax are 0, I assume circular field, else rectangular.

   for {set x $xmin} {$x <= $xmax} {incr x} {
	for {set y $ymin} {$y <= $ymax} {incr y} {
	   set xfract [expr (1.*$x)/max(1., max(abs($xmin),abs($xmax)))]
	   set yfract [expr (1.*$y)/max(1., max(abs($ymin),abs($ymax)))]
	   if {$clip == "yes" && $xfract*$xfract + $yfract*$yfract > 1} \
		continue
	   if {$nzero < 2} {
		set weight 1.
	   } else {
		set weight [expr sqrt($xfract*$xfract + $yfract*$yfract)]
		if {$x == 0 && $y == 0} {
		   set weight [expr 1./$nspot]
		   }
		}
	   lappend _optic(spots) [list $xfract $yfract $weight $lscale]
	   if {[verbose]} {
		echo Spot at [format %.2f $xfract] [format %.2f $yfract] \
		 weight [format %.2f $weight]
		}
	   }
	}
   return
   }

##################################################################
#Clear spots in cache

proc clearSpot {} {
   global _optic
   if {[info exists _optic(spots)]} {
	set _optic(spots) ""	
	}
   return
   }

##################################################################
#Add a new spot to the cache.
proc addSpot {x y {weight 1} {lscale 0}} {
   global _optic
   if {![info exists _optic(spots)]} {set _optic(spots) ""}
   lappend _optic(spots) [list $x $y $weight $lscale]
   return
   }

##################################################################
#Echo flags
#In principle this is all in the lstsq structure, and I am just reinventing
#the wheel here.

proc flagList {optic} {
   global _optic
   global solve
   if {[info exists _optic(flags)]} {
	set flags $_optic(flags)
   } elseif {[info exists _optic(flagcache)]} {
	set flags $_optic(flagcache)
   } else {
	return
	}
   if {[info exists solve]} {unset solve}
   set surflist ""
   set solve(0) ""

#Lists of link colors and params
   set linkcolor ""
   set linkparam ""
   foreach combo $flags {
	set surfid [lindex $combo 0]
	set isurf [surfIndex $optic $surfid]
	set param [lindex $combo 1]
	set flag [lindex $combo 2]
	if {$isurf !=0} {
	   if {[lsearch $surflist $isurf] < 0} {
		lappend surflist $isurf
		set solve($isurf) ""
		}
	   if {[lsearch $solve($isurf) $param] < 0} {
		lappend solve($isurf) $param
		set solve($isurf,surfid) $surfid
		}
	   set solve($isurf,$param,flag) $flag
	   set solve($isurf,$param,val) [showparam $optic $surfid $param]
	   if {$flag > 1 && [lsearch $linkparam $flag] < 0} {
		lappend linkparam $flag
		}
	} else {
	   if {[lsearch $solve(0) $param] < 0} {
		lappend solve(0) $param
		set solve(0,surfid) $surfid
		}
	   set solve(0,$param,flag) $flag
	   if {$flag > 1 && [lsearch $linkcolor $flag] < 0} {
		lappend linkcolor $flag
		}
	   }
	}
   if {[llength $surflist] > 1} {
	set surflist [lsort -integer $surflist]
	}
   if {[llength $linkcolor] > 1} {
	set linkcolor [lsort -integer $linkcolor]
	}
   if {[llength $linkparam] > 1} {
	set linkparam [lsort -integer $linkparam]
	}
   if {[llength $solve(0)] > 1} {
	set solve(0) [lsort -integer $solve(0)]
	}
   foreach isurf $surflist {
	if {[llength $solve($isurf)] > 1} {
	   set solve($isurf) [lsort -integer $solve($isurf)]
	   }
	}

#Now echo results
   echo
   foreach color $solve(0) {
	if {$solve(0,$color,flag) == 1} {
	   echo Adjusting color $color
	   }
	}
   foreach flag $linkcolor {
	set line "Linking (flag=$flag) colors: "
	foreach color $solve(0) {
	   if {$solve(0,$color,flag) == $flag} {
		lappend line $color
		}
	   }
	echo $line
	}
   echo
   foreach isurf $surflist {
	foreach param $solve($isurf) {
	   if {$solve($isurf,$param,flag) ==1} {
		echo Solving surface $solve($isurf,surfid) param $param
		}
	   }
	}
   foreach flag $linkparam {
	echo "Linking (flag=$flag):"
	foreach isurf $surflist {
	   foreach param $solve($isurf) {
		if {$solve($isurf,$param,flag) == $flag} {
		   echo "        surface $solve($isurf,surfid) param $param"
		   }
		}
	   }
	}
   echo
   set solve(surflist) $surflist
   set solve(linkcolor) $linkcolor
   set solve(linkparam) $linkparam
   return
   }

##################################################################
#List rms.

proc rmsList {optic lstsq} {
   foreach var "rmsx rmsy nxy rmscx rmscy ncxy rmsi ni" {
	set $var [exprGet $lstsq.$var]
	}
   if {$nxy == 0} {
	echo No spots made it to the focal plane
	return
	}
   if {$nxy < 0} {
	echo Weird! nxy < 0
	return
	}
   if {$rmsx < 0. || $rmsy < 0.} {
	echo Gak! rmsx or rmsy is negative!
	return
	}
   set fract [exprGet $lstsq.fract]
   puts stdout [format "Fractional increment: %.3f" $fract]

#Gak!  I often have the wrong scale factor for filter 1 if I am working
#with other filters.  Oh well, another way to go wrong.
   set scale [expr abs([showScale $optic 1])]
   set rmsx [expr sqrt($rmsx/$nxy)*$scale]
   set rmsy [expr sqrt($rmsy/$nxy)*$scale]
   puts stdout [format \
"Image size residuals in arcsec: Rms(X) = %8.3f   Rms(Y) = %8.3f" \
	$rmsx $rmsy]
   if {$ncxy > 0} {
	if {$rmscx < 0. || $rmscy < 0.} {
	   echo Gak! rmscx or rmscy is negative!
	   return
	   }
	set rmscx [expr sqrt(1.*$rmscx/$ncxy)*$scale]
	set rmscy [expr sqrt(1.*$rmscy/$ncxy)*$scale]
	puts stdout [format \
   "Image position residuals in arcsec: Rms(X) = %8.3f   Rms(Y) = %8.3f" \
	   $rmscx $rmscy]
	}
   if {$ni > 0} {
	if {$rmsi < 0.} {
	   echo Gak! rmsi is negative!
	   return
	   }
	set rmsi [expr sqrt(1.*$rmsi/$ni)*57.3]
	puts stdout [format \
   "Incidence angle residuals in deg: Rms = %8.3f" $rmsi]
	}
   return
   }

##################################################################
#List updated parameters after least squares.
#In principle this is all in the lstsq structure, and I am just reinventing
#the wheel here.

proc updateList {optic} {
   global _optic
   global solve
   if {![info exists solve]} return
   set surflist $solve(surflist)

#Lists of link colors and params
   set linkcolor $solve(linkcolor)
   set linkparam $solve(linkparam)

#Now echo results
   echo
   foreach isurf $surflist {
	foreach param $solve($isurf) {
	   set oldval $solve($isurf,$param,val)
	   set val [showparam $optic $solve($isurf,surfid) $param]
	   set solve($isurf,$param,val) $val
	   set diff [expr $val-$oldval]
	   
	   echo Surface $solve($isurf,surfid) param $param \
		old value [format %10.3e $oldval] \
		increment [format %10.3e $diff]
	   }
	}
   echo
   return
   }

##################################################################
#Normally fastLstsq is run after we have run one iteration of lstsq.
#However, in some cases we like to run it right from the beginning.
#I need to initialize some values that are normally cached.
#Assume that flags have already been set.

proc fastLstsqInit {{lscale ""}} {
   global _optic

#The following default values are normally what is wanted
   set _optic(lscale) $lscale
   set _optic(fractmax) 1
   set _optic(niter) 1

#Need to create a flags cache.  I assume that I have already initialized
#the flags using setColorFlag, etc.
#Cache flags
   if {[info exists _optic(flags)]} cacheFlags
   return
   }

##################################################################
#Repeat a least squares fit based on cached parameters from a run using
#lstsq above.

proc fastLstsq {optic} {
   global _optic

#Cache current optic structure
   if {![info exists _optic(opticache)] || \
	   [catch {handleGet $_optic(opticache)}]} {
	set _optic(opticache) [opticNew]
	}
   opticCopy $optic $_optic(opticache)

   resetFlags $optic

#Set flags in optic structure
   flagSet $optic

#Cache flags
   cacheFlags

   set lscale $_optic(lscale)
   set fractmax $_optic(fractmax)

#Initialize
   set lstsq [genericNew LSTSQ]

   spotSet $lstsq

#Backward compatibility:
#   If lscale is supplied, set all spot flags to this value.
   if {$lscale != ""} {
	loop i 0 [exprGet $lstsq.nspotdef] {
	   handleSet $lstsq.lscale<$i> $lscale
	   }
	}

   lstsq1 $optic $lstsq
   lstsq1a $lstsq
   lstsq1b $lstsq

#Automatically run 1 iteration to set spot centers.  This will not be
#counted in the iteration count.
   lstsq2 $optic $lstsq

#First iteration computes the spot centers
#The following is a bit of a hack
   handleSet $lstsq.fractmax 1.
   lstsq3 $optic $lstsq

   set niter $_optic(niter)

#There was a bug here in the old version: After running lstsq3, I did not
#have the latest rms computed.  In plstsq, I ran 1 extra iteration and bailed
#before the final update.  With the new lstsq3, I get the latest rms right
#away, so this bug is squashed.
   loop i 0 $niter {

#The goal is to parallelize the next statement, which is where all the CPU
#goes.
	lstsq2 $optic $lstsq

#The following is a bit of a hack
	handleSet $lstsq.fractmax $fractmax
	lstsq3 $optic $lstsq

#Check for convergence
	set fract [exprGet $lstsq.fract]
	if {$fract == 0} break
	}
   lstsq2 $optic $lstsq
   lstsq4 $lstsq

#More caching
   foreach var "rmsx rmsy nxy rmscx rmscy ncxy rmsi ni" {
	set $var [exprGet $lstsq.$var]
	}
   set _optic(rmsx) $rmsx
   set _optic(rmsy) $rmsy
   set _optic(nxy) $nxy
   set _optic(rmscx) $rmscx
   set _optic(rmscy) $rmscy
   set _optic(ncxy) $ncxy
   set _optic(rmsi) $rmsi
   set _optic(ni) $ni
   genericDel $lstsq
   return
   }
