#Recompute least squares fit based on cached variables.
#Figure of merit is the rms radius stdev.

proc fom {optic} {
   global _optic

#Reset flags
   if {![info exists _optic(niter)]} {
	error "No cached least squares"
	}
   fastLstsq $optic
   if {![info exists _optic(fomType)]} {
	set _optic(fomType) spot
	}
   if {$_optic(fomType) == "spot"} {
	foreach var "rmsx rmsy nxy rmscx rmscy ncxy" {
	   set $var $_optic($var)
	   }
	set fom 0.

#rmsx, etc are actually the sum of squares.
	if {$nxy > 0} {
	   set fom [expr sqrt(($rmsx + $rmsy)/$nxy)]
	   }
	if {$ncxy > 0} {
	   set fom [expr sqrt($fom*$fom + ($rmscx + $rmscy)/$ncxy)]
	   }
   } else {
	set fom 0.
	set n 0

#Which color?  Look through the flags, pick the first color found
	set color 0
	foreach set $_optic(flags) {
	   set surf [lindex $set 0]
	   if {$surf != 0} continue
	   set color [lindex $set 1]
	   break
	   }
	if {$color == 0} {
	   error "No color found in flags"
	   }
	foreach spot $_optic(spots) {
	   set xfract [lindex $spot 0]
	   set yfract [lindex $spot 1]
	   set weight [lindex $spot 2]
	   set xrad [showFocal $optic 1 xrad]
	   set yrad [showFocal $optic 1 yrad]
	   set xmm [expr abs($xrad)*$xfract]
	   set ymm [expr abs($yrad)*$yfract]
	   set rmsErr [waveFrontRms $optic $xmm $ymm $color]
	   set fom [expr $fom + $rmsErr * $rmsErr * $weight]
	   set n [expr $n + $weight]
	   }
	if {$n == 0} {error "No valid spots/weights"}
	set fom [expr sqrt($fom/$n)]
	set wave [showWave $optic $color]
	set fom [expr $fom*$wave*1.e3]
	}
   resetOptic $optic
   return $fom
   }

#######################################################################
#What type of figure-of-merit?  Allow either spot size or
#waveFrontError in nm.  

proc fomType {{type ""}} {
   global _optic
   set valid [list spot wave]
   if {$type == ""} {
	echo Types are [join $valid ", "]
	return
	}
   if {[lsearch $valid $type] >= 0} {
	set _optic(fomType) $type
   } else {
	error "Invalid type \"$type\""
	}
   return
   }

#######################################################################
#Increment a parameter and see what happens to the figure of merit.
#fract is used to scale the value of the increment
#We allow for the possibility of incrementing parameters on multiple
#surfaces.  We do not yet allow for the possibility of incrementing
#multiple parameters as well - not clear if this would be useful.
#Well, it is (for refraction index)

proc dither {optic surfs params fract} {

#If param is an integer, treat it as a refraction index.
   set fom0 [fom $optic]
   set optic1 [opticNew]
   opticCopy $optic $optic1
   if {[ctype digit [lindex $params 0]]} {
	set showCmd showIndex
	set setCmd setIndex
	set showIncCmd showIndexInc
   } else {
	set showCmd showSurf
	set setCmd setSurf
	set showIncCmd showSurfInc
	}

#Use increment from first surface.
   set inc [$showIncCmd $optic [lindex $surfs 0] [lindex $params 0]]
   foreach surf $surfs {
	foreach param $params {
	   set val0 [$showCmd $optic $surf $param]
	   set val1 [expr $val0 + $fract*$inc]
	   $setCmd $optic1 $surf $param $val1
	   }
	}

   set fom1 [fom $optic1]

   opticCopy $optic $optic1

   foreach surf $surfs {
	foreach param $params {
	   set val0 [$showCmd $optic $surf $param]
	   set val2 [expr $val0 - $fract*$inc]
	   $setCmd $optic1 $surf $param $val2
	   }
	}
   set fom2 [fom $optic1]

   opticDel $optic1
   set diff1 [expr $fom1-$fom0]
   set diff2 [expr $fom2-$fom0]
   return [list $fom0 $diff1 $diff2]
   }

#####################################################
#Find the increment that changes the figure of merit by either a fixed
#fraction or up to a fixed amount

proc ditherStudy {optic surfs params {fractinc 1} {fixed 0} {converge .02}} {

#First, base figure of merit
   set fom0 [fom $optic]

#What is my target?  Take the bigger of fract*fom and fixed
   set fomtarg [expr max($fom0*(1.+$fractinc), $fixed)]
   if {$fomtarg < $fom0} {
	error "Bad criteria for figure of merit study"
	}

#Express target as the increase in fom
   set target [expr $fomtarg - $fom0]
   set factor(0) 6
   set n 0

#Nloop is the maximum number of iterations.
   set nloop 10
   loop i 0 $nloop {
	set list [dither $optic $surfs $params $factor($i)]
	set diff1 [lindex $list 1]
	set diff2 [lindex $list 2]
	set avg [expr ($diff1+$diff2)/2.]

#Assume diff varies as square of factor
#New guess at factor
	if {$avg == 0} {
	   echo Average is 0.
	   break
	   }
	incr n
	if {$avg > 0.} {
	   set factor($n) [expr $factor($i)*sqrt($target/$avg)]
	} elseif {$n > 1} {
	   set factor($n) [expr ($factor([expr $n-1]) + $factor([expr $n-2])) \
		/2.]
	} else {
	   set factor($n) 0
	   break
	   }
	if {[verbose]} {echo FACTOR $factor($n)}
	if {abs(($factor($n) - $factor($i))/$factor($n)) \
	   <= $converge} break
	}
   loop i 0 [expr $n+1] {
	if {[verbose]} {echo $i $factor($i)}
	}
   if {[ctype digit [lindex $params 0]]} {
	set showIncCmd showIndexInc
   } else {
	set showIncCmd showSurfInc
	}
   set inc [$showIncCmd $optic [lindex $surfs 0] [lindex $params 0]]
   set tolerance [expr abs($factor($n)*$inc)]
   return $tolerance
   }

########################################################################
#Clear out all tolerancing parameters.  Not strictly necessary if I run
#tolSurf, tolParam, and tolSetup, but needed to clean out old cruft in case
#I want to run tolDynamic, tolStatic, tolTarget alone, and also helps
#prevent errors if some steps are omitted.

proc tolInit {} {
   global _optic
   global tolTable
   if {[info exists _optic(tolSurf)]} {unset _optic(tolSurf)}
   if {[info exists _optic(tolParam)]} {unset _optic(tolParam)}
   if {[info exists tolTable]} {unset tolTable}
   return
   }

########################################################################
#Define surfaces for tolerancing.  Input list is a list of surface ids to
#tolerance.

proc tolSurf {args} {
   global _optic
   set _optic(tolSurf) ""
   foreach item $args {
	lappend _optic(tolSurf) $item
	}
   return
   }

######################################################################
#Define parameters for tolerancing.  Input is a list of parameter names.
#Each name can optionally also be a list giving multiple parameters.
#In this case, testing will be done only for the first of these parameters,
#but the results will be applied to both parameters.  This is useful if
#there are multiple degrees of freedom (e.g., x and y offsets) that
#have equal impact on the optical performance.

proc tolParam {args} {
   global _optic
   set _optic(tolParam) ""
   foreach item $args {
	lappend _optic(tolParam) $item
	}
   return
   }

#########################################################################
#Combine surface and parameter combinations into tolTable.
#I do not support linking surfaces.  Linking of parameters can be done but
#only by making them a list of lists.  This is actually true for tolSetup
#and all my other tol.... commands.  E.G., I need to say
#   tolDynamic h0 4 [list "1 2 3"]


proc tolSetup {} {
   global _optic
   global tolTable
   if {![info exists _optic(tolSurf)]} {
	error "Must run tolSurf first!"
	}
   if {![info exists _optic(tolParam)]} {
	error "Must run tolParam first!"
	}
   if {[info exists tolTable]} {unset tolTable}
   set tolList ""

#No linked surfaces allowed (yet) in tolSetup, although I think I could add.
   foreach surf $_optic(tolSurf) {
	if {$surf <= 0} continue
	foreach param $_optic(tolParam) {
	   set weight [llength $param]
	   set param [lindex $param 0]
	   lappend tolList [list $surf $param]
	   set tolTable($surf,$param,weight) $weight

#Default tolerance type is "dynamic"
	   set tolTable($surf,$param,type) dynamic
	   set tolTable($surf,$param,tol) 0
	   set tolTable($surf,$param,rms) 0
	   }
	}
   set tolTable() $tolList
   return
   }

#########################################################################
#Set a surface/parameter combo to be static (fixed tolerance) rather
#than dynamic.

#Param can be a list of parameters - the first is toleranced and weighted
#by the number of params.  Value is the static tolerance value in units
#appropriate to the parameter.

proc tolStatic {surf param value} {
   global tolTable
   if {![info exists tolTable]} {
	set tolTable() ""
	}
   set weight [llength $param]
   set param [lindex $param 0]
   if {![info exists tolTable($surf,$param,type)]} {
	lappend tolTable() [list $surf $param]
	set tolTable($surf,$param,weight) $weight
	}
   set tolTable($surf,$param,type) static
   set tolTable($surf,$param,tol) $value
   set tolTable($surf,$param,rms) 0
   return
   }

#########################################################################
#Set a surface/parameter combo to be dynamic. (f.o.m. determined dynamically)

#Param can be a list of parameters - the first is toleranced and weighted
#by the number of params

proc tolDynamic {surf param} {
   global tolTable
   if {![info exists tolTable]} {
	set tolTable() ""
	}
   set weight [llength $param]
   set param [lindex $param 0]
   if {![info exists tolTable($surf,$param,type)]} {
	lappend tolTable() [list $surf $param]
	set tolTable($surf,$param,weight) $weight
	}
   set tolTable($surf,$param,type) dynamic
   set tolTable($surf,$param,tol) 0
   set tolTable($surf,$param,rms) 0
   return
   }

#########################################################################
#Set a surface/parameter combo to have fixed incremental f.o.m.

#Param can be a list of parameters - the first is toleranced and weighted
#by the number of params

#Input value (rms) is added in quadrature with fom0 to get final fom1.
#This is changed from my previous strategy, where I would input fom1-fom0.

proc tolTarget {surf param rms} {
   global tolTable
   if {![info exists tolTable]} {
	set tolTable() ""
	}
   set weight [llength $param]
   set param [lindex $param 0]
   if {![info exists tolTable($surf,$param,type)]} {
	lappend tolTable() [list $surf $param]
	set tolTable($surf,$param,weight) $weight
	}
   set tolTable($surf,$param,type) target
   set tolTable($surf,$param,tol) 0
   set tolTable($surf,$param,rms) $rms
   return
   }

#########################################################################
#Determine tolerance limits for all surfaces and parameters specified in
#the setup.  I input a fractional tolerance limit and an absolute tolerance
#limit on the figure of merit.  Either one (but not both) can be null,
#in which case they will not be used.  The algorithm is to distribute
#tolerances equally among all degrees of freedom.  Parameters that are
#used as compensators (e.g., focus) are not included.
#It is necessary to have runs one pass of least squares with the compensators
#set.  The flags, spot positions, and weights from the fit will be used
#for tolerancing runs.

#fractinc is the total target fractional increase in FOM.
#(Thus, fomlim is fractinc*fom0)
#fixed is, alternatively, the total target FOM.  TH
#(fomtarget = max(fom0 + fomlim, fixed) = max(fom0*(1+fractinc), fixed)

proc tolLimit {optic {fractinc 1} {rmsinc 0} {converge .01}} {
   global _optic
   global tolTable
   if {![info exists tolTable]} {
	error "Must run tolSetup first"
	}
   if {![info exists _optic(flagcache)]  && ![info exists _optic(flags)]} {
	error "Must run least squares first!"
	}

#Always run with verbosity = 0
   set verbosity [verbose]
   if {$verbosity != 0} {
	echo Setting verbosity to 0
	verbose 0
	}

#Use flags first, then flagcache
   if {[info exists _optic(flags)]} {
	set flags $_optic(flags)
   } else {
	set flags $_optic(flagcache)
	}
   set tolTable(name) [lindex [exprGet $optic.name] 0]

#Cache info used by least squares
   set tolTable(flags) $flags
   set tolTable(spots) $_optic(spots)
   set tolTable(lscale) $_optic(lscale)
   set tolTable(niter) $_optic(niter)


#List of lstsqs flags includes linking parameters.  For purposes here,
#remake list to have just surface id and parameter, but keep track of
#flag.
#
#I allow for multiple surface id's being toleranced as a unit.  This is
#useful when I have multiple configurations where the same physical surface
#is duplicated in different configurations.
   foreach combo $flags {
	set surf [lindex $combo 0]
	set param [lindex $combo 1]
	set flag($surf,$param) [lindex $combo 2]
	}
   set ndyn 0
   set nstat 0
   set ntarg 0
   set tolDyn ""
   set tolStat ""
   set tolTarg ""
   foreach combo $tolTable() {

#No tolerancing of focal plane parameters for now.  Also, no tolerancing
#of index of refraction, although this will eventually be needed.  I
#think I should enter a list of surfaces and indices in a separate proc
#for that case.
	set surfs [lindex $combo 0]
	set params [lindex $combo 1]
	set weight $tolTable($surfs,$params,weight)
	set type $tolTable($surfs,$params,type)

#Check on surfaces - don't tolerance focal plane surfaces, and don't
#tolerance surface/parameters that are compensation variables unless they
#are linked parameters.  If I have compensation variables that are linked,
#I will allow the tolerancing of any one of them.  If I have
#linked surfaces being toleranced, I require that at least one surface
#qualify.  If I have the same linked surfaces being toleranced and used as
#compensation variables, I will have trouble.  I could try to detect,
#but the logic gets messy.
	set iflag 0
	foreach surf $surfs {
	   if {$surf >= 0} {
		foreach param $params {
		   set list [surfTrans $surf $param]
		   set index [lindex $list 1]

#Check tolerances variable against lstsqs flags.  If variable is not
#flagged, or if it is but flag is != 1 (and thus it is linked), allow
#tolerancing to proceed.
		   if {![info exists flag($surf,$index)]} {
			set iflag 1
		   } else {
			if {$flag($surf,$index) != 1} {set iflag 1}
			}
		   }
		}
	   }
	if {$iflag == 0} continue
	if {[verbose]} {
	   echo Tolerancing surface $surfs parameter $params weight $weight
	   }
	if {$type == "dynamic"} {
	   lappend tolDyn [list $surfs $params]
	   set ndyn [expr $ndyn+$weight]
	} elseif {$type == "target"} {
	   lappend tolTarg [list $surfs $params]
	   set ntarg [expr $ntarg+$weight]
	} else {
	   lappend tolStat [list $surfs $params]
	   set nstat [expr $nstat+$weight]
	   }
	}

   if {$ndyn == 0 && $nstat == 0 && $ntarg == 0} {
	echo No items to tolerance
	verbose $verbosity
	return
	}

#Cache some info.
   set tolTable(dyn) $tolDyn
   set tolTable(stat) $tolStat
   set tolTable(targ) $tolTarg
   set tolTable(fractinc) $fractinc
   set tolTable(rmsinc) $rmsinc
   set tolTable(converge) $converge
   set tolTable(ndyn) $ndyn
   set tolTable(nstat) $nstat
   set tolTable(ntarg) $ntarg
   set tolTable(fom) [fom $optic]
   set tolTable(fomType) $_optic(fomType)

#fom0 is unperturbed f.o.m.
   set fom0 $tolTable(fom)


#Now loop through, run "dither" for static elements.
#Don't fork, because evaluation is relatively fast and I need to know
#the cumulative change in fom to set targets for dynamic parameters.
   set cum 0.
   foreach item $tolStat {
	set surfs [lindex $item 0]
	set params [lindex $item 1]
	set weight $tolTable($surfs,$params,weight)

#Tolerance value is static
	set tol $tolTable($surfs,$params,tol)

	if {[ctype digit [lindex $params 0]]} {
	   set showIncCmd showIndexInc
	} else {
	   set showIncCmd showSurfInc
	   }

#Increment is always takes from the first surface if multiple are being
#linked.
	set inc [$showIncCmd $optic [lindex $surfs 0] [lindex $params 0]]
	set fract [expr $tol/$inc]
	set list [dither $optic $surfs $params $fract]

#We use average change in fom.
	set diff1 [lindex $list 1]
	set diff2 [lindex $list 2]

#avg is the increment in fom: avg = fom1 - fom0.
	set avg [expr ($diff1+$diff2)/2.]

#I now want the difference done in quadrature:
	set rms [expr sqrt(pow($fom0+$avg,2) - pow($fom0,2))]
	set tolTable($surfs,$params,rms) $rms

#Just for kicks, save the difference as well.  I ought to add code to save
#difference down for other types of tolerances as well.  In some cases, I get
#a strong asymmetry.  Asymmetry reflects the fact that we may not be at an
#exact minimum w.r.t. this parameter in the default design.
	set tolTable($surfs,$params,asym) [expr $diff2-$diff1]

#Cum is the sum of the increments in chi^2.
	set cum [expr $cum + $weight * ($rms*$rms)]
	}

#Run ditherStudy on "target" parameters - these have a fixed incremental
#f.o.m.  Implement so it is parallelizable.
   forkEach item $tolTarg {
	set surfs [lindex $item 0]
	set params [lindex $item 1]
	set weight $tolTable($surfs,$params,weight)
	if {[verbose]} {echo Tolerancing $item}

#Target increment in f.o.m. for this parameter.
	set rms $tolTable($surfs,$params,rms)
	set target [expr sqrt($fom0*$fom0 + $rms*$rms)]
	set delta [expr $target - $fom0]
	if {[catch {set tolerance [ditherStudy $optic $surfs $params \
		0 $target $converge]} msg]} {
	   echo Disaster for $item $msg
	   set file error-[join $surfs #]-[join $params #]
	   set fid [open $file w]
	   puts $fid "Surface(s) $surfs parameter(s) $params, \
		error in ditherStudy"
	   puts $fid "$msg"
	   close $fid
	} else {

#Need to separate surfs and param names in case param is a number (refr. index)
	   set file tol-[join $surfs #]-[join $params #]
	   set fid [open $file w]
	   puts $fid [list $surfs $params $tolerance]
	   close $fid
	   }
	}
   foreach item $tolTarg {
	set surfs [lindex $item 0]
	set params [lindex $item 1]

#Target increment in f.o.m. for this parameter.  Add to "cum"
	set rms $tolTable($surfs,$params,rms)
	set weight $tolTable($surfs,$params,weight)
	set cum [expr $cum + $weight * ($rms*$rms)]
	set file tol-[join $surfs #]-[join $params #]
	set fid [open $file]
	set line [gets $fid]
	close $fid
	set tolerance [lindex $line 2]
	set tolTable($surfs,$params,tol) $tolerance
	if {[verbose]} {
	   echo Surface(s) $surfs $params tolerance \
		$tolTable($surfs,$params,tol)
	   }
	exec rm $file
	}

#Run "ditherStudy" on dynamic parameters and cache results.
#Adjust fractinc, rmsinc to be per parameter.
#If no dynamic params, just set minimal defaults and return.
   if {$ndyn == 0} {
	set tolTable(fomtarg) [expr sqrt($fom0*$fom0 + $cum)]
	verbose $verbosity
	return
	}

#Target worst-case f.o.m.
   set fomtarg [expr $fom0*(1.+$fractinc)]
   set fomtarg [expr max($fomtarg, sqrt($fom0*$fom0 + $rmsinc*$rmsinc))]
   set tolTable(fomtarg) $fomtarg

#dyntarg is target figure-of-merit excluding static and "target" contributions.
#Technically, it is the target chi^2.  cum is the sum of the incremental
#chi^2 of the static and "target" contributions.
   set dyntarg2 [expr $fomtarg*$fomtarg - $cum - $fom0*$fom0]
   if {$dyntarg2 < 0.} {
	echo Static contributions exceed dynamic target
	verbose $verbosity
	return
	}

#Target incremental f.o.m. for each parameter.  Use this as a "fixed"
#target in ditherStudy
#I will try an approximation that should work for small and large fractional
#increments.  For small increments, delta varies as 1./ndyn.  For large
#increments, delta varies as 1./sqrt(ndyn)
   set rmstarg [expr sqrt($dyntarg2/$ndyn)]

#Implement so it is parallelizable.
   forkEach item $tolDyn {
	set surfs [lindex $item 0]
	set params [lindex $item 1]
	set weight $tolTable($surfs,$params,weight)
	if {[verbose]} {echo Tolerancing $item}

#delta is target increment in f.o.m. for this parameter.
	set target [expr sqrt($rmstarg*$rmstarg + $fom0*$fom0)]
	if {[catch {set tolerance [ditherStudy $optic $surfs $params \
	   0 $target $converge]} msg]} {
	   echo Disaster for $item $msg
	   set file error-[join $surfs #]-[join $params #]
	   set fid [open $file w]
	   puts $fid "Surface(s) $surfs parameter(s) $params, error in \
		ditherStudy"
	   puts $fid "$msg"
	   close $fid
	} else {
	   set file tol-[join $surfs #]-[join $params #]
	   set fid [open $file w]
	   puts $fid [list $surfs $params $tolerance]
	   close $fid
	   }
	}
   foreach item $tolDyn {
	set surfs [lindex $item 0]
	set params [lindex $item 1]
	set file tol-[join $surfs #]-[join $params #]
	set fid [open $file]
	set line [gets $fid]
	close $fid
	set tolerance [lindex $line 2]
	set tolTable($surfs,$params,tol) $tolerance

#Need to set this outside the forkEach loop.  delta is same for all params.
	set tolTable($surfs,$params,rms) $rmstarg
	if {[verbose]} {
	   echo Surface(s) $surfs $params tolerance \
		$tolTable($surfs,$params,tol)
	   }
	exec rm $file
	}
   verbose $verbosity
   return
   }

#########################################################################
#List tolerances
proc tolList {{file ""}} {
   global tolTable _optic
   if {![info exists tolTable]} return
   if {$file == ""} {
	set stdout stdout
   } else {
	set stdout [open $file w]
	}
   set tolList $tolTable()
   puts $stdout "Baseline F.O.M.  : [format %.4f $tolTable(fom)]"
   puts $stdout "Target F.O.M.    : [format %.4f $tolTable(fomtarg)]"
   puts $stdout "Fractional incr. : [format %.4f $tolTable(fractinc)]"
#   puts $stdout "Fixed target     : [format %.4f $tolTable(fixed)]"
   puts $stdout "Convergence      : [format %.4f $tolTable(converge)]"
   puts $stdout "N(dyn)           : [format %d $tolTable(ndyn)]"
   puts $stdout "N(static)        : [format %d $tolTable(nstat)]"
   puts $stdout "N(target)        : [format %d $tolTable(ntarg)]"
   puts $stdout "FOM Types:       : $tolTable(fomType)"
   puts $stdout "Adjust params    : $tolTable(flags)"
   puts $stdout "Number of spots  : [llength $tolTable(spots)]"

   set format "%-8s %-8s %-8s %-8s %-8s %-8s"
   puts $stdout [format $format Surface Param Type "F.O.M." \
Tolerance "RMS F.O.M."]
   puts $stdout [format $format "" "" "" "Increment" "" ""]
   foreach item $tolList {
	set surfs [lindex $item 0]
	set params [lindex $item 1]
	set rms $tolTable($surfs,$params,rms)
	set fom0 $tolTable(fom)
	set fom1 [expr sqrt($rms*$rms + $fom0*$fom0)]
	set fomlim [expr $fom1-$fom0]
	foreach surf $surfs {
	   foreach param $params {
		puts $stdout "[format $format $surf \
		   $param $tolTable($surfs,$params,type) \
		[format %.2g $fomlim] \
		[format %.2g $tolTable($surfs,$params,tol)] \
		[format %.2g $rms]]"
		}
	   }
	}
   if {$file != ""} {close $stdout}
   return
   }

#########################################################################
#Scale tolerances.  RMS scales like scale

proc tolScale {scale} {
   global tolTable _optic
   set tolList $tolTable()
   foreach item $tolList {
	set surfs [lindex $item 0]
	set params [lindex $item 1]
	set tolTable($surfs,$params,tol) [expr $tolTable($surfs,$params,tol) \
	   * $scale]
	set tolTable($surfs,$params,rms) [expr \
	   $tolTable($surfs,$params,rms) * $scale]
	}
   return
   }

#########################################################################
#Copy least squares parameters from tolTable to _optic.  This is useful
#if I read tolerance table back in and want to know what the compensation
#procedure is.

proc tolToOptic {} {
   global tolTable _optic
   if {![info exists tolTable(flags)] || ![info exists tolTable(lscale)] \
	|| ![info exists tolTable(spots)] || ![info exists tolTable(niter)] \
	|| ![info exists tolTable(fomType)]} {
	echo "Missing one of flags, lscale, spots, niter, fomType \
	   from tolTable"
	return
	}
   set _optic(flags) $tolTable(flags)
   set _optic(flagcache) $tolTable(flags)
   set _optic(lscale) $tolTable(lscale)
   set _optic(spots) $tolTable(spots)
   set _optic(fomType) $tolTable(fomType)
   set _optic(niter) $tolTable(niter)
   return
   }

#########################################################################
#Run a Monte Carlo simulation to determine the distribution function of the
#Figure of Merit given a table of max. tolerances.
#The algorithm is to perturb all inputs randomly, adjust all compensators,
#and recompute the F.O.M.  Repeat a zillion times.

#If niter = 0, just dither parameters once and return optics structure.

proc monteCarlo {optic {niter 0}} {
   global tolTable
   set command forkLoop
   set save 0
   if {$niter == 0} {
	set command loop
	set niter 1
	set save 1
	}
   $command i 0 $niter {
	set opticinc [opticNew]
	opticCopy $optic $opticinc

#Need to randomize the random number generator!
	loop j 0 [expr $i*[llength $tolTable()]] {
	   expr rand()
	   }
	foreach item $tolTable() {
	   set surfs [lindex $item 0]
	   set params [lindex $item 1]
	   set tolerance $tolTable($surfs,$params,tol)

#Increment is +/- max(tolerance)
	   set inc [expr 2.*(rand()-.5)*$tolerance]

#If param is a digit, this is a refraction index.
	   foreach surf $surfs {
		foreach param $params {
		   if {[ctype digit $param]} {
			setIndex $opticinc $surf $param [expr \
		   	   [showIndex $opticinc $surf $param] + $inc]
		   } else {
			setSurf $opticinc $surf $param [expr \
		   	   [showSurf $opticinc $surf $param] + $inc]
			}
		   }
		}
	   }
	if {$save == 1} {
	   fastLstsq $opticinc
	   return $opticinc
	   }
	set fom [fom $opticinc]
	set file mc$i
	set fid [open $file w]
	puts $fid "$fom"
	close $fid
	opticDel $opticinc
	}
   set fomlist ""
   loop i 0 $niter {
	set file mc$i
	set fid [open $file]
	set fom [gets $fid]
	close $fid
	exec rm $file
	lappend fomList $fom
	}
   set tolTable(monteCarlo) $fomList
   return
   }

###################################################################
#Plot results of Monte Carlo
#Also compute mean, std dev

proc monteCarloPlot {} {
   global tolTable
   if {![info exists tolTable(monteCarlo)]} return
   plotInit a
   set list [lsort -real $tolTable(monteCarlo)]
   set ymax [lindex $list end]
   set ymin [lindex $list 0]
   set fom $tolTable(fom)
   set sum 0.
   set sum2 0.

#Minimum
   set ymin [expr min($ymin, $fom)]

#Target
   set target $tolTable(fomtarg)
   set ymax [expr max($ymax, $target)]
   pgSci 1
   set n [llength $list]
   pgEnv 0 $n $ymin $ymax 0 0
   loop i 0 $n {
	set val [lindex $list $i]
	pgPoint $i $val 3
	set sum [expr $sum + $val]

#Stdev is relative to f.o.merit
	set sum2 [expr $sum2 + ($val-$fom)*($val-$fom)]
	}
   pgSci 7
   pgLine "0 $n" "$fom $fom"
   pgSci 6
   pgLine "0 $n" "$target $target"
   set mean [expr $sum/$n]
   set stdev [expr sqrt($sum2/$n)]
   pgSci 1
   pgLabel "N" "Figure of Merit" "Monte Carlo: $tolTable(name)"
   pgSci 5
   echo Mean [format %.4f $mean] Stdev [format %.4f $stdev] \
	Target [format %.4f $target]
   pgText [expr .05*$n] [expr $ymin + .93*($ymax - $ymin)] \
	   "Base [format %.4f $fom]"
   pgText [expr .05*$n] [expr $ymin + .86*($ymax - $ymin)] \
	   "Mean [format %.4f $mean]"
   pgText [expr .05*$n] [expr $ymin + .79*($ymax - $ymin)] \
	   "Stdev [format %.4f $stdev]"
   pgText [expr .05*$n] [expr $ymin + .72*($ymax - $ymin)] \
	   "Target [format %.4f $target]"

   return
   }

#######################################################################
#Compute covariances between parameters.  I will assume that I have
#a tolerance table where we have NO compensators.

#To do covariances, make sure we use a healthy fractinc (e.g., 1) so
#fractincs for each variable are sizable.
#Also, make sure we run lstsqs to completion.
proc covar {optic} {
   global tolTable

   set ncombo [llength $tolTable()]
   set opticinc [opticNew]
   loop i 0 $ncombo {
	set combo1 [lindex $tolTable() $i]
	set surfs1 [lindex $combo1 0]
	set params1 [lindex $combo1 1]
	loop j [expr $i+1] $ncombo {
	   set combo2 [lindex $tolTable() $j]
	   set surfs2 [lindex $combo2 0]
	   set params2 [lindex $combo2 1]
	   set inc1 [expr .707*$tolTable($surfs1,$params1,tol)]
	   set inc2 [expr .707*$tolTable($surfs2,$params2,tol)]

#Increment is +/- .707*max(tolerance)

#+x +y
	   opticCopy $optic $opticinc
	   foreach surf1 $surfs1 {
		foreach param1 $params1 {
		   setSurf $opticinc $surf1 $param1 [expr \
			[showSurf $opticinc $surf1 $param1] + $inc1]
		   }
		}
	   foreach surf2 $surfs2 {
		foreach param2 $params2 {
		   setSurf $opticinc $surf2 $param2 [expr \
			[showSurf $opticinc $surf2 $param2] + $inc2]
		   }
		}
	   set fom11 [fom $opticinc]
#+x -y
	   opticCopy $optic $opticinc
	   foreach surf1 $surfs1 {
		foreach param1 $params1 {
		   setSurf $opticinc $surf1 $param1 [expr \
			[showSurf $opticinc $surf1 $param1] + $inc1]
		   }
		}
	   foreach surf2 $surfs2 {
		foreach param2 $params2 {
		   setSurf $opticinc $surf2 $param2 [expr \
			[showSurf $opticinc $surf2 $param2] - $inc2]
		   }
		}
	   set fom12 [fom $opticinc]
#-x y
	   opticCopy $optic $opticinc
	   foreach surf1 $surfs1 {
		foreach param1 $params1 {
		   setSurf $opticinc $surf1 $param1 [expr \
			[showSurf $opticinc $surf1 $param1] - $inc1]
		   }
		}
	   foreach surf2 $surfs2 {
		foreach param2 $param2 {
		   setSurf $opticinc $surf2 $param2 [expr \
			[showSurf $opticinc $surf2 $param2] + $inc2]
		   }
		}
	   set fom21 [fom $opticinc]
#-x -y
	   opticCopy $optic $opticinc
	   foreach surf1 $surfs1 {
		foreach param1 $params1 {
		   setSurf $opticinc $surf1 $param1 [expr \
			[showSurf $opticinc $surf1 $param1] - $inc1]
		   }
		}
	   foreach surf2 $surfs2 {
		foreach param2 $params2 {
		   setSurf $opticinc $surf2 $param2 [expr \
			[showSurf $opticinc $surf2 $param2] - $inc2]
		   }
		}
	   set fom22 [fom $opticinc]
	   set ds [expr (($fom11 + $fom22) - ($fom21 + $fom12))/4.]
	   set tolTable($params1,$surfs1,$params2,$surfs2,covar) \
		[expr $ds/(($tolTable(fomtarg)-$tolTable(fom)))]
	   echo $combo1 $combo2 [format %.4g \
		$tolTable($params1,$surfs1,$params2,$surfs2,covar)]
	   }
	}
   opticDel $opticinc
   return
   }

#######################################################################
#Save a tolerance table in a form that can be reread.

proc tolSave {{n 1}} {
   set var tolTable$n
   global $var
   if {[info exists $var]} {unset $var}
   global tolTable
   array set $var [array get tolTable]
   return
   }

###########################################################
#Merge tolerance tables
proc tolMerge {{n 1}} {
   global tolTable
   set var tolTable$n
   global $var
   foreach combo [set ${var}()] {
	lappend tolTable() $combo
	set surfs [lindex $combo 0]
	set params [lindex $combo 1]
	set tolTable($surfs,$params,tol) [set ${var}($surfs,$params,tol)]
	}
   set tolTable(fractinc) [expr $tolTable(fractinc) + \
	[set ${var}(fractinc)]]
   set tolTable(ntol) [expr $tolTable(ntol) + [set ${var}(ntol)]]

#Don't reset tolTable(fixed) yet - too much logic!
   return
   }
