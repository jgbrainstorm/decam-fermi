#Procedures for reading and computing glass properties.

###################################################################
#Sellmeier equation - used by Schott

proc sellmeier {wave b1 b2 b3 c1 c2 c3} {
   set w2 [expr pow($wave,2)]
   if {[catch {set n2 [expr (1. + $b1*$w2/($w2 - $c1) + \
	       $b2*$w2/($w2 - $c2) + \
	       $b3*$w2/($w2 - $c3))]}]} {return 0.}
   if {$n2 < 0} {
	set n -1.
   } else {
	set n [expr sqrt($n2)]
	}
   return $n
   }

###################################################################
#Cauchy equation - used by Ohara, Hoya.
#Note documentation is extremely sketchy.  I pulled this eqn from
#a Pilkington PDF.  Don't know if it works for Hoya
#It seems the powers of lambda can vary in different implementations of 
#the equation - probably need to check with vendors to figure out
#exactly how it is done.

proc cauchy {wave a0 a1 a2 a3 a4 a5} {
   set w2 [expr pow($wave,2)]

   set n2 [expr $a0 + $a1*$w2 + $a2/$w2 + $a3*pow($w2,-2) + $a4*pow($w2,-3) \
	+ $a5*pow($w2,-4)]
   return [expr sqrt($n2)]
   }

###################################################################
#Read dispersion coefficients

proc dispRead {} {
   global disp env
   if {[info exists disp]} {unset disp}
   set disp() ""

#Data file is DISPERSION.TXT
   set file DISPERSION.TXT
   if {![info exists $file]} {
	set file $env(CRAY_DIR)/etc/$file
	}
   if {![file exists $file]} {
	error "Cannot find file $file!"
	}
   echo Reading $file
   set fid [open $file]
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set name [lindex $line 0]

#Now have a proper SIO2 table.
	if {$name == "Quartz"} {set name quartz}
	if {$name == "vacuum"} {set name vac}
	if {$name == "Fluorite"} {set name CaF2}
	set maker [lindex $line 1]
	set type [lindex $line 2]
	set coeffs [lrange $line 3 end]
	if {[lsearch $disp() $name] < 0} {
	   lappend disp() $name
	   set disp($name,maker) $maker
	   set disp($name,type) $type
	   loop i 0 6 {
		set disp($name,c$i) [lindex $coeffs $i]
		}
	   }
	}
   close $fid
   return
   }

###################################################################
#Compute dispersion function for the specified glass
#Return 0 if glass not found, but don't print error message.
proc disp {glass wave} {
   global disp
   if {![info exists disp]} {dispRead}
   if {[lsearch $disp() $glass] < 0} {
	return 0
	}
   set type $disp($glass,type)
   loop i 0 6 {
	set c$i $disp($glass,c$i)
	}
   if {$type == "c"} {
	set n [cauchy $wave $c0 $c1 $c2 $c3 $c4 $c5]
	return $n
	}
   if {$type == "s"} {
	set n [sellmeier $wave $c0 $c1 $c2 $c3 $c4 $c5]
	return $n
	}
   return 0
   }

###################################################################
#See what we have in OPUSGLASS.TXT
proc glassRead {} {
   global glass env
   if {[info exists glass]} {unset glass}
   set glass() ""

#Data file is OPUSGLASS.TXT
   set file OPUSGLASS.TXT
   if {![file exists $file]} {
	set file $env(CRAY_DIR)/etc/$file
	}
   if {![file exists $file]} {
	error "Cannot find file $file!"
	}
   echo Reading $file
   set fid [open $file]
   set name ""
   set maker ""
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {$line == ""} continue
	set func [lindex $line 0]
	if {$func == "name"} {
	   set name [lindex $line 2]

#Now have a proper SIO2 table.
	   if {$name == "Quartz"} {set name quartz}
	   if {$name == "vacuum"} {set name vac}
	   if {$name == "Fluorite"} {set name CaF2}
	   }
	if {$func == "maker"} {
	   set maker [lindex $line 2]
	   if {[lsearch $glass() $name] < 0} {
		lappend glass() $name
		set glass($name) ""
		set glass($name,maker) $maker
		set glass($name,wave) ""
		set glass($name,index) ""
		}
	   }
	if {$func == "n"} {
	   if {$name == "" || $maker == ""} continue
	   set wave [lindex $line 1]
	   set n [lindex $line 3]
	   lappend glass($name,wave) $wave
	   lappend glass($name,index) $n
	   }
	}
   close $fid
   return
   }

######################################################################
#Read transmission data.  Data are for 25 mm thickness glass.
proc transRead {} {
   global transmit env
   if {[info exists transmit]} {unset transmit}
   set transmit() ""

#Data file is TRANSMIT.TXT
   set file TRANSMIT.TXT
   if {![file exists $file]} {
	set file $env(CRAY_DIR)/etc/$file
	}
   if {![file exists $file]} {
	error "Cannot find file $file!"
	}
   echo Reading $file
   set fid [open $file]
   set name ""
   set maker ""
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {$line == ""} continue
	set func [lindex $line 0]
	if {$func == "name"} {
	   set name [lindex $line 2]
	   }
	if {$func == "maker"} {
	   set maker [lindex $line 2]
	   if {[lsearch $transmit() $name] < 0} {
		lappend transmit() $name
		set transmit($name) ""
		set transmit($name,maker) $maker
		set transmit($name,wave) ""
		set transmit($name,trans) ""
		}
	    }
	if {$func == "t"} {
	   if {$name == "" || $maker == ""} continue
	   set wave [lindex $line 1]
	   set n [lindex $line 3]
	   lappend transmit($name,wave) $wave
	   lappend transmit($name,trans) $n
	   }
	}
   close $fid
   return
    }

####################################################################
#Interpolate in glass table from opus.  This code is depracated because I
#really want to use dispersion formulae (plus, OPUS must have computed its
#table from those formulae).  However, I keep it here for a) legacy uses;
#b) as a cross-check.
#Input is in microns.  Opus uses angstroms, so
#I need to multiply by 1.e4 to convert.
#Return 0 if glass not found.

proc opus {name wave} {
   global glass
   if {![info exists glass]} {
	glassRead
	}
   if {[lsearch $glass() $name] < 0} {return 0}

   set wave [expr $wave*1.e4]
   set i 1
   loop j 1 [llength $glass($name,wave)] {
	if {$wave > [lindex $glass($name,wave) [expr $j-1]] && \
	   $wave <= [lindex $glass($name,wave) $j]} {
	   set i $j
	   break
	   }
	}
   if {$wave > [lindex $glass($name,wave) end]} {
	set i [expr [llength $glass($name,wave)] - 1]
	}
   set w0 [expr 1.*[lindex $glass($name,wave) [expr $i-1]]]
   set w1 [expr 1.*[lindex $glass($name,wave) $i]]
   set n0 [lindex $glass($name,index) [expr $i-1]]
   set n1 [lindex $glass($name,index) $i]
   set n [expr $n0 + ($wave-$w0) * ($n1-$n0) / ($w1-$w0)]
   return $n
   }

####################################################################
#Interpolate in glass table.  Input is in microns.  Opus uses angstroms, so
#I need to multiply by 1.e4 to convert.

proc glass {name wave} {
   global glass GLASSMODE
   if {![info exists GLASSMODE]} {
	set GLASSMODE glass
	}

#If GLASSMODE is not glass, try dispersion formulae first.
   if {$GLASSMODE == "glass"} {
	set n [opus $name $wave]
	if {$n == 0} {
	   set n [disp $name $wave]
	   if {$n == 0} {
		error "No glass $name"
		}
	   echo Using dispersion relation
	   }
   } else {
	set n [disp $name $wave]

#Don't use opus files if we are in "disp" mode.
#	if {$n == 0} {
#	   set n [opus $name $wave]
#	   if {$n == 0} {
#		error "No glass $name"
#		}
#	   echo Using interpolation
#	   }
	if {$n == 0} {
	   error "No glass $name"
	   }
	}
   return $n
   }

####################################################################
#Interpolate in transmit table.  Input is in microns.  Opus uses angstroms, so
#I need to multiply by 1.e4 to convert.  Transmission is for 25 mm reference
#thickness.

proc transmit {name wave} {
   global transmit
   if {![info exists transmit]} {
	transRead
	}
   set wave [expr $wave*1.e4]

#Not all glasses have transmission.  I will return 0 if not so there is an
#obvious indicator of a problem.
   if {[lsearch $transmit() $name] < 0} {
	return 0.
	}
   set i 1
   loop j 1 [llength $transmit($name,wave)] {
	if {$wave > [lindex $transmit($name,wave) [expr $j-1]] && \
	   $wave <= [lindex $transmit($name,wave) $j]} {
	   set i $j
	   break
	   }
	}
   if {$wave > [lindex $transmit($name,wave) end]} {
	set i [expr [llength $transmit($name,wave)] - 1]
	}
   set w0 [expr 1.*[lindex $transmit($name,wave) [expr $i-1]]]
   set w1 [expr 1.*[lindex $transmit($name,wave) $i]]
   set n0 [lindex $transmit($name,trans) [expr $i-1]]
   set n1 [lindex $transmit($name,trans) $i]
   set n [expr $n0 + ($wave-$w0) * ($n1-$n0) / ($w1-$w0)]
   return $n
   }

#########################################################################
#Set the name of a glass for a surface.
#setparam in read.c does not yet handle glass names, so deal with innards
#of data structures explicitly.


###################################################################
#Read densities

proc densityRead {} {
   global density env
   if {[info exists density]} {unset density}
   set density() ""

#Data file is DENSITY.TXT
   set file DENSITY.TXT
   if {![info exists $file]} {
	set file $env(CRAY_DIR)/etc/$file
	}
   if {![file exists $file]} {
	error "Cannot find file $file!"
	}
   echo Reading $file
   set fid [open $file]
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set name [lindex $line 0]

#Now have a proper SIO2 table.
	if {$name == "Quartz"} {set name quartz}
	if {$name == "vacuum"} {set name vac}
	if {$name == "Fluorite"} {set name CaF2}
	set maker [lindex $line 2]
	set den [lindex $line 1]
	if {[lsearch $density() $name] < 0} {
	   lappend density() $name
	   set density($name,maker) $maker
	   set density($name,den) $den
	   }
	}
   close $fid
   return
   }

###############################################################
proc density {glass} {
   global density
   if {![info exists density]} {densityRead}
   if {[lsearch $density() $glass] < 0} {
	echo Using default density for $glass
	return $density(sio2,den)
	}
   return $density($glass,den)
   }

#################################################################
#################################################################
#Do not set the indices themselves - the choice of which index to set and
#what sign can be complex, depending on the situation.
proc setGlass {optic id name} {
   set isurf [surfIndex $optic $id]
   handleSet $optic.optic<$isurf>->glass $name
   return
   }

################################################################
proc showGlass {optic id} {
   set isurf [surfIndex $optic $id]
   return [lindex [exprGet $optic.optic<$isurf>->glass] 0]
   return
   }

#################################################################
#Set the comment field - I will use this to give surfaces/elements names.

proc setComment {optic id name} {
   set isurf [surfIndex $optic $id]
   handleSet $optic.optic<$isurf>->comment $name
   return
   }

################################################################
proc showComment {optic id} {
   set isurf [surfIndex $optic $id]
   return [lindex [exprGet $optic.optic<$isurf>->comment] 0]
   return
   }

#################################################################
#Alias for setComment

proc setName {optic id name} {
   setComment $optic $id $name
   return
   }

################################################################
proc showName {optic id} {
   return [showComment $optic $id]
   }

###############################################################
#Change the type of a glass in an existing optic structure.
#I loop through all colors for a specified surface and change the
#index for all non-zero indices.

#I need to know the sign of the index of the previous surface.
#I will be manic and fetch it for each color, even though it is a bit
#slow for lots of colors.

proc surfGlassSwitch {optic surf glass} {
   set ncolor [exprGet $optic.ncolor]
   for {set i 1} {$i <= $ncolor} {incr i} {
	set index [showIndex $optic $surf $i]
	if {$index == 0} continue
	set wave [showFocal $optic $i wave]
	set index [glass $glass $wave]
	setIndex $optic $surf $i $index
	}
   setGlass $optic $surf $glass
   return
   }

##################################################################
#Create a global array with standard wavelengths
proc waveStd {} {
   global waveStd waveName
   if {[info exists waveStd]} {unset waveStd}
   set waveStd(t) 1.01398
   set waveStd(A') .768195
   set waveStd(r)  .706519
   set waveStd(C)  .656273
   set waveStd(C') .643847
   set waveStd(D) .589294
   set waveStd(d) .587562
   set waveStd(e) .546074
   set waveStd(F) .486133
   set waveStd(F') .479991
   set waveStd(g) .435834
   set waveStd(h) .404656
   set waveStd(i) .365015

   foreach name [array names waveStd] {
	set waveName($waveStd($name)) $name
	}
   return
   }

#################################################################
#List a glass at standard wavelengths
proc glassList {glass} {
   global waveStd waveName
   waveStd
   puts stdout [format "%-4s %7s %8s %8s" Line Wave Index Eqn]
   foreach wave [lsort -real [array names waveName]] {
	set name $waveName($wave)
	set n [glass $glass $wave]
	set nd [disp $glass $wave]
	puts stdout [format "%-4s %7.4f %8.5f %8.5f" $name $wave $n $nd]
	}
   return
   }

########################################################################
#Wavelengths are normally ordered shortest to longest.
#Which wavelengths, exactly, do I want?  Schott uses d line.
#Others use the e line.

proc abbeList {glass {wave1 ""} {wave2 ""} {wave3 ""}} {
   global waveStd
   waveStd
   if {$wave1 == ""} {set wave1 $waveStd(F')}
   if {$wave2 == ""} {set wave2 $waveStd(e)}
   if {$wave3 == ""} {set wave3 $waveStd(C')}
   set ne [glass $glass $wave2]
   set nF [glass $glass $wave1]
   set nC [glass $glass $wave3]
   if {$nF == $nC} {return 0.}
   set abbe [expr ($ne-1.)/($nF - $nC)]
   return [format %.2f $abbe]
   }

#######################################################################
#Make a plot of index, dispersion for a given maker
#Additional arguments are specific glasses to highlight

proc abbePlot {maker args} {
   global glass waveStd
   eval abbeHandlePlot [list "" "" "" ""] $maker $args
   return
   }

#######################################################################
#Make a plot of index, dispersion for a given maker.
#Additional arguments are specific glasses to highlight
#Wavelengths are taken from the handle to an optic structure - this way
#I can use non-standard wavelenghths.

proc abbeHandlePlot {hndl ifil1 ifil2 ifil3 maker args} {
   global GLASSMODE
   global $GLASSMODE waveStd
   if {![info exists waveStd]} waveStd
   if {![info exists $GLASSMODE]} {eval ${GLASSMODE}Read}

   if {$hndl != ""} {
	set wave1 [showWave $hndl $ifil1]
	set wave2 [showWave $hndl $ifil2]
	set wave3 [showWave $hndl $ifil3]
	set waved $wave2
   } else {
	set wave1 ""
	set wave2 ""
	set wave3 ""
	set waved $waveStd(d)
	}
   plotInit a
   pgSci 1

#Compute abbe number of BK7 for these wavelengths
   set bk7 [abbeList SIO2 $wave1 $wave2 $wave3]
   set xmax [expr $bk7*1.5]
   set xmin [expr $bk7/3.]

#Range covers only common glasses for big optics
   pgEnv $xmax $xmin 1.4 1.7 0 0
   pgLabel "Abbe # $wave1 $wave2 $wave3" "Index(d)" "Dispersons for $maker"
   update

#Make a cache of all glasses that are plotted, so we can recall using the
#pointer later.
   global _glasscache
   if {[info exists _glasscache]} {unset _glasscache}
   set _glasscache() ""

#First, plot the specific requested glasses, regardless of type.
   foreach type [set ${GLASSMODE}()] {
	set i 1
	foreach arg $args {
	   incr i
	   if {$i > 7} {set i 2}
	   if {$arg == "$type"} {
		set abbe [abbeList $type $wave1 $wave2 $wave3]
		set n [glass $type $waved]
		pgSci $i
		pgPoint $abbe $n 2
		pgPoint $abbe $n 3
		pgPoint $abbe $n 4
		pgSch -1
		pgText [expr $abbe-1.] $n $arg
		pgSch 1
		pgSci 1
		if {[lsearch $_glasscache() $type] < 0} {
		   lappend _glasscache() $type
		   set _glasscache($type,n) $n
		   set _glasscache($type,abbe) $abbe
		   }
		}
	   }
	update
	}

#Next, plot all glasses for this maker.
   pgSch .2
   foreach type [set ${GLASSMODE}()] {
	if {[string compare [set ${GLASSMODE}($type,maker)] $maker] != 0} \
	   continue
	set abbe [abbeList $type $wave1 $wave2 $wave3]
	set n [glass $type $waved]
	pgPoint $abbe $n 2
	pgPoint $abbe $n 3
	pgPoint $abbe $n 4
	update
	if {[lsearch $_glasscache() $type] < 0} {
	   lappend _glasscache() $type
	   set _glasscache($type,n) $n
	   set _glasscache($type,abbe) $abbe
	   }
	}
   pgSch 1.
   return
   }

#####################################################################
#Get nearest glass to a point in the abbePlot diagram

proc nearestGlass {} {

#Click on a point
   echo Click on a point on the Abbe plot ...
   set list [pgCoord]
   set abbe [lindex $list 0]
   set n [lindex $list 1]

#Get plot limits
   set limits [pgQwin]
   set xl [lindex $limits 0]
   set xr [lindex $limits 1]
   set yb [lindex $limits 2]
   set yt [lindex $limits 3]

#If this is not an abbePlot, I can't really tell.  I would need to attach
#a name or other identifier.
   if {$abbe > $xl || $abbe < $xr} return
   if {$n < $yb || $n > $yt} return
   set dx [expr abs($xr - $xl)]
   set dy [expr abs($yt - $yb)]
   if {$dx == 0 || $dy == 0} return
   global _glasscache
   if {![info exists _glasscache]} return
   set distmin 999.
   set typemin ""
   foreach type $_glasscache() {
	set a1 $_glasscache($type,abbe)
	set n1 $_glasscache($type,n)
	set dist [expr pow(($a1-$abbe)/$dx,2) + pow(($n1-$n)/$dy,2)]
	if {$dist < $distmin} {
	   set typemin $type
	   set abbemin $a1
	   set nmin $n1
	   set distmin $dist
	   }
	}
   if {$typemin == ""} return
   pgSci 6
   pgPoint $abbemin $nmin 2
   pgPoint $abbemin $nmin 3
   pgPoint $abbemin $nmin 4
   pgSch -1
   pgText [expr $abbemin-1.] $n $typemin
   pgSch 1
   pgSci 1
   return $typemin
   }

#####################################################################
#Plot effective abbe number from specified surface of a design.

proc surfAbbe {hndl surfid ifil1 ifil2 ifil3} {

   set nF [expr abs([showIndex $hndl $surfid $ifil1])]
   set ne [expr abs([showIndex $hndl $surfid $ifil2])]
   set nC [expr abs([showIndex $hndl $surfid $ifil3])]
   set abbe [expr ($ne-1.)/($nF - $nC)]
   echo abbe = $abbe
   pgSci 6
   pgPoint $abbe $ne 2
   pgPoint $abbe $ne 3
   pgPoint $abbe $ne 4
   pgSch -1
   pgText [expr $abbe-1.] $ne "Surf $surfid"
   pgSch 1
   pgSci 1
   return
   }
