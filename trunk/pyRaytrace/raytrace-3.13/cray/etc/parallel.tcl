#TCL proc for running least squares in parallel mode
#
#Procedures:
#	plstsq <optic> [lscale]
#	pfastLstsq <optic>

proc plstsq {optic {lscale ""}} {
   global _optic
   global FORK

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
#lstsq is the master template.

   set lstsq [genericNew LSTSQ]

#Use a separate proc to set spot positions
   spotSet $lstsq
   set nspotdef [exprGet $lstsq.nspotdef]

#Backward compatibility:
#   If lscale is supplied, set all spot flags to this value.
   if {$lscale != ""} {
        loop i 0 [exprGet $lstsq.nspotdef] {
           handleSet $lstsq.lscale<$i> $lscale
           }
        }

#Set up master array.  I will copy to working arrays later.
#CAUTION!!! The LSTSQ structure contains the parameterized spot positions
#that are updated by least squares, plus the iterator count.  These
#need to be preserved in successive iterations!  I do not want to copy
#pointers for temporary storage, however, which are allocated in lstsq1a.
   lstsq1 $optic $lstsq
   if {[verbose]} {flagList $optic}

#Run one iteration to adjust the spot centers.  This is not parallelized.
   lstsq1a $lstsq
   lstsq1b $lstsq
   lstsq2 $optic $lstsq

#The following is a bit of a hack
   handleSet $lstsq.fractmax 1.
   lstsq3 $optic $lstsq
   lstsq4 $lstsq

#Parallelize.

   set niter 0

   set totiter 1
   while {1} {

	if {[verbose]} {
	   puts stdout "Iteration $niter"
	   rmsList $optic $lstsq
	   }
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

#Clear out any temporary storage in master template from prior iterations,
#copy to working structures, then allocate temporary storage in all.
#Clean up working structures when done.
	lstsq4 $lstsq

#Clear out old rms sums
	lstsq1b $lstsq

#Allocate working lstsq structures
	loop i 0 $nspotdef {
	   set array($i) [genericNew LSTSQ]

#Copy master template
	   handleSetFromHandle $array($i) $lstsq
	   }

#Allocate temporary arrays in master template.
	lstsq1a $lstsq

#Loop through all spots and working structures.
	forkLoop i 0 $nspotdef {
	   lstsq1a $array($i)

#This step is probably not needed.
	   lstsq1b $array($i)

#Zero out weights on all spots except i.
	   loop j 0 $nspotdef {
		if {$j != $i} {handleSet $array($i).weight<$j> 0}
		}
#There might be spots that don't make it to the focal plane.
	   if {[file exists lstsq$i.dat]} {exec rm lstsq$i.dat}
	   if {![catch {lstsq2 $optic $array($i)}]} {
		if {[verbose]} {echo Computing spot $i}
		set fid [open lstsq$i.dat w]
		binWrite $fid $array($i)
		set ndim [exprGet $array($i).ndim]
		set ncmat [expr $ndim*($ndim+1)/2]
		binWrite $fid $array($i).rhs<0> [expr $ndim*8]
		binWrite $fid $array($i).cmat<0> [expr $ncmat*8]
		close $fid
		}
	   }

	loop i 0 $nspotdef {
	   if {[verbose]} {echo Adding in spot $i}

#Here is where I read in results of forkLoop.
	   lstsq4 $array($i)
	   if {![file exists lstsq$i.dat]} {
		genericDel $array($i)
		continue
		}
	   set fid [open lstsq$i.dat]

#binRead seems to have a bizarre behavior that it skips to the end of file
#regardless of the amount read.  Use seek to reset
	   set nread [binRead $fid $array($i)]
	   seek $fid $nread start
	   set ndim [exprGet $array($i).ndim]
	   lstsq1a $array($i)
	   set ncmat [expr $ndim*($ndim+1)/2]
	   set n [binRead $fid $array($i).rhs<0> [expr $ndim*8]]
	   set nread [expr $nread+$n]
	   seek $fid $nread start
	   binRead $fid $array($i).cmat<0> [expr $ncmat*8]
	   close $fid
	   exec rm lstsq$i.dat
	   lstsqSum $lstsq $array($i)

#Delete temparary arrays
	   lstsq4 $array($i)
	   genericDel $array($i)
	   }

	incr niter

#The following is a bit of a hack
        handleSet $lstsq.fractmax $fractmax
	lstsq3 $optic $lstsq
	if {[verbose]} {updateList $optic}
	set fract [exprGet $lstsq.fract]
	if {$fract == 0.} break
	}

   lstsq4 $lstsq

#More caching
   foreach var "rmsx rmsy nxy rmscx rmscy ncxy fractmax" {
	set $var [exprGet $lstsq.$var]
	}
   set _optic(niter) $niter
   set _optic(rmsx) $rmsx
   set _optic(rmsy) $rmsy
   set _optic(nxy) $nxy
   set _optic(rmscx) $rmscx
   set _optic(rmscy) $rmscy
   set _optic(ncxy) $ncxy
   set _optic(fractmax) $fractmax
   genericDel $lstsq
   return
   }

########################################################################
#Repeat a least squares fit based on cached parameters from a run using
#plstsq above.

proc pfastLstsq {optic} {
   global _optic
   global FORK

#Cache current optic structure
   if {![info exists _optic(opticache)] || \
	   [catch {handleGet $_optic(opticache)}]} {
	set _optic(opticache) [opticNew]
	}
   opticCopy $optic $_optic(opticache)

   resetFlags $optic
   flagSet $optic
   cacheFlags
   set lscale $_optic(lscale)
   set fractmax $_optic(fractmax)

#Initialize
#lstsq is the master template.

   set lstsq [genericNew LSTSQ]

   spotSet $lstsq
   set nspotdef [exprGet $lstsq.nspotdef]

#Backward compatibility:
#   If lscale is supplied, set all spot flags to this value.
   if {$lscale != ""} {
        loop i 0 [exprGet $lstsq.nspotdef] {
           handleSet $lstsq.lscale<$i> $lscale
           }
        }

#Set up master array.  I will copy to working arrays later.
#CAUTION!!! The LSTSQ structure contains the parameterized spot positions
#that are updated by least squares, plus the iterator count.  These
#need to be preserved in successive iterations!  I do not want to copy
#pointers for temporary storage, however, which are allocated in lstsq1a.
   lstsq1 $optic $lstsq

#Run one iteration to adjust the spot centers.  This is not parallelized.
   lstsq1a $lstsq
   lstsq1b $lstsq
   lstsq2 $optic $lstsq

#The following is a bit of a hack
   handleSet $lstsq.fractmax 1.
   lstsq3 $optic $lstsq
   lstsq4 $lstsq

   set niter $_optic(niter)
   loop iter 0 $niter {

#Clear out any temporary storage in master template from prior iterations,
#copy to working structures, then allocate temporary storage in all.
#Clean up working structures when done.
	lstsq4 $lstsq

#Clear out old rms sums
	lstsq1b $lstsq

#Allocate working lstsq structures
	loop i 0 $nspotdef {
	   set array($i) [genericNew LSTSQ]

#Copy master template
	   handleSetFromHandle $array($i) $lstsq
	   }

#Allocate temporary arrays in master template.
	lstsq1a $lstsq

#Loop through all spots and working structures.
	forkLoop i 0 $nspotdef {
	   lstsq1a $array($i)

#This step is probably not needed.
	   lstsq1b $array($i)

#Zero out weights on all spots except i.
	   loop j 0 $nspotdef {
		if {$j != $i} {handleSet $array($i).weight<$j> 0}
		}
#There might be spots that don't make it to the focal plane.
	   if {[file exists lstsq$i.dat]} {exec rm lstsq$i.dat}
	   if {![catch {lstsq2 $optic $array($i)}]} {
		set fid [open lstsq$i.dat w]
		binWrite $fid $array($i)
		set ndim [exprGet $array($i).ndim]
		set ncmat [expr $ndim*($ndim+1)/2]
		binWrite $fid $array($i).rhs<0> [expr $ndim*8]
		binWrite $fid $array($i).cmat<0> [expr $ncmat*8]
		close $fid
		}
	   }

	loop i 0 $nspotdef {

#Here is where I read in results of forkLoop.
	   lstsq4 $array($i)
	   if {![file exists lstsq$i.dat]} {
		genericDel $array($i)
		continue
		}
	   set fid [open lstsq$i.dat]

#binRead seems to have a bizarre behavior that it skips to the end of file
#regardless of the amount read.  Use seek to reset
	   set nread [binRead $fid $array($i)]
	   seek $fid $nread start
	   set ndim [exprGet $array($i).ndim]
	   lstsq1a $array($i)
	   set ncmat [expr $ndim*($ndim+1)/2]
	   set n [binRead $fid $array($i).rhs<0> [expr $ndim*8]]
	   set nread [expr $nread+$n]
	   seek $fid $nread start
	   binRead $fid $array($i).cmat<0> [expr $ncmat*8]
	   close $fid
	   exec rm lstsq$i.dat
	   lstsqSum $lstsq $array($i)

#Delete temparary arrays
	   lstsq4 $array($i)
	   genericDel $array($i)
	   }
	lstsq3 $optic $lstsq
	set fract [exprGet $lstsq.fract]
	if {$fract == 0.} break
	}

   lstsq4 $lstsq

#More caching
   foreach var "rmsx rmsy nxy rmscx rmscy ncxy" {
	set $var [exprGet $lstsq.$var]
	}
   set _optic(niter) $niter
   set _optic(rmsx) $rmsx
   set _optic(rmsy) $rmsy
   set _optic(nxy) $nxy
   set _optic(rmscx) $rmscx
   set _optic(rmscy) $rmscy
   set _optic(ncxy) $ncxy
   genericDel $lstsq
   return
   }


	
	
