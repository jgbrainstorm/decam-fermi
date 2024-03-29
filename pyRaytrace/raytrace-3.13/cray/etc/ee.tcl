#Code to compute encircled energy vs. radius and make a plot.
#Personally, I think these plots are only marginally useful, but ray trace
#programs love to produce them, so I guess I will too.
#
#This procedure works from a file, but we specify filter as input.
#The file is generated by psfMap.
#
#If icolor is not numeric, treat it as the file name itself.  (E.G., I can
#use spotMap to create geometric profile
#
#How to use - we first create a .fit file of the PSF from either psfMap
#(for diffraction) or spotMap (for geometric).  For diffraction, we just
#input the color number.  For geometric, we input the file name "spot.fit".

proc eePlot {icolor {overlay n}} {
   global XMIN XMAX
   if {[ctype digit $icolor]} {
	set file fft$icolor.fit
   } else {
	set file $icolor
	}

   if {![file exists $file]} {
	error "No file $file!"
	}
   set reg [regReadFromFits $file]

#Get scale factors
   set flux [hdrGetAsDbl $reg.hdr FLUX]
   set pixmm [hdrGetAsDbl $reg.hdr PIXMM]
   set xmm [format %.2f [hdrGetAsDbl $reg.hdr XMM]]
   set ymm [format %.2f [hdrGetAsDbl $reg.hdr YMM]]
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]

#Centroid is in the middle of a pixel
   set row [expr $nrow/2+0.5]
   set col [expr $ncol/2+0.5]
   set radmax [min [expr $nrow/2] [expr $ncol/2]]
#   set chain [annulus $reg $row $col $radmax]

#Try annflux - better flux preservation at expense of more ragged profile.
   set chain [annflux $reg $row $col $radmax]
   set mmin 0
   set mmax [expr $radmax*abs($pixmm)]
   if {[info exists XMIN]} {set mmin $XMIN}
   if {[info exists XMAX]} {set mmax $XMAX}
   if {$overlay == "n"} {
	plotInit a
	pgSci 1
	pgEnv $mmin $mmax 0 1 0 0
	pgLabel "Radius (mm)" "Encircled Energy" "Filter $icolor xmm \
$xmm ymm $ymm"
	}
   set xpts ""
   set ypts ""
   global img
   chainEach hndl $chain {
	set rad [exprGet $hndl.rad]
	set cum [exprGet $hndl.cum]
	set mmrad [expr $rad*abs($pixmm)]
	set fract [expr $cum/$flux]
	lappend xpts $mmrad
	lappend ypts $fract
	}

#Keep current color index if we are overlaying.
   if {$overlay == "n"} {pgSci 5}
   pgLine $xpts $ypts
   genericChainDel $chain
   regDel $reg
   return
   }

########################################################################
proc eeManyPlot {colors} {
   loop i 0 [llength $colors] {
	set color [lindex $colors $i]
	if {$i == 0} {
	   set overlay n
	} else {
	   set overlay y
	   }
	eePlot $color $overlay
	}
   return
   }

########################################################################
#Like above, but print out analysis

proc eePrint {icolor} {
   if {[ctype digit $icolor]} {
	set file fft$icolor.fit
   } else {
	set file $icolor
	}

   if {![file exists $file]} {
	error "No file $file!"
	}
   set reg [regReadFromFits $file]

#Get scale factors
   set flux [hdrGetAsDbl $reg.hdr FLUX]
   set pixmm [hdrGetAsDbl $reg.hdr PIXMM]
   set xmm [format %.2f [hdrGetAsDbl $reg.hdr XMM]]
   set ymm [format %.2f [hdrGetAsDbl $reg.hdr YMM]]
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]

#Centroid is in the middle of a pixel
   set row [expr $nrow/2+0.5]
   set col [expr $ncol/2+0.5]
   set radmax [min [expr $nrow/2] [expr $ncol/2]]
#   set chain [annulus $reg $row $col $radmax]

#Try annflux - better flux preservation, less accurate profile.
   set chain [annflux $reg $row $col $radmax]
   set mmax [expr $radmax*$pixmm]
   set format "%6s %7s %10s %15s %10s"
   puts stdout [format $format Step Radius Counts Cumulative Fract]
   set i 0
   set step 1
   set arg1 ""
   if {$arg1 != ""} {set step $arg1}
   set max $flux
   if {$max == 0} {set max 1.}
   loop i 0 [chainSize $chain] $step {
	set hndl [chainGet $chain $i]
	if {$hndl == ""} break
	set radius [format %.4f [expr $i*$pixmm]]
	set flux [format %.2f [exprGet $hndl.cnts]]
	set cum [exprGet $hndl.cum]
	set fract [format %.3f [expr $cum/$max]]
	set cum [format %.2f $cum]
	puts stdout [format $format $i. $radius $flux $cum $fract]
	handleDel $hndl
	if {$fract == 1.0} break
	}
   genericChainDel $chain
   regDel $reg
   return
   }
