#Simulate a PSF from a spot diagram (geometric optics).
#Use same parameters are in psfMap to control image size, resolution.
#	set PIXMM ...
#	set NPIX ...
#However SIGMA (the convolution width of a Gaussian) is in arcsec, not
#pixels, since arcsec is more natural to me at this point.
#
#We allow multiple colors and guarantee that the centroid is at center of
#a pixel (r = nrow/2+0.5, etc) so we can use eePlot and eePrint to analyze.
####################################################################
#Compute PSF map using spot diagram

proc spotMap {hndl xmm ymm colors {stopcheck 1}} {
   global NPIX SCALE PIXMM SIGMA FWHM
   if {[info exists NPIX]} {
	set npix $NPIX
   } else {
	set npix 128
	}

#Telescope params.
   opticInfo $hndl [lindex $colors 0]

#Focal length
   set fl [expr abs([showFocal $hndl [lindex $colors 0] fl])]

#Mirror diameter - this is twice the radius of the aperture stop.
   set diam [telDiam $hndl]

#Final focal ratio
   set fr [expr $fl/$diam]

#Arcsec per pixel.
#This default is nuts - showFocal returns arcsec/mm, not pixel.
   set scale [showFocal $hndl [lindex $colors 0] scale]
   if {[info exists SCALE]} {
	set scale $SCALE
	}

#Pixel size in mm.  Note that the defaults for SCALE and PIXMM are meaningless,
#so one or the other of SCALE and PIXMM should be specified.
#Since "scale" is usually used to mean arcsec/mm, its use here can be
#confusing, so I should be sure to track what I mean.
   set pixmm [expr $scale*$fl/206265.]
   if {[info exists PIXMM]} {
	set pixmm $PIXMM
	set scale [expr $pixmm*206265./$fl]
	}
   set field [expr $scale*$npix]
   set fieldmm [expr $pixmm*$npix]
   echo Field size = [format %.2f $field] arcsec = [format %.3f $fieldmm] mm.
   echo NPIX = $npix
   if {[info exists SIGMA]} {
	}

#Gaussian SIGMA is inputted as arcsec; convert to pixels.
   if {[info exists SIGMA]} {
	echo Convolving with Gaussian, sigma = [format %.4f $SIGMA] arcsec
	set sigma $SIGMA
	set sig [expr $sigma/$scale]
   } elseif {[info exists FWHM]} {
	echo Convolving with Gaussian, FWHM = [format %.4f $FWHM] arcsec
	set sigma [expr $FWHM/2.35]
	set sig [expr $sigma/$scale]
   } else {
	set sig 1.
	set sigma [expr $sig*$scale]
	}
   set reg [regNew $npix $npix]

#Compute weighted centroid
   set wsum 0.
   set xcen 0.
   set ycen 0.
   foreach icolor $colors {
	set weight [showFocal $hndl $icolor weight]

#Assume that a weight of 0 is a weight of 1.
	if {$weight == 0.} {set weight 1.}
	set wsum [expr $wsum+$weight*1.]
	rtrace $hndl $xmm $ymm $icolor $stopcheck
	set xcen [expr $xcen + [exprGet $hndl.diagram->xcen]*$weight]
	set ycen [expr $ycen + [exprGet $hndl.diagram->ycen]*$weight]
	}
   set xcen [expr $xcen/$wsum]
   set ycen [expr $ycen/$wsum]

   foreach icolor $colors {
	rtrace $hndl $xmm $ymm $icolor $stopcheck
	handleSet $hndl.diagram->xcen $xcen
	handleSet $hndl.diagram->ycen $ycen

#I pass length of $filters so we can normalize properly.
	set weight [showFocal $hndl $icolor weight]
	if {$weight == 0.} {set weight 1.}

#Adjust weight so sum over all colors is 1.
	set weight [expr 1.*$weight/$wsum]
	set wgt($icolor) $weight
	diagramMap $hndl $reg $pixmm $sig $weight
	}
   hdrInsertLine $reg.hdr 1 \
"COMMENT PSF Map from CRAY Spot Diagram"
   global env
   hdrInsWithAscii $reg.hdr VERSION \
	[file tail $env(CRAY_DIR)] "CRAY version"
   hdrInsWithAscii $reg.hdr DESIGN [exprGet $hndl.name]
   loop i 0 [llength $colors] {
	set icolor [lindex $colors $i]
	set weight [format %.2f $wgt($icolor)]
	set wave [showFocal $hndl [lindex $colors $i] wave]
	hdrInsWithDbl $reg.hdr WAVE$icolor $wave "Filter index $icolor\
wavelength (microns)"
	hdrInsWithDbl $reg.hdr WEIGHT$icolor $weight "Filter index $icolor\
weight"
	}
   hdrInsWithDbl $reg.hdr SCALE [format %.5f $scale] "Arcsec/pixel"
   hdrInsWithDbl $reg.hdr PIXMM [format %.4f $pixmm] "Pixel size (mm)"
   hdrInsWithDbl $reg.hdr SIGMA [format %.4f $sigma] "Gaussian blur (arcsec)"
   hdrInsWithDbl $reg.hdr XMM $xmm \
	"Target x position (mm) in focal plane"
   hdrInsWithDbl $reg.hdr YMM $ymm \
	"Target y position (mm) in focal plane"


#Location of center - use FITS convention.  We have a square image.
   set pixcen [format %.1f [expr $npix/2.+1.]]
   hdrInsWithDbl $reg.hdr XPIXCEN $pixcen "NAXIS1 center (FITS convention)"
   hdrInsWithDbl $reg.hdr YPIXCEN $pixcen "NAXIS2 center (FITS convention)"

   hdrInsWithDbl $reg.hdr XCEN [format %.4f $xcen] \
	"Actual X position (mm) in focal plane"
   hdrInsWithDbl $reg.hdr YCEN [format %.4f $ycen] \
	"Actual y position (mm) in focal plane"

#Total flux in region.  At the moment, this is hardwired.
   hdrInsWithDbl $reg.hdr FLUX 1000. "Total design flux in region (ADU)"
   regWriteAsFits $reg spot.fit
   regDel $reg
   return
   }

######################################################################
#Plot the latest spot diagram in an OPTIC structure
#Points are plotted relative to centroid
proc diagramMap {hndl reg pixmm sig weight} {
   set nray [handleShow $hndl.diagram->nray]
   if {$nray <=1} return
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]
   set radmax [expr max($nrow,$ncol)*1.]

#radmax doesn't need to be so huge - just go out to 10-sigma
   set radmax [expr min($radmax,10.*$sig)]

#Center in pixels
   set rcen [expr $nrow/2.]
   set ccen [expr $ncol/2.]

#Ahh, eePlot assumes that center is in the middle of a pixel.
#This matches the diffraction fft calculation.
   set rcen [expr $rcen + 0.5]
   set ccen [expr $ccen + 0.5]

#Center in mm.  This has been reset to match the ensemble over all filters.
   set xcen [handleShow $hndl.diagram->xcen]
   set ycen [handleShow $hndl.diagram->ycen]

#flux is integrated intensity of a Gaussian.  I assume 1000. ADU total for the
#image, and adjust flux accordingly to achieve this goal
echo weight $weight nray $nray
   set flux [expr 1000.*$weight/($nray-1.)]

#I delete i=0, because this is the chief ray and is not included in
#xcen, ycen, etc in rtrace.
   for {set i 1} {$i < $nray} {incr i} {
	set xpoint [handleShow $hndl.diagram->xpoint<$i>]
	set ypoint [handleShow $hndl.diagram->ypoint<$i>]

#Does x go left to right?  Let's make it so - that's the column direction.
	set x [expr ($xpoint-$xcen)/$pixmm + $ccen]
	set y [expr ($ypoint-$ycen)/$pixmm + $rcen]

	regGaussAdd $reg $y $x $flux $sig $radmax
	}
   return
   }

