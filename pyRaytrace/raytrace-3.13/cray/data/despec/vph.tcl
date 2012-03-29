#Repeat spreadsheet calculations using a TCL routine - easier to 
#add what-if logic.

proc vph {{wmin ""} {wmax ""} {fiber ""} {fin ""} {R ""} {alpha ""}} {

   echo -------------------------------------
   echo Input parameters:

#Min wavelength (microns)
   if {$wmin == ""} {
	set wmin 0.55
	}
   echo Minimum wavelength (microns) wmin $wmin

#Max wavelength (microns)
   if {$wmax == ""} {
	set wmax 1.05
	}
   echo Maximum wavelength (microns) wmax $wmax

#Fiber diameter (microns)
   if {$fiber == ""} {
	set fiber 80.
	}
   echo Fiber diameter (microns) fiber $fiber

#Input f/ratio
   if {$fin == ""} {
	set fin 2.8
	}
   echo Input focal ratio fin $fin

#Nominal Resolution - this is taken to be projected size of fiber diameter.
#The FHWM is actually 0.59 times this.
   if {$R == ""} {
	set R 3000.
	}
   echo Resolution R $R

#Input angle in medium (degrees)
   if {$alpha == ""} {
	set alpha 6.
	}
   echo Bragg angle (deg) alpha $alpha
   echo

#-----------------

#Global parameters
   global Q DETUNE RFWHM
   echo Globals:

#Nominal Q
   if {![info exists Q]} {set Q 10.}
   echo Q $Q

#Detune is used to obtain better balance in efficiencies.
   if {![info exists DETUNE]} {set DETUNE 1.}
   echo DETUNE $DETUNE

#Convert fiber diameter to FWHM - for resolution calculations
   if {![info exists RFWHM]} {set RFWHM .59}
   echo Diameter to FWHM conversion: RFWHM $RFWHM
   echo

#-----------------
#Default inputs

#Detector 1-D lenth (mm)
   set D 61.4

#Number of pixels
   set npix 4096

#Pixel size (micron)
   set pix [expr $D*1.e3/4096.]

#Mean index
   set n 1.5

#------------------

#Derived quantities
   echo Derived parameters:

#Mean wavelength (mm)
   set wmeanmm [expr ($wmin + $wmax)*1.e-3/2.]

#Wavelength interval
   set deltawmm [expr ($wmax-$wmin)*1.e-3]

#Fringe spacing (mm).  Use Bragg condition.
   set bmm [expr $wmeanmm/(2.*$n*sin($alpha/57.3))]
   echo Fringe space b (microns) [format %.2f [expr $bmm*1.e3]]

#Grating lines/mm
   set lines [expr 1./$bmm]
   echo Grating lines/mm lines [format %.1f $lines]

#Grating depth.  Use eqn for Q.
   set dmm [expr $Q*$n*$bmm*$bmm / (2.*3.14*$wmeanmm)]
   echo Grating depth (microns) d [format %.2f [expr $dmm*1.e3]]

#Index semi-amplitude
   set n1 [expr $DETUNE * $wmeanmm*cos($alpha/57.3)/(2.*$dmm)]
   echo Index semi-amplitude n1 [format %.3f $n1]

#Camera focal length
   set Fcam [expr $D*$bmm*cos($alpha/57.3)/$deltawmm]
   echo Camera focal length (mm) Fcam [format %.1f $Fcam]

#-----------------

#Efficiencies.  I do not account for the change in alpha with wavelength.
   set nu1 [expr 3.14*$n1*$dmm/(($wmeanmm-$deltawmm/2.)*cos($alpha/57.3))]
   set xi1 [expr ($deltawmm/2.)*3.14*$dmm/($n*$bmm*$bmm*2.*cos($alpha/57.3))]
   set eff1 [expr pow(sin(sqrt($nu1*$nu1+$xi1*$xi1)),2)/(1.+pow($xi1/$nu1,2))]
   echo Efficiency at short wavelength eff1 [format %.2f $eff1]

   set nu2 [expr 3.14*$n1*$dmm/(($wmeanmm+$deltawmm/2.)*cos($alpha/57.3))]
   set xi2 [expr ($deltawmm/2.)*3.14*$dmm/($n*$bmm*$bmm*2.*cos($alpha/57.3))]
   set eff2 [expr pow(sin(sqrt($nu2*$nu2+$xi2*$xi2)),2)/(1.+pow($xi2/$nu2,2))]
   echo Efficiency at long wavelength eff2 [format %.2f $eff2]

#-----------------

#More quantities

#Fiber diameter in mm
   set fibermm [expr $fiber*1.e-3]

#Pix size in mm
   set pixmm [expr $pix*1.e-3]

#Wavelength per pixel
   set wpixmm [expr $deltawmm/$npix]
   echo Wavelength per pixel (A) wpix [format %.2f [expr $wpixmm*1.e7]]

#Desired resolution in wavelength (mm) - apply to maximum wavelength.
   set rmm [expr ($wmeanmm + ($deltawmm/2.))/$R]
   echo Wavelength resolution (A) r [format %.2f [expr $rmm*1.e7]]

#Resolution in pixels
   set rpix [expr $rmm/$wpixmm]
   echo Resolution in pixels rpix [format %.2f $rpix]

#Demagnification - convert fiber diameter to fiber FWHM
   set demag [expr ($fibermm*$RFWHM)/($pixmm*$rpix)]
   echo Demagnification demag [format %.2f $demag]

#Focal ratio (camera)
   set fout [expr $fin/$demag]
   echo Camera focal ratio fout [format %.2f $fout]

#Collimator focal length
   set Fcoll [expr $Fcam*$demag]
   echo Collimator Focal Length (mm) Fcoll [format %.1f $Fcoll]

#Beam size (mm)
   set beam [expr $Fcam/$fout]
   echo Beam size (mm) beam [format %.2f $beam]

#Field of view
   set fov [expr $D/$Fcam*57.3]
   echo Field of View (deg) fov [format %.1f $fov] diagonal [format %.1f \
	[expr $fov*1.414]]

#Beam blockage for a Schmidt design
   set detarea [expr $D*$D]
   set beamarea [expr 3.14*$beam*$beam/4.]
   set fract [expr $detarea/$beamarea]
   echo Fraction of beam blocked by detector fract [format %.2f $fract]
   echo -------------------------------------
   return
   }
