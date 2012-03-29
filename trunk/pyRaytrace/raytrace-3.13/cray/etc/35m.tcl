#################################################################
#Polish off 3.5m design
proc 35mDesign {optic} {

#Arbitrary central wavelength of shortest wave filter.
   set wave0 .44
   setFocal $optic 1 wave .44

#Insert stops
   stopAdd $optic 0 -4838 382 4000
   strutAdd $optic 0 -4838 20 0
   strutAdd $optic 0 -4838 20 90
   strutAdd $optic 0 -4838 20 180
   strutAdd $optic 0 -4838 20 270

#Make field size bigger.
   setFocal $optic 1 xsize 75
   setFocal $optic 1 ysize 75

#Final touches
   colorcount $optic
   stopcomp $optic
   opticinc $optic 1
   return
   }

################################################################
proc 35mOpticPlot {hndl} {
   global XMIN XMAX YMIN YMAX
   set XMIN -5000
   set XMAX 5000
   set YMIN -7000
   set YMAX 3000
   opticPlot $hndl 0 10 1
   scaleSet 0
   return
   }
