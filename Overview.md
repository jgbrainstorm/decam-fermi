_**The general idea of this package is to make it portable and easy to use. Therefore, it does not require any installation and special packages.**_

## dependency ##
Python packages needed:
  * numpy
  * pyfits
  * scipy
  * matplotlib

## Standalone routines ##
(type the command to see the usage)
  * noise        (calculate the noise for a given extension)
  * noisenew     (calculate the noise for a given extension using two images)
  * noiseall     (calculate the noise using one image for all extensions)
  * noiseallnew  (calculate the noise using two images for all extensions)
  * noisechange  (calculate the change of the noise for two images)
  * dark         (calculate dark current for a given extension)
  * darkall      (calculate dark current for all extensions)
  * darkchange   (calculate the change of dark current for two images)
  * sub\_img      (subtract two images)
  * medianImg    (calculate the median image)
  * oscanCCD     (check the overscan)
  * DECamCCD\_def.py  (definition of the CCD and positions)
  * DECamCCD.py      (some defined functions)
  * xymove.py        (send socket command to move the xy-stage to given CCDs)
  * ptc.py           (calculate the photon transfer curve for a given set of images)