# Mosaicking
mosaic.py is a script to fill an irregular SAOImage ds9 polygon with regularly-spaced mosaic pointings

Key functions written by Taco at the SMA and adapted for use for the CMZoom Survey by Nimesh Patel

Before use:
* On line 122, update the temporary fits file name (Mppix.fits) with the fits file name that you will use to define the mosaic.
* On line 123, update with your primary beam size from 52.0" to your primary beam size
* If desired, also update the beam spacing on line 124.

The script will open the fits file of your choice, then will instruct you to click on the central position of the region, 
followed by a series of points to define the polygonal region to be mapped. It should also work if you open a different 
fits file from the ds9 window that is opened (or modify it, such as adding contours).

This will produce an SAOImage ds9 region file ds9.reg with regularly spaced, primary-beam-sized circular regions.
