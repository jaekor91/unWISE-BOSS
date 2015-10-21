# unWISE_tile_level_analysis.py
import BOSS_unWISE_conversion
import numpy as np
import fitsio




# Input: tilename, channelNumber, BOSSra, BOSSdec, pixel number on each side
# Output: Whole image, Image cube, cutout ra and dec pix positions, residuals in pix RA/DEC to be added 
def unWISE_BOSS_cutouts(tilename, channelNumber, BOSSra, BOSSdec, pixSize):
	RA = BOSSra.copy()
	DEC =BOSSdec.copy()
	tileName = tilename
	channel = channelNumber

	# x, y coordinates.
	x, y = BOSS_unWISE_conversion.unWISE2BOSS(tileName, channel, RA, DEC, plot=False)

	# Only selecting galaxies that are away from the boundary 
	# by at least pixSize/2 amount. 
	tol = pixSize/2 # This automatically rounds "down". So if I give 51, tol is 25.
	iBool = (x>(tol-1)) & (x<(2047-tol)) & (y>(tol-1)) & (y<(2047-tol))
	cutoutX = x[iBool]
	cutoutY = y[iBool]
	# I think the indices are correct here.
	
	intX = cutoutX.astype(int)
	intY = cutoutY.astype(int)

	diffX = cutoutX - intX
	diffY = cutoutY - intY

	#### Cutting out square pixels "centered" at the object, stacking them,
	# and returning them as an image cube.
	# Getting the image file address and importing the image.
	fileaddress = BOSS_unWISE_conversion.get_unwise_filename(tileName, channel)
	objs1 = fitsio.FITS(fileaddress)
	blockImage =objs1[0][:,:]

	imageCube=np.zeros((2*tol, 2*tol),dtype=float)
	# Cutting out and stacking are done here.
	if cutoutX.size > 0:
		for X, Y in zip(intX, intY):
			imageCube = np.dstack((imageCube, blockImage[Y-tol:Y+tol,X-tol:X+tol]))

	imageCube = imageCube[:,:,1:]

	return blockImage, imageCube, cutoutX, cutoutY, diffX, diffY, pixSize











