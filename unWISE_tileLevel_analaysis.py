# unWISE_tile_level_analysis.py
import BOSS_unWISE_conversion
import numpy as np
import fitsio




# Input: tilename, channelNumber, BOSSra, BOSSdec, pixel number on each side
# Output: Whole image, Image cube, cutout ra and dec positions as x and y.
def unWISE_BOSS_cutouts(tilename, channelNumber, BOSSra, BOSSdec, pixSize):
	RA = BOSSra.copy()
	DEC =BOSSdec.copy()
	tileName = tilename
	channel = channelNumber

	# x, y coordinates.
	# I noticed a potential "unit vector problem" with BOSS_unWISE_conversion but I think
	# as long as the problme is fixed there, it should be OK.
	x, y = BOSS_unWISE_conversion.unWISE2BOSS(tileName, channel, RA, DEC, plot=False)

	# Only selecting galaxies that are away from the boundary 
	# by at least pixSize/2 amount. 
	tol = pixSize/2 # This automatically rounds "down".
	iBool = (x>(tol-1)) & (x<(2047-tol)) & (y>(tol-1)) & (y<(2047-tol))
	cutoutX = x[iBool]
	cutoutY = y[iBool]
	# I might be off by one index but, oh well.

	#### Cutting out square pixels "centered" at the object, stacking them,
	# and returning them as an image cube.
	# Getting the image file address and importing the image.
	fileaddress = BOSS_unWISE_conversion.get_unwise_filename(tileName, channel)
	objs1 = fitsio.FITS(fileaddress)
	blockImage =objs1[0][:,:]
	# blockImage = blockImage[0,::-1]  #modifying the block image by reflecting. 

	# Cutting out and stacking are done here.
	if cutoutX.size > 0:
		imageCube = blockImage[cutoutX[0]-tol:cutoutX[0]+tol, cutoutY[0]-tol:cutoutY[0]+tol] # Not sure if there is away to avoid doing this. 
		for k in range(1,cutoutX.size):
			imageCube = np.dstack((imageCube, blockImage[cutoutX[k]-tol:cutoutX[k]+tol, cutoutY[k]-tol:cutoutY[k]+tol]))
	else:
		imageCube = []

	return blockImage, imageCube, cutoutX, cutoutY