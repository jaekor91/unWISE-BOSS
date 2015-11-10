# unWISE_tile_level_analysis.py
import BOSS_unWISE_conversion
import numpy as np
import fitsio
import matplotlib.pyplot as plt
import math



# Input: tilename, channelNumber, BOSSra, BOSSdec, pixel number on each side
# Output: Whole image, Image cube, cutout ra and dec pix positions, residuals in pix RA/DEC to be added 
def unWISE_BOSS_cutouts(tilename, channelNumber, BOSSra, BOSSdec, pixSize):
	RA = BOSSra
	DEC =BOSSdec
	tileName = tilename
	channel = channelNumber

	# x, y coordinates.
	x, y = BOSS_unWISE_conversion.unWISE2BOSS(tileName, channel, RA, DEC, plot=False)

	# Only selecting galaxies that are away from the boundary 
	# by at least pixSize/2 amount. 
	tol = pixSize/2 # This automatically rounds "down". So if I give 51, tol is 25.
	iBool = (x>(tol)) & (x<(2047-tol)) & (y>(tol)) & (y<(2047-tol)) # Here the indexing might be wrong.
	cutoutX = x[iBool]
	cutoutY = y[iBool]
	# I think the indices are correct here.
	
	intX = np.around(cutoutX)
	intY = np.around(cutoutY)

	diffX = cutoutX - intX
	diffY = cutoutY - intY

	#### Cutting out square pixels "centered" at the object, stacking them,
	# and returning them as an image cube.
	# Getting the image file address and importing the image.
	fileaddress = BOSS_unWISE_conversion.get_unwise_filename(tileName, channel)
	objs1 = fitsio.FITS(fileaddress)
	blockImage =objs1[0][:,:]

	imageCube=np.zeros((2*tol+1, 2*tol+1),dtype=float)
	# Cutting out and stacking are done here.
	if cutoutX.size > 0:
		for X, Y in zip(intX, intY):
			addBlock = blockImage[Y-tol:Y+tol+1, X-tol:X+tol+1]
			# print addBlock.shape #For Debugging purpose.
			imageCube = np.dstack((imageCube, addBlock))
		imageCube = imageCube[:,:,1:] #I am not sure what is the proper way to think about this. I think this is OK because it won't contribute to the sum but there might be a problem with normalization.

	return blockImage, imageCube, cutoutX, cutoutY, diffX, diffY, pixSize




def view_unWISE_cutouts(imageCube, iMin=0, iMax=1): 
# """Given an image cube, this function will plot cutouts as subplots in the range iMin to iMax."""
	cubeHeight = imageCube.shape[2]
	if (iMin > (cubeHeight-1)):
		# if iMin is greater than the cubeHeight, then only show 
		# the last element.
		iMin = (cubeHeight-2)
		iMax = (cubeHeight-1) #This is wrong but I am going with it anyway.
	elif (iMax > (cubeHeight-1)): 
		iMax = (cubeHeight-1)
		# the boundary cases should be better checked. 
	# print iMin, iMax

	# Calculating the number of rows and columns
	numPlots = (iMax-iMin)
	numRows = math.ceil(np.sqrt(numPlots))
	numCol = math.ceil(np.sqrt(numPlots))

	# print numRows, numCol

	pixSize = imageCube.shape[0]
	iCounter=iMin # You exit the while loop before when iCounter reaches iMax
	plotCounter = 1
	while (iCounter < iMax):
		print iCounter, plotCounter
		plt.subplot(numRows, numCol, plotCounter) # numrows, numcols, fignum. (fignum=1, presumably means everything is plotted at the same time.)
		plt.imshow(imageCube[:,:,iCounter], cmap='gray', vmin=-50, vmax=250, interpolation='nearest',origin='lower')
		plt.scatter(pixSize/2, pixSize/2, facecolors='none', edgecolors='r',s=500)
		# plt.scatter(pixSize/2+diffX[iCounter], pixSize/2+diffY[iCounter], facecolors='none', edgecolors='r',s=1000)		
		iCounter = iCounter+1
		plotCounter = plotCounter+1
	
	plt.show()
	return None




def weighted_sum_imageCube(imageCube, weight=0, iMin=0, iMax=1):
	if weight == 0: #That is, if there is no input for weight, then creates an weight image cube with all 1's. 
		weight = np.ones(imageCube.shape, dtype=float)

	weightedCube = imageCube*weight
	summedWeightedCube = np.sum(weightedCube[:,:,iMin:iMax], axis=2)
	normalization = np.sum(weight[:,:,iMin:iMax], axis=2)
	normalizedSummedWeightedCube = summedWeightedCube/normalization

	return normalizedSummedWeightedCube


