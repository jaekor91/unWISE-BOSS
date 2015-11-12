# unWISE_tile_level_analysis.py
import BOSS_unWISE_conversion
import numpy as np
import fitsio
import matplotlib.pyplot as plt
import math
from astrometry.util.util import Tan
from astrometry.util.starutil_numpy import degrees_between

# Input: The name of the tile
# Output: None. A plot of the tile. 
def view_tile(tName, channel):
	objs1 = fitsio.FITS(BOSS_unWISE_conversion.get_unwise_filename(tName, channel))
	blockImage =objs1[0][:,:]
	plt.imshow(blockImage, cmap='gray', vmin=-50, vmax=300, origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 
	plt.show()


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




# 11/11/2015 Addition
def unWISE_mask_map(tileName,channel, tmass_ra, tmass_dec, tmass_k,plot=True, plotsave=True, vmin=-50,vmax=300):
# Input: tileName, channel, and ra,dec,k of tmass objects that are near the tile.
# Output: mask_map is saved in my directory. Also, shows the masked map, 2 2plots.
	bins = np.arange(-5,22,1.0)
	inds = np.digitize(tmass_k, bins)
	fileaddress = BOSS_unWISE_conversion.get_unwise_filename(tileName,channel)
	objs1 = fitsio.FITS(fileaddress)
	blockImage =objs1[0][:,:] 
	rTolerance = np.array([40,35,30,20,15,10,8,5,4,3.5,2,2,2]) # corresponding to m=11 to m=23 (or, [5,6] through [17,18] bins)
	n = 2048 # Size of an unWISE tile.
# As long as r is less than 100, I can only consider objects within the tile.
	array = np.ones((n, n))

	tol = 2.00 # This is pretty large.
	if 'm' in tileName:
		ra = float(tileName.split('m')[0])/10.0
		dec = float(tileName.split('m')[1])/10.0
	else: 
		ra = float(tileName.split('p')[0])/10.0
		dec = float(tileName.split('p')[1])/10.0

	for m in range(11, 24):
		r = rTolerance[m-11]		
		iBool = degrees_between(ra, dec, tmass_ra[inds==m], tmass_dec[inds==m]) <tol

		# Getting x,y positions of the objects near by the center of the tile. 
		wcs = Tan(fileaddress)
		ok, x, y = wcs.radec2pixelxy(tmass_ra[inds==m][np.where(iBool==True)], tmass_dec[inds==m][np.where(iBool==True)])
		# As long as r is less than 100, I can only consider objects within the tile.
		a = np.isnan(x)
		b = np.isnan(y)

		x[a] = 0
		y[b] = 0

		# ibool = (0<x)&(x<500)&(0<y)&(y<500)#
		ibool = (1<x)&(x<2048)&(1<y)&(y<2048)
		x -= 1
		y -= 1

		x=x[ibool]
		y=y[ibool]

		for (a,b) in zip (y,x):
			Y,X = np.ogrid[-a:n-a, -b:n-b] # I guess this doesn't really matter.
			mask = X*X + Y*Y <= r*r
			array[mask] = 0

		print m,r, x.size #Printing the the bin number, mask radius, the size of x.

	if plot:
		# Holes filled with the median image. 
		plt.subplot(1,2, 1)
		plt.imshow(blockImage, cmap='gray', vmin=vmin, vmax=vmax, origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 

		plt.subplot(1,2,2)
		filtered_image = np.copy(blockImage)
		filtered_image[array==0] = np.median(filtered_image)
		plt.imshow(filtered_image, cmap='gray', vmin=vmin, vmax=vmax, origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 
		plt.savefig(tileName+'_masked_compare.eps', bbox_inches='tight',interpolation='nearest')		
		plt.show() 

		# Turning up the contrast
		plt.imshow(filtered_image, cmap='gray', vmin=filtered_image.min(), vmax=np.percentile(filtered_image,99), origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 
		plt.savefig(tileName+'_masked.eps', bbox_inches='tight',interpolation='nearest')		
		plt.show() 

	fits = fitsio.FITS(tileName+'_mask.fits','rw') #This can probably be fixed so that the program is more efficient. Say, save as integer type... I am not sure though.
	fits.write_image(array)
	fits.close()
	return array

