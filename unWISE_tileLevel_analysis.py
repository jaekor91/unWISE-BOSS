# unWISE_tile_level_analysis.py
import BOSS_unWISE_conversion
import numpy as np
import fitsio
import matplotlib.pyplot as plt
import math
from astrometry.util.util import Tan
from astrometry.util.starutil_numpy import degrees_between


# Input: Tile name (Must have an accompanying mask in the directory)
# Output: Juxtapose the image with the masked image
def view_tile_mask_compare(tName, channel, vmin=-50,vmax=300):
	objs1 = fitsio.FITS(BOSS_unWISE_conversion.get_unwise_filename(tName, channel))
	blockImage =objs1[0][:,:]
	plt.subplot(1,2, 1)
	plt.imshow(blockImage, cmap='gray', vmin=vmin, vmax=vmax, origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 

	objs2 = fitsio.FITS(get_TMASS_mask_filename(tName))
	array = objs2[0][:,:]
	plt.subplot(1,2,2)
	filtered_image = np.copy(blockImage)
	filtered_image[array==0] = np.median(filtered_image)
	plt.imshow(filtered_image, cmap='gray', vmin=vmin, vmax=vmax, origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 
	plt.show() 

	# Turning up the contrast
	plt.imshow(filtered_image, cmap='gray', vmin=filtered_image.min(), vmax=np.percentile(filtered_image,99), origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 
	plt.show()

# INput:unWISE tile name
# Output: the file address for the corresponding mask map based on 2mass
def get_TMASS_mask_filename(tName):
	return '/global/homes/j/jaehyeon/unWISE-BOSS/2MASS_masks'+'/'+tName+'_mask.fits'


# Input: The name of the tile
# Output: None. A plot of the tile. 
def view_tile(tName, channel):
	objs1 = fitsio.FITS(BOSS_unWISE_conversion.get_unwise_filename(tName, channel))
	blockImage =objs1[0][:,:]
	plt.imshow(blockImage, cmap='gray', vmin=-50, vmax=300, origin='lower',interpolation='nearest') # 10/8/2015: Becareful about the orientation of the matrix. 
	plt.show()

# Given RA/DEC returns x,y,z on the unit sphere.
def radec2xyz(ra,dec):
	x = np.cos(ra/180.0*np.pi)*np.cos(dec/180.0*np.pi)
	y = np.sin(ra/180.0*np.pi)*np.cos(dec/180.0*np.pi) # Spotted an error here.
	z = np.sin(dec/180.0*np.pi)
	return x,y,z	

# # Given ra, dec of a tile coordinates and TMASS_X,Y,Z return a boolean array that indicates objects that are nearby.
def iBoolNearby(ra, dec, TMASS_X, TMASS_Y, TMASS_Z, tol=2.0): 
	x,y,z = radec2xyz(ra,dec)
	cosAngle = x*TMASS_X+y*TMASS_Y+z*TMASS_Z
	iBool = cosAngle>np.cos(tol*np.pi/180.0)
	return iBool

# Input: tilename, channelNumber, BOSSra, BOSSdec, pixel number on each side
# Output: Whole image, Image cube, cutout ra and dec pix positions, residuals in pix RA/DEC to be added 
def unWISE_BOSS_cutouts(tilename, channelNumber, BOSSra, BOSSdec, pixSize,TMASS_RA, TMASS_DEC, TMASS_K, TMASS_X, TMASS_Y, TMASS_Z,saveCubes=True, maskSave=True):
	RA = BOSSra
	DEC =BOSSdec
	tileName = tilename
	channel = channelNumber

	# x, y coordinates.
	x, y = BOSS_unWISE_conversion.unWISE2BOSS(tileName, channel, RA, DEC, plot=False)
	if x.size == 0:
		return 

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

	mask = unWISE_mask_map(tileName,channel, TMASS_RA, TMASS_DEC, TMASS_K, TMASS_X, TMASS_Y, TMASS_Z,plot=False, plotsave=maskSave) # The mask will be saved separately.

	imageCube=np.zeros((2*tol+1, 2*tol+1),dtype=float)
	maskCube = np.zeros((2*tol+1, 2*tol+1),dtype=float)

	# Cutting out and stacking are done here.
	if cutoutX.size > 0:
		for X, Y in zip(intX, intY):
			addBlock = blockImage[Y-tol:Y+tol+1, X-tol:X+tol+1]
			addBlock2 = mask[Y-tol:Y+tol+1, X-tol:X+tol+1]
			# print addBlock.shape #For Debugging purpose.
			imageCube = np.dstack((imageCube, addBlock))
			maskCube = np.dstack((maskCube, addBlock2))
		imageCube = imageCube[:,:,1:] #I am not sure what is the proper way to think about this. I think this is OK because it won't contribute to the sum but there might be a problem with normalization.
		maskCube = maskCube[:,:,1:]

	if saveCubes:
		fits = fitsio.FITS(get_imageCube_filename(tileName,channel),'rw', clobber=True) # clobber=True is there to overwrrite the existing file.
		fits.write(imageCube)
		fits.write(maskCube)
		fits.close() 

	return blockImage, mask, imageCube, maskCube, cutoutX, cutoutY, diffX, diffY, pixSize


def get_imageCube_filename(tileName, channel):
	return '/global/homes/j/jaehyeon/unWISE-BOSS/unWISE_imageCubes'+'/'+tileName+'_icubes.fits'



def view_unWISE_cutouts(imageCube, iMin=0, iMax=1, vminPercentile=0, vmaxPercentile=95): 
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
	vmin = np.percentile(imageCube[:,:,:],vminPercentile)
	vmax = np.percentile(imageCube[:,:,:],vmaxPercentile)
	while (iCounter < iMax):
		print iCounter, plotCounter
		plt.subplot(numRows, numCol, plotCounter) # numrows, numcols, fignum. (fignum=1, presumably means everything is plotted at the same time.)
		plt.imshow(imageCube[:,:,iCounter], cmap='gray', vmin=vmin, vmax=vmax, interpolation='nearest',origin='lower')
		plt.scatter(pixSize/2, pixSize/2, facecolors='none', edgecolors='r',s=500)
		# plt.scatter(pixSize/2+diffX[iCounter], pixSize/2+diffY[iCounter], facecolors='none', edgecolors='r',s=1000)		
		iCounter = iCounter+1
		plotCounter = plotCounter+1
	
	plt.show()
	return None




def weighted_sum_imageCube(imageCube, weight=0, iMin=0, iMax=1):
	if type(weight)==int: #That is, if there is no input for weight, then creates an weight image cube with all 1's. 
		weight = np.ones(imageCube.shape, dtype=float)

	weightedCube = imageCube*weight
	summedWeightedCube = np.sum(weightedCube[:,:,iMin:iMax], axis=2)
	normalization = np.sum(weight[:,:,iMin:iMax], axis=2)
	normalizedSummedWeightedCube = summedWeightedCube/normalization

	return normalizedSummedWeightedCube




# 11/11/2015 Addition
def unWISE_mask_map(tileName,channel, TMASS_RA, TMASS_DEC, TMASS_K, TMASS_X, TMASS_Y, TMASS_Z,plot=True, plotsave=True, vmin=-50,vmax=300):
# Input: tileName, channel, and ra,dec,k of tmass objects that are near the tile.
# Output: mask_map is saved in my directory. Also, shows the masked map, 2 2plots.
	fileaddress = BOSS_unWISE_conversion.get_unwise_filename(tileName,channel)

	tol = 2.00 # This is pretty large.
	if 'm' in tileName:
		ra = float(tileName.split('m')[0])/10.0
		dec = float(tileName.split('m')[1])/10.0
	else: 
		ra = float(tileName.split('p')[0])/10.0
		dec = float(tileName.split('p')[1])/10.0

	iBool = iBoolNearby(ra, dec, TMASS_X, TMASS_Y, TMASS_Z, tol=tol)
	tmass_ra = TMASS_RA[iBool]
	tmass_dec = TMASS_DEC[iBool]
	tmass_k = TMASS_K[iBool]

	bins = np.arange(-5,22,1.0)
	inds = np.digitize(tmass_k, bins)	

	n = 2048 # Size of an unWISE tile.
	array = np.ones((n, n),dtype=int)

	rTolerance = np.array([1000, 800, 500, 350, 180, 110, 90, 70, 40,35,30,20,15,10,8,5,4,3.5,2,2,2])*1.25 # 11/15/2015: This change was made to make the mask sizes larger.
	for m in range(3, 24):
		r = rTolerance[m-3]	
		# iBool = degrees_between(ra, dec, tmass_ra[inds==m], tmass_dec[inds==m]) <tol # Not necessary.

		# Getting x,y positions of the objects near by the center of the tile. 
		wcs = Tan(fileaddress)
		ok, x, y = wcs.radec2pixelxy(tmass_ra[inds==m], tmass_dec[inds==m])
		# This is an insurance	
		a = np.isnan(x)
		b = np.isnan(y)
		x[a] = -5000 
		y[b] = -5000 

		ibool = ((-r-1)<x)&(x<(2048+r))&((-r-1)<y)&(y<(2048+r))  # This saves so much time! 
		x -= 1
		y -= 1

		x=x[ibool]
		y=y[ibool]

		for (a,b) in zip (y,x):
			Y,X = np.ogrid[-a:n-a, -b:n-b] # I guess this doesn't really matter.
			mask = X*X + Y*Y <= r*r
			array[mask] = 0

		print m, '(',bins[m-1],',',bins[m],')', r, x.size #Printing the the bin number, mask radius, the size of x.

	if plot:
		# Holes filled with the median image.
		objs1 = fitsio.FITS(fileaddress)
		blockImage =objs1[0][:,:] 		
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

	if plotsave:
		fits = fitsio.FITS(get_TMASS_mask_filename(tileName), 'rw', clobber=True) #
		fits.write(array)
		fits.close()
	return array


# # This function is retired as of 11/15/2015. See nearby function above 
# def tmass_nearby_tile(tileName, TMASS_RA, TMASS_DEC, TMASS_K):
# 	if 'm' in tileName:
# 		ra = float(tileName.split('m')[0])/10.0
# 		dec = float(tileName.split('m')[1])/10.0
# 	else: 
# 		ra = float(tileName.split('p')[0])/10.0
# 		dec = float(tileName.split('p')[1])/10.0
# 	iBool = (TMASS_RA>(ra-2.0))&(TMASS_RA<(ra+2.0))&(TMASS_DEC>(dec-2.0))&(TMASS_DEC<(dec+2.0))
# 	return TMASS_RA[iBool], TMASS_DEC[iBool], TMASS_K[iBool]
