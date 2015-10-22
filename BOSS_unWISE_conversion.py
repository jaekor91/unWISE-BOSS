# BOSS_unWISE_conversion.py

import fitsio
import numpy as np
import matplotlib.pyplot as plt
from astrometry.util.util import Tan
from astrometry.util.starutil_numpy import degrees_between

# Already available variables: BOSS RA/DEC. 
# Input: unWISE tile name, channel number
# Output: BOSS RA/DEC within the til
def unWISE2BOSS(tilename, channelNumber, BOSSra, BOSSdec, plot=True):
	# Input tile name
	RA = BOSSra.copy()
	DEC =BOSSdec.copy()
	tileName = tilename
	channel = channelNumber
	# Contructing file address
	fileaddress = get_unwise_filename(tileName, channel)

	tol = 2.00
	if 'm' in tileName:
		ra = float(tileName.split('m')[0])/10.0
		dec = float(tileName.split('m')[1])/10.0
	else: 
		ra = float(tileName.split('p')[0])/10.0
		dec = float(tileName.split('p')[1])/10.0

	iBool = degrees_between(ra, dec, RA,DEC) <tol

	# Getting x,y positions of the objects near by the center of the tile. 
	wcs = Tan(fileaddress)
	ok, x, y = wcs.radec2pixelxy(RA[iBool], DEC[iBool])

	a = np.isnan(x)
	b = np.isnan(y)

	x[a] = 0
	y[b] = 0

	iBool = (1<x)&(x<2048)&(1<y)&(y<2048)
	x=x[iBool]
	y=y[iBool]
	x -= 1
	y -= 1	

	#From the result, I choose to use 1504p196
	if plot:
		objs1 = fitsio.FITS(fileaddress)
		blockImage =objs1[0][:,:] 
		plt.imshow(blockImage, cmap='gray', vmin=-50, vmax=300, origin='lower') # 10/8/2015: Becareful about the orientation of the matrix. 
		plt.scatter(x, y, facecolors='none', edgecolors='r',s=100)
		plt.show() 

	#Return RA/DEC
	return x, y



def get_unwise_filename(tileName, channel):
	return '/project/projectdirs/cosmo/data/unwise/unwise-coadds/'+tileName[0:3]+'/'+tileName+'/unwise-'+tileName+'-'+channel+'-img-m.fits'


# Already available variables: BOSS RA/DEC. 
# Input: a BOSS galaxy RA/DEC and tileID, tileDEC, tileRA
# Output: The name of the tiles that contain the object. Note that there could be multiple tiles.
# To test: 
# import fitsio
# objs2 =fitsio.FITS('/project/projectdirs/cosmo/data/unwise/unwise-coadds/allsky-atlas.fits')
# blockID=objs2[1]['coadd_id'][:] 
# blockRA = objs2[1]['ra'][:] 
# blockDEC= objs2[1]['dec'][:] 
# BOSS2unWISE(150.0, 19.0, blockID, blockRA, blockDEC)
# 

def BOSS2unWISE(ra, dec, tileID, tileRA, tileDEC):
	tol = 1.56
	iBool = degrees_between(ra, dec, tileRA,tileDEC) <tol
	tilesCandidates = tileID[iBool] 
	tiles = []
	for tName in tilesCandidates:
		channel = 'w1'
		fileaddress = get_unwise_filename(tName, channel)
		wcs = Tan(fileaddress)
		ok, x, y = wcs.radec2pixelxy(ra, dec)
		x -= 1
		y -= 1
		iTile= (0<x) &(x<2047)&(0<y)&(y<2047)
		if iTile:
			tiles = np.append(tiles, [tName])

	return tiles


