# BOSS_unWISE_conversion.py

import fitsio
import numpy as np
import matplotlib.pyplot as plt
from astrometry.util.util import Tan


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

	# Getting the RA/DEC bounds corresponding to the tile.
	wcs = Tan(fileaddress)
	#Figuering out the block boundary 
	xbound = np.array([0,2047])
	ybound = np.array([0, 2047])
	ra, dec = wcs.pixelxy2radec(xbound, ybound)

	# Getting the RA/DEC of objects that lie within the tile boundary.
	iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
	RA = RA[iRADEC]
	DEC = DEC[iRADEC]
	RA.size
	ok, x, y = wcs.radec2pixelxy(RA, DEC)
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

def BOSS2unWISE(RA, DEC, tileID, tileRA, tileDEC, plot=True):
	tol = 1.56
	angdist = np.arccos(np.sin(tileDEC)*np.sin(DEC)+np.cos(tileDEC)*np.cos(DEC)*np.cos(tileRA-RA))
	iBool = angdist < tol
	tilesCandidates = tileID[iBool]
	tiles = []
	for x in tilesCandidates:
		channel = 'w1'
		tileName = x
		fileaddress = get_unwise_filename(x, channel)
		wcs = Tan(fileaddress)
		xbound = np.array([0,2047])
		ybound = np.array([0, 2047])
		ra, dec = wcs.pixelxy2radec(xbound, ybound)
		iTile = (RA< ra[0]) & (RA>ra[1]) & (DEC>dec[0]) & (DEC <dec[1])
		if iTile:
			tiles = np.append(tiles, [x])	

	return tiles




###################################################################################### 
# This portion of the code is used when developing the code.
# Loading BOSS RA/DEC here
# import fitsio
# import numpy as np
# import matplotlib.pyplot as plt
# from astrometry.util.util import Tan

# objs = fitsio.FITS('specObj-BOSS-dr10.fits')
# RA = objs[1]['plug_ra'][:]
# DEC = objs[1]['plug_dec'][:]
# #For now I am implementing only the basic mask conditions
# iMASK = (objs[1]['zwarning'][:] == 0) & (objs[1]['class'][:]=="GALAXY")
# # redZ = objs[1]['Z'][:]
# # iZ = (redZ>0.5)&(redZ<0.7)
# RA = RA[iMASK]
# DEC = DEC[iMASK]

# # Input tile name
# tileName = '1504p196'
# channel = 'w1'

# # Contructing file address
# fileaddress = '/project/projectdirs/cosmo/data/unwise/unwise-coadds/'+tileName[0:3]+'/'+tileName+'/unwise-'+tileName+'-'+channel+'-img-m.fits'

# # Block information can be found here. Perhaps irrelevant.
# objs2 =fitsio.FITS('/project/projectdirs/cosmo/data/unwise/unwise-coadds/allsky-atlas.fits')
# # # The following contains all the RA DEC information. 
# blockID=objs2[1]['coadd_id'][:] 
# blockRA = objs2[1]['ra'][:] 
# blockDEC= objs2[1]['dec'][:] 


# # Understnading what the tile name means. 
# tileName = '1535p242'
# channel = 'w1'
# fileaddress = '/project/projectdirs/cosmo/data/unwise/unwise-coadds/'+tileName[0:3]+'/'+tileName+'/unwise-'+tileName+'-'+channel+'-img-m.fits'
# wcs = Tan(fileaddress)
# RA = float(tileName[0:3]+'.'+tileName[3])
# DEC = float(tileName[5:7]+'.'+tileName[7])
# ok, x, y = wcs.radec2pixelxy(RA, DEC)
# x -= 1
# y -= 1

# objs1 = fitsio.FITS(fileaddress)
# blockImage =objs1[0][:,:]
# plt.imshow(blockImage, cmap='gray', vmin=-50, vmax=300, origin='lower')
# plt.scatter(x, y, facecolors='none', edgecolors='r',s=100)
# plt.show()
# # so the tile names are pointing close to the center location but not exactly
# # As it appears to be truncated. So perhaps the best thing to do for the function
# # is to find all the tiles that come within 1.56/2 * 2 degree, get their boundaries,
# # and see if the RA/DEC of the object of interest lies within it.
# # Test RA/DEC
# RA = 151.0
# DEC = 20.0
# tol = 1.56
# iBool = (blockRA>(RA-tol))&(blockRA<(RA+tol))&(blockDEC>(DEC-tol))&(blockDEC<(DEC+tol))
# tilesCandidates = blockID[iBool]
# tiles = []
# for x in tilesCandidates:
# 	channel = 'w1'
# 	tileName = x
# 	fileaddress = '/project/projectdirs/cosmo/data/unwise/unwise-coadds/'+tileName[0:3]+'/'+tileName+'/unwise-'+tileName+'-'+channel+'-img-m.fits'
# 	wcs = Tan(fileaddress)
# 	xbound = np.array([0,2047])
# 	ybound = np.array([0, 2047])
# 	ra, dec = wcs.pixelxy2radec(xbound, ybound)
# 	iTile = (RA< ra[0]) & (RA>ra[1]) & (DEC>dec[0]) & (DEC <dec[1])
# 	if iTile:
# 		tiles = np.append(tiles, [x])

# # From this it looks like there is a unique tile for each object ra,dec. 
# # But I will double check it to make sure.
# RA = 151.0
# DEC = np.arange(19,21,0.001)
# tol = 1.56
# for y in DEC:
# 	iBool = (blockRA>(RA-tol))&(blockRA<(RA+tol))&(blockDEC>(y-tol))&(blockDEC<(y+tol))
# 	tilesCandidates = blockID[iBool]
# 	tiles = []
# 	counter = 0
# 	for x in tilesCandidates:
# 		channel = 'w1'
# 		tileName = x
# 		fileaddress = '/project/projectdirs/cosmo/data/unwise/unwise-coadds/'+tileName[0:3]+'/'+tileName+'/unwise-'+tileName+'-'+channel+'-img-m.fits'
# 		wcs = Tan(fileaddress)
# 		xbound = np.array([0,2047])
# 		ybound = np.array([0, 2047])
# 		ra, dec = wcs.pixelxy2radec(xbound, ybound)
# 		iTile = (RA< ra[0]) & (RA>ra[1]) & (y>dec[0]) & (y <dec[1])
# 		if iTile:
# 			counter += 1
# 	if (counter > 1):
# 		print "More than 1 tile for an object."
# 		print RA, y, counter


## Judging from this result. It looks like the boundaries of neighboring tiles overlap.
# More than 1 tile for an object.                                                                                                       
# 151.0 20.418                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.419                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.42                                                                                                                           
# More than 1 tile for an object.                                                                                                       
# 151.0 20.421                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.422                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.423                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.424                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.425                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.426                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.427                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.428                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.429                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.43                                                                                                                           
# More than 1 tile for an object.                                                                                                       
# 151.0 20.431                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.432                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.433                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.434                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.435                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.436                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.437                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.438                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.439                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.44                                                                                                                           
# More than 1 tile for an object.                                                                                                       
# 151.0 20.441                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.442                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.443                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.444                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.445                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.446                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.447                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.448                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.449                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.45                                                                                                                           
# More than 1 tile for an object.                                                                                                       
# 151.0 20.451                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.452                                                                                                                          
# More than 1 tile for an object.                                                                                                       
# 151.0 20.453
# More than 1 tile for an object.
# 151.0 20.454
# More than 1 tile for an object.
# 151.0 20.455
# More than 1 tile for an object.
# 151.0 20.456
# More than 1 tile for an object.
# 151.0 20.457
# More than 1 tile for an object.
# 151.0 20.458
# More than 1 tile for an object.
# 151.0 20.459
# More than 1 tile for an object.
# 151.0 20.46
# More than 1 tile for an object.
# 151.0 20.461
# More than 1 tile for an object.
# 151.0 20.462
# More than 1 tile for an object.
# 151.0 20.463
# More than 1 tile for an object.
# 151.0 20.464
# More than 1 tile for an object.
# 151.0 20.465
# More than 1 tile for an object.
# 151.0 20.466

			




# '1534p378', '1535m243', '1535m319', '1535m349', '1535m743',
#        '1535p242', '1535p318', '1535p348', '1535p742', '1536m137',
#        '1536m197', '1536p136', '1536p196', '1537m425', '1537m652',
#        '1537p424', '1537p651', '1538m107', '1538m182', '1538m273',
#        '1538p106', '1538p181', '1538p272', '1539m440', '1539p439',
#        '1540m076', '1540m167', '1540p075', '1540p166', '1542m016',
#        '1542m031', '1542m288', '1542m409', '1542m682', '1542m773',
#        '1542p000', '1542p015', '1542p030', '1542p287',








