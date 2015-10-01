# unWISE2BOSS.py
# Already available variables: BOSS RA/DEC. 
# Input: unWISE tile name, channel number
# Output: BOSS RA/DEC within the til


def unWISE2BOSS(tilename, channelNumber, BOSSra, BOSSdec):
	# Input tile name
	RA = BOSSra
	DEC =BOSSdec
	tileName = tilename
	channel = channelNumber
	# Contructing file address
	fileaddress = '/project/projectdirs/cosmo/data/unwise/unwise-coadds/'+tileName[0:3]+'/'+tileName+'/unwise-'+tileName+'-'+channel+'-img-m.fits'

	# Getting the RA/DEC bounds corresponding the tile.
	from astrometry.util.util import Tan
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
	ok, RA, DEC = wcs.radec2pixelxy(RA, DEC)
	RA -= 1
	DEC -= 1


	#From the result, I choose to use 1504p196
	objs1 = fitsio.FITS(fileaddress)
	blockImage =objs1[0][:,:]
	plt.imshow(blockImage, cmap='gray', vmin=-50, vmax=300, origin='lower')
	plt.scatter(RA, DEC, facecolors='none', edgecolors='r',s=100)
	plt.show()

	#Return RA/DEC
	return RA, DEC


################################# 
# # This portion of the code is used when developing the code.
# # Loading BOSS RA/DEC here
# import fitsio
# import numpy as np
# import matplotlib.pyplot as plt

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


# # Getting the RA/DEC bounds corresponding the tile.
# from astrometry.util.util import Tan
# wcs = Tan(fileaddress)
# #Figuering out the block boundary 
# xbound = np.array([0,2047])
# ybound = np.array([0, 2047])
# ra, dec = wcs.pixelxy2radec(xbound, ybound)


# # Getting the RA/DEC of objects that lie within the tile boundary.
# iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
# RA = RA[iRADEC]
# DEC = DEC[iRADEC]
# RA.size
# ok, RA, DEC = wcs.radec2pixelxy(RA, DEC)
# RA -= 1
# DEC -= 1


# #From the result, I choose to use 1504p196
# objs1 = fitsio.FITS(fileaddress)
# blockImage =objs1[0][:,:]
# plt.imshow(blockImage, cmap='gray', vmin=-50, vmax=300, origin='lower')
# plt.scatter(RA, DEC, facecolors='none', edgecolors='r',s=100)
# plt.show()





# # # Block information can be found here. Perhaps irrelevant.
# # objs2 =fitsio.FITS('/project/projectdirs/cosmo/data/unwise/unwise-coadds/allsky-atlas.fits')
# # # #The following contains all the RA DEC information. 
# # blockID=objs2[1]['coadd_id'][:] 
# # blockRA = objs2[1]['ra'][:] 
# # blockDEC= objs2[1]['dec'][:] 