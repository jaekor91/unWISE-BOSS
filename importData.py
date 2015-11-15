# importData.py
# This utility is used for importing the following data:
# BOSS galaxies RA/DEC with appropriate mask
# unWISE tile names


import fitsio
import numpy as np

def import2MASS():
	objs = fitsio.FITS('2MASS_Reduced_RADEC_XYZ.fits')
	RA = objs[1]['RA_TMASS'][:]
	DEC = objs[1]['DEC_TMASS'][:]
	K_Band = objs[1]['K_TMASS'][:]
	X = objs[1]['X_TMASS'][:]
	Y = objs[1]['Y_TMASS'][:]
	Z = objs[1]['Z_TMASS'][:]
	objs.close()		
	return RA, DEC, K_Band, X, Y, Z

def importBOSSradec():
	objs = fitsio.FITS('specObj-BOSS-dr10.fits')
	RA = objs[1]['plug_ra'][:]
	DEC = objs[1]['plug_dec'][:]
	redZ = objs[1]['Z'][:]
	iZ = (redZ>0.45)&(redZ<0.7)
	iMASK = (objs[1]['zwarning'][:] == 0) & (objs[1]['class'][:]=="GALAXY") & iZ
	objs.close()
	return RA[iMASK], DEC[iMASK]


def importUNWISEtiles():
	objs2 =fitsio.FITS('/project/projectdirs/cosmo/data/unwise/unwise-coadds/allsky-atlas.fits')
	# # The following contains all the RA DEC information. 
	blockID=objs2[1]['coadd_id'][:] 
	blockRA = objs2[1]['ra'][:] 
	blockDEC= objs2[1]['dec'][:]
	objs2.close()
	return blockDEC, blockRA, blockID

def importUNWISEtiles_w_MASK():
	objs2 = fitsio.FITS('allsky-atlas-mask.fits')
	blockID=objs2[1]['coadd_id'][:] 
	blockRA = objs2[1]['ra'][:] 
	blockDEC= objs2[1]['dec'][:]
	blockMASK = objs2[1]['mask'][:]
	objs2.close()
	return blockDEC, blockRA, blockID, blockMASK


