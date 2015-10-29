# importData.py
# This utility is used for importing the following data:
# BOSS galaxies RA/DEC with appropriate mask
# unWISE tile names


import fitsio
import numpy as np


def importBOSSradec():
	objs = fitsio.FITS('specObj-BOSS-dr10.fits')
	RA = objs[1]['plug_ra'][:]
	DEC = objs[1]['plug_dec'][:]
	redZ = objs[1]['Z'][:]
	iZ = (redZ>0.45)&(redZ<0.7)
	iMASK = (objs[1]['zwarning'][:] == 0) & (objs[1]['class'][:]=="GALAXY") & iZ
	return RA[iMASK], DEC[iMASK]


def importUNWISEtiles():
	objs2 =fitsio.FITS('/project/projectdirs/cosmo/data/unwise/unwise-coadds/allsky-atlas.fits')
	# # The following contains all the RA DEC information. 
	blockID=objs2[1]['coadd_id'][:] 
	blockRA = objs2[1]['ra'][:] 
	blockDEC= objs2[1]['dec'][:] 
	return blockDEC, blockRA, blockID



