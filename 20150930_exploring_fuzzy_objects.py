# 9/30/2015
# Exploring the previously noted problem in the unwise-1504p196 tile. Some BOSS objects at lower redshifts do not show.

# Dwonloading the BOSS catalog:
#wget http://data.sdss3.org/sas/dr10/env/SPECTRO_REDUX/specObj-BOSS-dr10.fits

# Using python to explores the BOSS catalog:
source gopython.sh
python

import fitsio
# from fitsio import FITS,FITSHDR
import numpy as np
import matplotlib.pyplot as plt

objs = fitsio.FITS('specObj-BOSS-dr10.fits')

#Plotting RA/DEC
RA = objs[1]['plug_ra'][:]
DEC = objs[1]['plug_dec'][:]

fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax1.scatter(RA[0:1000], DEC[0:1000], s=0.01)
ax1.scatter(RA, DEC, s=0.005,color='black')
plt.axis('equal')
plt.show()

# #Find a near by block
# objs =fitsio.FITS('/project/projectdirs/cosmo/data/unwise/unwise-coadds/allsky-atlas.fits')
# #The following contains all the RA DEC information. 
# blockID=objs[1]['coadd_id'][:] 
# blockRA = objs[1]['ra'][:] 
# blockDEC= objs[1]['dec'][:] 

# blockRA150DEC20 = blockID[(blockRA<150+1.5)&(blockRA>150-1.5)&(blockDEC>20-1.5)&(blockDEC<20+1.5)]


#From the result, I choose to use 1504p196
objs1 = fitsio.FITS('/project/projectdirs/cosmo/data/unwise/unwise-coadds/150/1504p196/unwise-1504p196-w1-img-m.fits')
blockImage =objs1[0][:,:]
plt.imshow(blockImage, cmap='gray', vmin=-50, vmax=300, origin='lower')
plt.show()


iMASK = (objs[1]['zwarning'][:] == 0) & (objs[1]['class'][:]=="GALAXY")
redZ = objs[1]['Z'][:]
#Astrometry: from astrometry.util.util import Tan
from astrometry.util.util import Tan
wcs = Tan('/project/projectdirs/cosmo/data/unwise/unwise-coadds/150/1504p196/unwise-1504p196-w1-img-m.fits')
#Figuering out the block boundary 
xbound = np.array([0,750])
ybound = np.array([0, 750])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.5)&(redZ<0.7)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iMASK & iRADEC & iZ
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-zwarning-GALAXY-red0p5to0p7-1000pixSq.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()



# 9/30/2015: Systematically exploring previously noticed problem.
##########################################################################
# Redshift range 0 to 0.1
# All galaxies and zwarning == 0
xbound = np.array([0,750])
ybound = np.array([0, 750])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.0)&(redZ<0.1)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iMASK & iRADEC & iZ
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-zwarning-GALAXY-red0to0p1-750pixSq.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()





# Redshift range 0.1 to 0.2
# All galaxies and zwarning == 0
xbound = np.array([0,750])
ybound = np.array([0, 750])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.1)&(redZ<0.2)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iMASK & iRADEC & iZ
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-zwarning-GALAXY-red0p1to0p2-750pixSq.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()




# Redshift range 0.2 to 0.4
# All galaxies and zwarning == 0
xbound = np.array([0,750])
ybound = np.array([0, 750])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.2)&(redZ<0.4)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iMASK & iRADEC & iZ
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-zwarning-GALAXY-red0p2to0p4-750pixSq.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()



# Redshift range 0.4 to 0.6
# All galaxies and zwarning == 0
xbound = np.array([0,750])
ybound = np.array([0, 750])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.4)&(redZ<0.6)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iMASK & iRADEC & iZ
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-zwarning-GALAXY-red0p4to0p6-750pixSq.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()





# Redshift range 0.0 to 0.1
# All galaxies with no zwarning == 0
xbound = np.array([0,750])
ybound = np.array([0, 750])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.0)&(redZ<0.1)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iRADEC & iZ & (objs[1]['class'][:]=="GALAXY") #(objs[1]['zwarning'][:] == 0) & 
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-GALAXY-red0to0p1-750pixSq.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()


# Redshift range 0.0 to 0.1
# All types with zwarning == 0
xbound = np.array([0,750])
ybound = np.array([0, 750])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.0)&(redZ<0.1)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iRADEC & iZ & (objs[1]['zwarning'][:] == 0) #&  (objs[1]['class'][:]=="GALAXY") #
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-zwarning-red0to0p1-750pixSq.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()


# Redshift range 0.0 to 0.2
# All types with zwarning == 0
xbound = np.array([0,2047])
ybound = np.array([0, 2047])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.0)&(redZ<0.2)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iRADEC & iZ & (objs[1]['zwarning'][:] == 0) #&  (objs[1]['class'][:]=="GALAXY") #
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-zwarning-red0to0p2.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()




# Redshift range 0.0 to 0.2
# Galaxy types with zwarning == 0
xbound = np.array([0,2047])
ybound = np.array([0, 2047])
ra, dec = wcs.pixelxy2radec(xbound, ybound)
# You don't need to subtract 1.

# Finding BOSS objects in the boundary as well as meet certain mask conditions. 
#Finding BOSS objs in the above frame
iZ = (redZ>0.0)&(redZ<0.2)
iRADEC = (RA>ra[1])&(RA<ra[0])&(DEC>dec[0])&(DEC<dec[1])
iBool = iRADEC & iZ & iMASK
RA150 = RA[iBool]
DEC20 = DEC[iBool]
RA150.size
ok, RA150, DEC20 = wcs.radec2pixelxy(RA150, DEC20)
RA150 -= 1
DEC20 -= 1

#Showing the plot and circles around the objects
implot = plt.imshow(blockImage[xbound[0]:xbound[1],ybound[0]:ybound[1]],cmap='gray', vmin=-50, vmax=200, origin='lower')
plt.scatter(x=RA150, y=DEC20, facecolors='none', edgecolors='r',s=100)
title = "unwise-1504p196-galaxy-zwarning-red0to0p2.eps"
plt.savefig(title, bbox_inches='tight')
plt.show()







# This line of code was used to check the astrometry.
# # #Plotting a cicrlce around a known star. TYC 4099-2033-1 
# centerX = 150.24768
# centerY = 19.88319
# ok, ra, dec = wcs.radec2pixelxy([centerX],[centerY])
# ra -= 1
# dec -= 1
# # dec = 2048-dec+1
# implot = plt.imshow(blockImage,cmap='gray', vmin=-50, vmax=300, origin='lower')
# plt.scatter(x=ra, y=dec, facecolors='none', edgecolors='r',s=100)
# plt.show()


