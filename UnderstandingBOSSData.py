# UnderstandingBOSSData
# This is a compliation of commands used to explore 'specObj-BOSS-dr10.fits' data.

# Loading python
source gopython.sh
python

# Importing the file
import fitsio
import numpy as np
import matplotlib.pyplot as plt

objs = fitsio.FITS('specObj-BOSS-dr10.fits')

# Reference http://data.sdss3.org/datamodel/files/SPECTRO_REDUX/specObj.html

# Basic info


######################################################################################################################################################################################################################################################################################################
# 9/29/2015-Jae
# Going through each column
### SURVEY	str	Survey that this object is part of
list(objs[1]['Survey'][:]).count('boss')
# Returns 1515000. All objects come from BOSS.



### INSTRUMENT	str	Instrument that this spectrum was observed with (SDSS or BOSS)
list(objs[1]['Instrument'][:]).count('BOSS')
# Returns 1515000. All objects come from BOSS.



### CHUNK	str	Name of tiling chunk that this spectrum was tiled in (boss1, boss2, etc), important for tracking large-scale structure samples
# Counting the occurence of each type.
from collections import Counter
counts = Counter(objs[1]['chunk'][:])

# Plotting the counts as is. 
import pylab as pl
X = np.arange(len(counts))
pl.bar(X, counts.values(), align='center', width=0.5)
pl.xticks(X, counts.keys())
ymax = max(counts.values()) + 1
pl.ylim(0, ymax)
pl.show()
# The counts are unordered. Also, the labels are hard to see.

# Ordering the dictionary: http://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
import operator
sorted_counts = sorted(counts.items(), key=operator.itemgetter(1)) #0 is sorted by name and 1 is sorted by # of occurences. 
X = np.arange(len(sorted_counts))
xlabels = []
for x in range(34):
    xlabels.append(sorted_counts[x][0])

yvals = []
for y in range(34):
    yvals.append(sorted_counts[y][1])

pl.bar(X, yvals, align='center', width=0.5)
pl.xticks(X, xlabels)
ymax = max(yvals) + 1
pl.ylim(0, ymax)
pl.show()
# You have to zoom into see the xlabels.





### PROGRAMNAME	str	Program within each survey that the plate was part of
from collections import Counter
counts = Counter(objs[1]['PROGRAMNAME'][:])

# Plotting the counts as is. 
import pylab as pl
X = np.arange(len(counts))
pl.bar(X, counts.values(), align='center', width=0.5)
pl.xticks(X, counts.keys())
ymax = max(counts.values()) + 2000
pl.ylim(0, ymax)
pl.show()
# The counts are unordered. Also, the labels are hard to see.




### PLATERUN	str	Drilling run that this plate was drilled in
from collections import Counter
counts = Counter(objs[1]['platerun'][:])

# Plotting the counts as is. 
import pylab as pl
X = np.arange(len(counts))
pl.bar(X, counts.values(), align='center', width=0.5)
pl.xticks(X, counts.keys())
ymax = max(counts.values()) + 2000
pl.ylim(0, ymax)
pl.show()



### PLATEQUALITY	str	Final quality of plate this spectrum came from ("good", "marginal" or "bad")
from collections import Counter
counts = Counter(objs[1]['platequality'][:])

# Plotting the counts as is. 
import pylab as pl
X = np.arange(len(counts))
pl.bar(X, counts.values(), align='center', width=0.5)
pl.xticks(X, counts.keys())
ymax = max(counts.values()) + 2000
pl.ylim(0, ymax)
pl.show()


### PLATESN2	float32	Overall signal-to-noise-squared measure for plate (only for SDSS spectrograph plates)
objs[1]['PLATESN2'][0:50]
min(objs[1]['PLATESN2'][:])
max(objs[1]['PLATESN2'][:])
plt.hist(objs[1]['PLATESN2'][:], bins=np.arange(0,62, 1))
plt.show()
# Not sure if I understand what "only for SDSS spectrograph plates" means.

### DEREDSN2	float32	Dereddened overall signal-to-noise-squared measure for plate (only for BOSS spectrograph plates)
objs[1]['DEREDSN2'][0:50]
min(objs[1]['DEREDSN2'][:])
max(objs[1]['DEREDSN2'][:])
plt.hist(objs[1]['DEREDSN2'][:], bins=np.arange(0,62, 1))
plt.show()
# Not sure if I understand what "only for SDSS spectrograph plates" means.


### LAMBDA_EFF	float32	Effective wavelength drilling position was optimized for (Angstroms)
objs[1]['LAMBDA_EFF'][0:50]
min(objs[1]['LAMBDA_EFF'][:])
max(objs[1]['LAMBDA_EFF'][:])
plt.hist(objs[1]['LAMBDA_EFF'][:], bins=np.arange(3800,5400,100))
plt.show()
# How doe ths drilling position affect the optimal effective wavelength?



### BLUE_FIBER	int32	Set to 1 if this was requested to be a "blue fiber" target, 0 if it was a "red fiber" (in BOSS high redshift LRGs are requested to be on red fibers)
plt.hist(objs[1]['BLUEFIBER'][:], bins=np.arange(-1, 2, 0.25))
plt.show()




### ZOFFSET	float32	Washer thickness used (for backstopping BOSS quasar targets, so they are closer to 4000 Angstrom focal plan (microns)
objs[1]['ZOFFSET'][0:50]
min(objs[1]['ZOFFSET'][:])
max(objs[1]['ZOFFSET'][:])
plt.hist(objs[1]['ZOFFSET'][:], bins=np.arange(0,310, 10))
plt.show()


### SNTURNOFF	float32	Signal to noise measure for MS turnoff stars on plate (-9999 if not appropriate)
objs[1]['SNTURNOFF'][0:50]
min(objs[1]['SNTURNOFF'][:]) #Returns 0
max(objs[1]['SNTURNOFF'][:]) #Returns 0
# What is MS turnoff stars? 


### NTURNOFF	int32	Number of stars used for SNTURNOFF determination
objs[1]['NTURNOFF'][0:50]
min(objs[1]['NTURNOFF'][:]) #Returns 0
max(objs[1]['NTURNOFF'][:]) #Returns 0

### SPECPRIMARY	int1	Set to 1 for primary observation of object, 0 otherwise
plt.hist(objs[1]['SPECPRIMARY'][:], bins=np.arange(-1, 2, 0.25))
plt.show()
# I am not sure what it means that all bits are set to 0. 

###  SPECSDSS	int1	Set to 1 for primary SDSS spectrograph observation of object, 0 otherwise
plt.hist(objs[1]['SPECSDSS'][:], bins=np.arange(-1, 2, 0.25))
plt.show()
# I am not sure what it means that all bits are set to 0. 


### SPECLEGACY	int1	Set to 1 for primary SDSS Legacy program observation of object, 0 otherwise
plt.hist(objs[1]['SPECLEGACY'][:], bins=np.arange(-1, 2, 0.25))
plt.show()
# I am not sure what it means that all bits are set to 0. 


### SPECSEGUE	int1	Set to 1 for primary SDSS SEGUE program observation of object (including SEGUE-1 and SEGUE-2), 0 otherwise
plt.hist(objs[1]['SPECSEGUE'][:], bins=np.arange(-1, 2, 0.25))
plt.show()
# I am not sure what it means that all bits are set to 0. 

### SPECSEGUE1	int1	Set to 1 for primary SDSS SEGUE-1 program observation of object, 0 otherwise
### SPECSEGUE2	int1	Set to 1 for primary SDSS SEGUE-2 program observation of object, 0 otherwise

### SPECBOSS	int1	Set to 1 for primary BOSS spectrograph observation of object, 0 otherwise
plt.hist(objs[1]['SPECBOSS'][:], bins=np.arange(-1, 2, 0.25))
plt.show()
# Returns
# (array([       0.,        0.,        0.,        0.,   116919.,        0.,
#               0.,        0.,  1398081.,        0.,        0.]), array([-1.  , -0.75, -0.5 , -0.25,  0.  ,  0.25,  0.5 ,  0.75,  1.  ,
#         1.25,  1.5 ,  1.75]), <a list of 11 Patch objects>)




### BOSS_SPECOBJ_ID	int32	Identification number internal to BOSS for SPECBOSS=1 objects
plt.scatter(np.arange(len(objs[1]['BOSS_SPECOBJ_ID'][:])), objs[1]['BOSS_SPECOBJ_ID'][:], s=0.1) 
plt.show()
# Looking at the plot, I am not sure how the id's are assigned.

# Here, plotting only those with specboss = 1.
yvals = (objs[1]['BOSS_SPECOBJ_ID'][:])[objs[1]['SPECBOSS'][:]==1]
plt.scatter(np.arange(len(yvals)), yvals, s=0.1) 
plt.show()
# Only plotting ones with specboss = 1 does not help me see if there are duplicates.

# Let's use counter to see if there are any duplicates.
import collections
dups = [item for item, count in collections.Counter(objs[1]['BOSS_SPECOBJ_ID'][:]).items() if count > 1]
# This tells me that there are 103520 objects that are duplicates. I am not sure where they are coming from.
plt.scatter(np.arange(len(dups)), dups, s=0.1) 
plt.show()


# I think for the present purpose, I better just "white-out" the duplicate points by drawing horizontal lines.
plt.scatter(np.arange(len(objs[1]['BOSS_SPECOBJ_ID'][:])), objs[1]['BOSS_SPECOBJ_ID'][:], s=0.1) 
for x in dups[0:1000]:
    plt.axhline(y=x,color='white', linewidth=0.5)
plt.show()
# This turned out to be a disaster. It takes a very long time to do this too.  I am going back to eliminating duplicats and plotting scatter plots approach again.

# Plotting objs[1]['BOSS_SPECOBJ_ID'][:] only if they do not have duplicates.
# For this I will call the original set, A, and the dups, B. 
# I will create a set A-B using the bitwise manipulation. 
# 2**1 corresponds to "in A"
# 2**1 = 2 correspond to "in B"
BitColumn = np.ones(len(objs[1]['BOSS_SPECOBJ_ID'][:]),dtype=int)
for x in range(len(objs[1]['BOSS_SPECOBJ_ID'][:])):
    if objs[1]['BOSS_SPECOBJ_ID'][x] in dups:
        BitColumn[x] +=2
# This is also very slow.... 
yavls = (objs[1]['BOSS_SPECOBJ_ID'][:])[BitColumn<2]
plt.scatter(np.arange(len(yvals)), yvals, s=0.1) 
plt.show()



### SPECOBJID	str	Unique database ID based on PLATE, MJD, FIBERID, RUN2D (same as SkyServer version)
objs[1]['SPECOBJID'][0:50]
import collections
dups = [item for item, count in collections.Counter(objs[1]['SPECOBJID'][:]).items() if count > 1]
plt.scatter(np.arange(len(dups)), dups, s=0.1) 
plt.show()
#This tells me that there are no duplicates.


### FLUXOBJID	str	Unique database ID of flux-based photometric match based on RUN, RERUN, CAMCOl, FIELD, ID (same as SkyServer version)
objs[1]['FLUXOBJID'][0:5]
from collections import Counter
c = Counter(objs[1]['FLUXOBJID'][:])
print c.most_common(100)
# [('                   ', 80097), ('1237663784214069708', 13), ('1237663784750940246', 12), ('1237666408458944620', 12), ('1237666406848397527', 12), ('1237663783677263927', 12), ('1237657070089470153', 12), ('1237663784213938912', 12), ('1237666407922271004', 12), ('1237663783677264888', 12), ('1237666407922270536', 12), ('1237657069553123640', 11), ('1237657070089798657', 11), ('1237666408458944876', 11), ('1237663783139606626', 11), ('1237663783139672221', 11), ('1237663784750153874', 11), ('1237666408458092604', 11), ('1237663783140262862', 11), ('1237663783139475869', 11), ('1237678617434325682', 11), ('1237666407922139693', 11), ('1237663784750875264', 11), ('1237663784749957251', 11), ('1237666406847742556', 11), ('1237663783140131835', 11), ('1237663784213938640', 11), ('1237657587096551994', 11), ('1237663784750743724', 11), ('1237657587095830650', 11), ('1237663238196822833', 11), ('1237657070090125814', 11), ('1237663783140262363', 11), ('1237663783677264698', 11), ('1237657070089994728', 11), ('1237657070089929928', 11), ('1237663783139475790', 11), ('1237663783140262255', 11), ('1237678617434456802', 11), ('1237657587095962118', 11), ('1237663784214135086', 11), ('1237666407922139465', 11), ('1237657070089273745', 11), ('1237663784750482274', 11), ('1237657587096027461', 11), ('1237666406848398320', 11), ('1237663782602736484', 11), ('1237663784750220273', 11), ('1237663783676543671', 11), ('1237666407922139542', 11), ('1237678617433866567', 11), ('1237678617434391038', 11), ('1237663783140262179', 11), ('1237678617434391111', 11), ('1237663784214200901', 11), ('1237663784213873029', 11), ('1237663784750023665', 11), ('1237663782603391823', 11), ('1237678437015618256', 11), ('1237663783139673126', 11), ('1237663783139738006', 11), ('1237663784214070310', 11), ('1237663784214070577', 11), ('1237678617434259473', 11), ('1237663783677199129', 11), ('1237663783677329691', 11), ('1237666408459207639', 11), ('1237657070089994919', 11), ('1237666408459076085', 11), ('1237663784213086826', 11), ('1237678617434456391', 11), ('1237666408459010165', 11), ('1237657587096748512', 11), ('1237663784213545037', 11), ('1237666407922336467', 11), ('1237666408459076351', 11), ('1237666407921746077', 11), ('1237657587095831471', 11), ('1237663784750088644', 11), ('1237663784213217284', 11), ('1237666407922401727', 11), ('1237663784750219491', 11), ('1237663783140392976', 11), ('1237663783676346861', 11), ('1237663784750416193', 11), ('1237663783676215553', 11), ('1237666408458551414', 11), ('1237663784213807898', 11), ('1237663237122883614', 11), ('1237657070090125588', 11), ('1237666408458289409', 11), ('1237657070090060789', 11), ('1237657587096813739', 11), ('1237657070090125967', 11), ('1237663783140327430', 11), ('1237678437015355521', 11), ('1237666406848069830', 11), ('1237678617433670041', 11), ('1237666407922270593', 11), ('1237663784213086320', 11)]
# I am not sure what the output means but there are a bunch of objects who do not have a unique id.


### BESTOBJID	str	Unique database ID of (recommended) position-based photometric match based on RUN, RERUN, CAMCOl, FIELD, ID (same as SkyServer version)
objs[1]['BESTOBJID'][0:50]
from collections import Counter
c = Counter(objs[1]['BESTOBJID'][:])
print c.most_common(100)
# [('                   ', 147210), ('1237663784214069708', 13), ('1237666408458944620', 12), ('1237666406848397527', 12), ('1237663783677263927', 12), ('1237657070089470153', 12), ('1237663784213938912', 12), ('1237666407922271004', 12), ('1237663783677264888', 12), ('1237666407922270536', 12), ('1237657069553123640', 11), ('1237657070089798657', 11), ('1237666408458944876', 11), ('1237663783139606626', 11), ('1237663783139672221', 11), ('1237663784750153874', 11), ('1237666408458092604', 11), ('1237663783140262862', 11), ('1237663783139475869', 11), ('1237678617434325682', 11), ('1237666407922139693', 11), ('1237663784750875264', 11), ('1237663784749957251', 11), ('1237666406847742556', 11), ('1237663783140131835', 11), ('1237663784213938640', 11), ('1237657587096551994', 11), ('1237663784750743724', 11), ('1237657587095830650', 11), ('1237663238196822833', 11), ('1237657070090125814', 11), ('1237663783140262363', 11), ('1237663783677264698', 11), ('1237657070089994728', 11), ('1237657070089929928', 11), ('1237663783139475790', 11), ('1237663783140262255', 11), ('1237678617434456802', 11), ('1237657587095962118', 11), ('1237663784214135086', 11), ('1237666407922139465', 11), ('1237657070089273745', 11), ('1237663784750482274', 11), ('1237666406848398320', 11), ('1237663784750220273', 11), ('1237663783676543671', 11), ('1237666407922139542', 11), ('1237678617433866567', 11), ('1237678617434391038', 11), ('1237663783140262179', 11), ('1237678617434391111', 11), ('1237663784750940249', 11), ('1237663784214200901', 11), ('1237663784213873029', 11), ('1237663784750023665', 11), ('1237663782603391823', 11), ('1237678437015618256', 11), ('1237663783139673126', 11), ('1237657587096027462', 11), ('1237663783139738006', 11), ('1237663784214070310', 11), ('1237663784214070577', 11), ('1237678617434259473', 11), ('1237663783677199129', 11), ('1237663783677329691', 11), ('1237666408459207639', 11), ('1237657070089994919', 11), ('1237666408459076085', 11), ('1237663784213086826', 11), ('1237678617434456391', 11), ('1237666408459010165', 11), ('1237657587096748512', 11), ('1237663784213545037', 11), ('1237666407922336467', 11), ('1237666408459076351', 11), ('1237666407921746077', 11), ('1237663784750088644', 11), ('1237663784213217284', 11), ('1237666407922401727', 11), ('1237663784750219491', 11), ('1237663783140392976', 11), ('1237663783676346861', 11), ('1237663784750416193', 11), ('1237663783676215553', 11), ('1237666408458551414', 11), ('1237663784213807898', 11), ('1237663237122883614', 11), ('1237657070090125588', 11), ('1237666408458289409', 11), ('1237657070090060789', 11), ('1237657587096813739', 11), ('1237657070090125967', 11), ('1237663783140327430', 11), ('1237678437015355521', 11), ('1237666406848069830', 11), ('1237678617433670041', 11), ('1237666407922270593', 11), ('1237663784213086320', 11), ('1237678617433800942', 11), ('1237666408458223652', 11)]
# >>>
# I am not sure what the output means but there are a bunch of objects who do not have a unique id.



### TARGETOBJID	str	Unique database ID of targeting object based on RUN, RERUN, CAMCOl, FIELD, ID (same as SkyServer version)
objs[1]['TARGETOBJID'][0:50]
from collections import Counter
c = Counter(objs[1]['TARGETOBJID'][:])
print c.most_common(100)
# [('                      ', 159031), ('   1188107536708731266', 1240), ('   1190370396063730577', 494), ('   1188954766252704783', 460), ('   1188108984649581267', 435), ('   1188111948177016197', 193), ('   1237666407922271004', 12), ('   1164192265454420047', 12), ('   1237666408458944620', 12), ('   1237666407922139693', 11), ('   1237666406848004112', 11), ('   1237663784214135086', 11), ('   1237663783140131835', 11), ('   1164476262125863808', 11), ('   1237666406847807721', 11), ('   1237657587096027462', 11), ('   1237657070089929928', 11), ('   1237663783140327430', 11), ('   1237663783139475790', 11), ('   1237666407922139465', 11), ('   1237657070089994728', 11), ('   1237663783140262255', 11), ('   1164192265454354627', 11), ('   1237657070090060789', 11), ('   1237663783676543671', 11), ('   1237663237122883614', 11), ('   1237663784213938912', 11), ('   1237663783677264857', 11), ('   1237663783676346861', 11), ('   1164192265990504635', 11), ('   1164476260514922695', 11), ('   1237663783140262363', 11), ('   1237678617434456802', 11), ('   1237663784750875264', 11), ('   1164192265453700086', 11), ('   1164476262125666824', 11), ('   1237657070090125588', 11), ('   1237657587096748512', 11), ('   1164476262125535338', 11), ('   1237663782602736484', 11), ('   1237663783139673126', 11), ('   1237666407922139542', 11), ('   1237663783677199129', 11), ('   1164192264379957794', 11), ('   1237663782603391823', 11), ('   1237666408458551414', 11), ('   1237663784750154004', 11), ('   1164192264379957491', 11), ('   1237663783140262179', 11), ('   1237663784214070310', 11), ('   1237678617434456391', 11), ('   1237663784214070577', 11), ('   1237657070089798139', 11), ('   1237666407922074678', 11), ('   1237666408459207639', 11), ('   1237663783139606626', 11), ('   1237663784750416193', 11), ('   1237663784213086320', 11), ('   1237663784213807898', 11), ('   1237663783676937053', 11), ('   1237657069553123640', 11), ('   1237678617434391038', 11), ('   1237657070089798657', 11), ('   1237663784214069708', 11), ('   1237678617434391111', 11), ('   1237657587095962118', 11), ('   1237663784750678945', 11), ('   1237663783677329691', 11), ('   1164192263843152035', 11), ('   1237666406847742556', 11), ('   1237663783677264698', 11), ('   1237663784213086826', 11), ('   1164476261052711152', 11), ('   1237663784213284115', 11), ('   1237666407922401727', 11), ('   1237663783140392976', 11), ('   1237663783676215553', 11), ('   1237663784214200901', 11), ('   1164476261589319938', 11), ('   1237663784213873029', 11), ('   1164476260515643648', 11), ('   1237657070089994919', 11), ('   1237663784213217633', 11), ('   1164192265453633826', 11), ('   1164476262125404441', 11), ('   1237666406848398320', 11), ('   1237657070089339608', 11), ('   1237657587095830650', 11), ('   1237666408459076085', 11), ('   1237666408458092604', 11), ('   1237663784213807662', 11), ('   1237657070089208330', 11), ('   1237678437015618256', 11), ('   1237666408459076351', 11), ('   1237657070090125814', 11), ('   1164192265453633590', 11), ('   1237663783139738006', 11), ('   1237678617433735806', 11), ('   1237657587096813835', 11), ('   1237666407921615135', 11)]



### PLATEID	str	Unique database ID of plate based on PLATE, MJD, RUN2D (same as SkyServer version)
objs[1]['PLATEID'][0:50]
from collections import Counter
c = Counter(objs[1]['PLATEID'][:])
print c.most_common(100)
# [('5648639928516878336', 1000), ('5103704375953858560', 1000), ('4512606916489650176', 1000), ('4343721929070747648', 1000), ('6068600598584238080', 1000), ('4756927197935443968', 1000), ('6060719299420889088', 1000), ('5756726323080208384', 1000), ('4711891201896620032', 1000), ('5045157577023168512', 1000), ('5826532117824544768', 1000), ('5225301567436365824', 1000), ('6117014294226149376', 1000), ('4580160912393379840', 1000), ('5180265565156417536', 1000), ('4532873114292723712', 1000), ('5770237120938909696', 1000), ('5210664863597469696', 1000), ('6073104198194831360', 1000), ('5505650640784072704', 1000), ('5782622019461193728', 1000), ('6075355997102546944', 1000), ('5843420618474004480', 1000), ('6368089972764188672', 1000), ('6502072062836088832', 1000), ('5769111221552160768', 1000), ('4537376719104253952', 1000), ('6791428339867721728', 1000), ('5404319648044163072', 1000), ('5051912980423647232', 1000), ('5972899106200625152', 1000), ('5216294362695475200', 1000), ('5334513853417267200', 1000), ('6749770043750752256', 1000), ('5772488921188802560', 1000), ('4795207795539845120', 1000), ('5412200947610165248', 1000), ('4336966529126375424', 1000), ('5704934924462989312', 1000), ('4558768816612843520', 1000), ('4360610427740495872', 1000), ('6053963899040309248', 1000), ('5193776364072083456', 1000), ('6172183390265417728', 1000), ('5596848533338988544', 1000), ('5461740543393800192', 1000), ('4913427291697455104', 1000), ('4535124919072464896', 1000), ('4321203936218718208', 1000), ('6820701737831505920', 1000), ('5299610956355477504', 1000), ('5219672067868598272', 1000), ('5491013941424693248', 1000), ('5752222722915966976', 1000), ('4451808326955966464', 1000), ('4606056611139952640', 1000), ('4763682597191950336', 1000), ('5071053273890693120', 1000), ('4912301386455457792', 1000), ('4597049410928910336', 1000), ('4448430627252215808', 1000), ('6139532292698546176', 1000), ('6182316489712214016', 1000), ('6821827637788680192', 1000), ('4536250818928975872', 1000), ('5305240457315753984', 1000), ('6012305603393101824', 1000), ('4937071188013096960', 1000), ('5200531763647356928', 1000), ('4166955648175972352', 1000), ('5961640107031535616', 1000), ('4468696819300704256', 1000), ('5321003055575343104', 1000), ('6359082773224235008', 1000), ('6022438701178953728', 1000), ('5970647306319831040', 1000), ('5778118420202921984', 1000), ('5060920174645223424', 1000), ('4472074523752407040', 1000), ('4609434310575267840', 1000), ('4901042387303145472', 1000), ('6122643794112684032', 1000), ('6826331237499936768', 1000), ('4147815344827146240', 1000), ('5188146870946766848', 1000), ('5460614647664484352', 1000), ('4983233090116001792', 1000), ('4530621319864524800', 1000), ('6035949500799262720', 1000), ('6745266443838169088', 1000), ('4746794099142959104', 1000), ('5364913151304671232', 1000), ('4839117892074479616', 1000), ('4339218329007169536', 1000), ('4428164423124131840', 1000), ('6115888394772291584', 1000), ('4908923686768484352', 1000), ('4268286639137497088', 1000), ('5662150727986192384', 1000), ('5058668375150305280', 1000)]
# 1000 fibers in each.



# 9/30/2015
######################################################################################################################################################################################################################################################################################################
### NSPECOBS	int16	Number of spectroscopic observations of this source
objs[1]['nspecobs'][0:50]
min(objs[1]['nspecobs'][:])
max(objs[1]['nspecobs'][:])
plt.hist(objs[1]['nspecobs'][:], bins=np.arange(0, 14, 0.5))
plt.show()
# Most objects have 1 or 2 spectra.


### FIRSTRELEASE	str	Name of first release this PLATE, MJD, FIBERID, RUN2D was associated with
### RUN2D	str	Spectroscopic 2D reduction (extraction of spectra) name
### RUN1D	str	Spectroscopic 1D reduction (redshift and classification) name
### DESIGNID	int32	Design identification number for plate
### CX	float64	Position of object on J2000 unit sphere
### CY	float64	Position of object on J2000 unit sphere
### CZ	float64	Position of object on J2000 unit sphere
### XFOCAL	float32	Hole position on plate (+X = +RA) in mm
### YFOCAL	float32	Hole position on plate (+Y = +DEC) in mm


### SOURCETYPE	str	String expressing type of source (similar to OBJTYPE in DR8 and earlier)
objs[1]['sourcetype'][0:50]
from collections import Counter
counts = Counter(objs[1]['sourcetype'][:])

import operator
sorted_counts = sorted(counts.items(), key=operator.itemgetter(1)) #0 is sorted by name and 1 is sorted by # of occurences. 
xlabels = []
yvals = []
for y in range(len(sorted_counts)):
    if sorted_counts[y][1] > 10000:
        yvals.append(sorted_counts[y][1])    
        xlabels.append(sorted_counts[y][0]) 

X = np.arange(len(yvals))
pl.bar(X, yvals, align='center', width=0.5)
pl.xticks(X, xlabels)
ymax = max(yvals) + 1
pl.ylim(0, ymax)
pl.show()


### TARGETTYPE	str	General type of target ("SCIENCE", "STANDARD" or "SKY")
objs[1]['TARGETTYPE'][0:50]
from collections import Counter
counts = Counter(objs[1]['TARGETTYPE'][:])

import operator
sorted_counts = sorted(counts.items(), key=operator.itemgetter(1)) #0 is sorted by name and 1 is sorted by # of occurences. 
xlabels = []
yvals = []
for y in range(len(sorted_counts)):
    if sorted_counts[y][1] > 10000:
        yvals.append(sorted_counts[y][1])    
        xlabels.append(sorted_counts[y][0]) 

X = np.arange(len(yvals))
pl.bar(X, yvals, align='center', width=0.5)
pl.xticks(X, xlabels)
ymax = max(yvals) + 1
pl.ylim(0, ymax)
pl.show()


### PRIMTARGET	int32	Deprecated version of primary (science) target flags (meanings highly overloaded)
### SECTARGET	int32	Deprecated version of secondary (calibration) target flags (meanings highly overloaded)
### LEGACY_TARGET1	int32	Primary (science) target flags for SDSS-I and SDSS-II Legacy survey
### LEGACY_TARGET2	int32	Secondary (calibration) target flags for SDSS-I and SDSS-II Legacy survey
### SPECIAL_TARGET1	int32	Primary (science) target flags for SDSS-I and SDSS-II special program targets
### SPECIAL_TARGET2	int32	Secondary (calibration) target flags for SDSS-I and SDSS-II special program targets
### SEGUE1_TARGET1	int32	Primary (science) target flags for SEGUE-1 targets
### SEGUE2_TARGET1	int32	Primary (science) target flags for SEGUE-2 targets
### SEGUE2_TARGET2	int32	Secondary (calibration) target flags for SEGUE-2 targets
### MARVELS_TARGET1	int32	Primary (science) target flags for MARVELS targets
### MARVELS_TARGET2	int32	Secondary (calibration) target flags for MARVELS targets
### BOSS_TARGET1	int64	Primary (science) target flags for BOSS targets
objs[1]['BOSS_TARGET1'][0:50]
min(objs[1]['BOSS_TARGET1'][:])
max(objs[1]['BOSS_TARGET1'][:])
plt.hist(objs[1]['BOSS_TARGET1'][:], bins=np.arange(0, 2**15, 2**12))
plt.show()
# I am not sure what these flags mean.

#BOSS_TARGET2	int64	Always set to zero (placeholder for BOSS target flags never used)
### ANCILLARY_TARGET1	int64	Target flags for BOSS ancillary targets
### ANCILLARY_TARGET2	int64	More target flags for BOSS ancillary targets
### SPECTROGRAPHID	int16	Which spectrograph (1 or 2)
### PLATE	int32	Plate number (each plate corresponds to an actual plug plate)
### TILE	int32	Tile number (each tile can have several plates drilled for it)
### MJD	int32	Modified Julian Day of observation
### FIBERID	int32	Fiber number



### OBJID	int32[5]	SDSS photometric object identification numbers (RUN, RERUN, CAMCOL, FIELD, ID)
objs[1]['OBJID'][0:50]


### PLUG_RA	float64	Right ascension of hole (J2000 deg)
### PLUG_DEC	float64	Declination of hole (J2000 deg)
#RA/DEC 
RA = objs[1]['plug_ra'][:]
DEC = objs[1]['plug_dec'][:]

#Plotting RA/DEC of all objects. 
fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax1.scatter(RA[0:1000], DEC[0:1000], s=0.01)
ax1.scatter(RA, DEC, s=0.005,color='black')
plt.axis('equal')
plt.show()



### CLASS	str	Best spectroscopic classification ("STAR", "GALAXY" or "QSO")
from collections import Counter
counts = Counter(objs[1]['CLASS'][:])

import operator
sorted_counts = sorted(counts.items(), key=operator.itemgetter(1)) #0 is sorted by name and 1 is sorted by # of occurences. 
xlabels = []
yvals = []
for y in range(len(sorted_counts)):
    if sorted_counts[y][1] > 10000:
        yvals.append(sorted_counts[y][1])    
        xlabels.append(sorted_counts[y][0]) 

X = np.arange(len(yvals))
pl.bar(X, yvals, align='center', width=0.5)
pl.xticks(X, xlabels)
ymax = max(yvals) + 1
pl.ylim(0, ymax)
pl.show()


#Classes in the catalog
objs[1]['class'][0:50] #To see examples. 
list(objs[1]['class'][:]).count('QSO   ') #318015
list(objs[1]['class'][:]).count('STAR  ') #203898
list(objs[1]['class'][:]).count('GALAXY') #993087

#Different objects and their RA/DEC distribution
# ibool = (objs[1]['class'][:]=='QSO   ')
ibool = (objs[1]['class'][:]=='STAR  ')
# ibool = (objs[1]['class'][:]=='GALAXY')
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(RA[ibool], DEC[ibool], s=0.005,color='black')
plt.axis('equal')
plt.show()
#Conclusion: It appears that stars, galaxies, and QSOs cover about the same region.




### SUBCLASS	str	Best spectroscopic subclassification
from collections import Counter
counts = Counter(objs[1]['subclass'][:])
import operator
sorted_counts = sorted(counts.items(), key=operator.itemgetter(1)) #0 is sorted by name and 1 is sorted by # of occurences. 
xlabels = []
yvals = []
for y in range(len(sorted_counts)):
    if sorted_counts[y][1] > 10000:
        yvals.append(sorted_counts[y][1])    
        xlabels.append(sorted_counts[y][0]) 

X = np.arange(len(yvals))
pl.bar(X, yvals, align='center', width=0.5)
pl.xticks(X, xlabels)
ymax = max(yvals) + 1
pl.ylim(0, ymax)
pl.show()



### Z	float32	Best redshift
#Different objects and their redshift distributions. 
objs[1]['z'][0:50]
min(objs[1]['z'][:])
max(objs[1]['z'][:])
plt.hist(objs[1]['z'][:], bins=np.arange(-0.5, 7.5, 0.1))
plt.show()

# Different objects and their color distribution. 
# Options: Class ("STAR", "GALAXY" or "QSO"), targettype ("SCIENCE", "STANDARD" or "SKY")
objs[1]['targettype'][0:50]
objs[1]['class'][0:50]

# ibool = (objs[1]['class'][:]=='STAR  ') #Centered around 0. Cuspy distribution.
# ibool = (objs[1]['class'][:]=="QSO   ")  # Pretty widely distributed from 0.002 to 7.02.
# ibool = (objs[1]['class'][:]=="GALAXY")
# ibool = (objs[1]['targettype'][:]=="SCIENCE ")
# ibool = (objs[1]['targettype'][:]=="SKY     ")
# ibool = (objs[1]['targettype'][:]=="STANDARD")
# ibool = (objs[1]['class'][:]=='STAR  ') & (objs[1]['ZWARNING'][:] ==0)
# ibool = (objs[1]['class'][:]=="QSO   ") & (objs[1]['ZWARNING'][:] ==0)
# ibool = (objs[1]['class'][:]=="GALAXY") & (objs[1]['ZWARNING'][:] ==0) 
# ibool = (objs[1]['targettype'][:]=="SCIENCE ") & (objs[1]['ZWARNING'][:] ==0)
ibool = (objs[1]['targettype'][:]=="SKY     ") & (objs[1]['ZWARNING'][:] ==0) # No objects. zwarning pretty much eliminates them all.
# ibool = (objs[1]['targettype'][:]=="STANDARD")
((objs[1]['z'][:])[ibool])[0:50]
min((objs[1]['z'][:])[ibool])
max((objs[1]['z'][:])[ibool])
plt.hist((objs[1]['z'][:])[ibool], bins=np.arange(min((objs[1]['z'][:])[ibool])-0.1, max((objs[1]['z'][:])[ibool])+0.1, 0.01))
plt.savefig('z-TARGETTYPE-SCIENCE-ZWARNING.eps', bbox_inches='tight')
plt.show()


### Z_ERR	float32	Error in best redshift



### RCHI2	float32	Reduced chi-squared of best fit
### DOF	int32	Number of degrees of freedom in best fit
### RCHI2DIFF	float32	Difference in reduced chi-squared between best and second best fit
### TFILE	str	File that best fit template comes from in idlspec2d product
### TCOLUMN	float32[10]	Columns of template files that correspond to each template
### NPOLY	int32	Number of polynomial terms in fit
### THETA	float32[10]	Template coefficients of best fit
### VDISP	float32	Velocity dispersion (km/s)
# VDISP_ERR	float32	Error in velocity dispersion (km/s)
# VDISPZ	float32	Redshift associated with best-fit velocity dispersion
# VDISPZ_ERR	float32	Error in redshift associated with best-fit velocity dispersion
# VDISP_DOF	int32	Number of degrees of freedom in velocity dispersion fit
# VDISPCHI2	float32	Chi-squared for best-fit velocity dispersion
# VDISPNPIX	int32	Number of pixels overlapping the templates used in the velocity dispersion fit



### WAVEMIN	float32	Minimum observed (vacuum) wavelength (Angstroms)
### WAVEMAX	float32	Maximum observed (vacuum) wavelength (Angstroms)




# WCOVERAGE	float32	Coverage in wavelength, in units of log10 wavelength
ZWARNING	int32	Bitmask of spectroscopic warning values; 0 means everything is OK
SN_MEDIAN_ALL	float32	Median signal-to-noise per pixel across full spectrum
SN_MEDIAN	float32[5]	Median signal-to-noise per pixel within each of the ugriz bandpasses
# CHI68P	float32	68-th percentile value of abs(chi) of the best-fit synthetic spectrum to the actual spectrum (around 1.0 for a good fit)
# FRACNSIGMA	float32[10]	Fraction of pixels deviant by more than N sigma relative to best-fit (for 1,2,..,10 sigma)
# FRACNSIGHI	float32	Fraction of pixels high by more than N sigma relative to best-fit (for 1,2,..,10 sigma)
# FRACNSIGLO	float32	Fraction of pixels low by more than N sigma relative to best-fit (for 1,2,..,10 sigma)
# SPECTROFLUX	float32[5]	Spectrum projected onto ugriz filters (nanomaggies)
SPECTROFLUX_IVAR	float32[5]	Inverse variance of spectrum projected onto ugriz filters (nanomaggies)
SPECTROSYNFLUX	float32[5]	Best-fit template spectrum projected onto ugriz filters (nanomaggies)
SPECTROSYNFLUX_IVAR	float32[5]	Inverse variance of best-fit template spectrum projected onto ugriz filters (nanomaggies)
SPECTROSKYFLUX	float32[5]	Sky flux in each of the ugriz imaging filters (nanomaggies)
ANYANDMASK	int32	For each bit, records whether any pixel in the spectrum has that bit set in its ANDMASK
ANYORMASK	int32	For each bit, records whether any pixel in the spectrum has that bit set in its ORMASK
# SPEC1_G	float32	Signal-to-noise squared for spectrograph #1, at g=20.20 for SDSS spectrograph spectra, g=21.20 for BOSS spectrograph spectra
# SPEC1_R	float32	Signal-to-noise squared for spectrograph #1, at r=20.25 for SDSS spectrograph spectra, r=20.20 for BOSS spectrograph spectra
# SPEC1_I	float32	Signal-to-noise squared for spectrograph #1, at i=19.90 for SDSS spectrograph spectra, i=20.20 for BOSS spectrograph spectra
# SPEC2_G	float32	Signal-to-noise squared for spectrograph #2, at g=20.20 for SDSS spectrograph spectra, g=21.20 for BOSS spectrograph spectra
# SPEC2_R	float32	Signal-to-noise squared for spectrograph #2, at r=20.25 for SDSS spectrograph spectra, r=20.20 for BOSS spectrograph spectra
# SPEC2_I	float32	Signal-to-noise squared for spectrograph #2, at i=19.90 for SDSS spectrograph spectra, i=20.20 for BOSS spectrograph spectra
# ELODIE_FILENAME	str	File name for best-fit ELODIE star
# ELODIE_OBJECT	str	Star name for ELODIE star
# ELODIE_SPTYPE	str	ELODIE star spectral type
# ELODIE_BV	float32	(B-V) color index for ELODIE star (mag)
# ELODIE_TEFF	float32	Effective temperature of ELODIE star (Kelvin)
# ELODIE_LOGG	float32	log10(gravity) of ELODIE star
# ELODIE_FEH	float32	Metallicity [Fe/H] of ELODIE star
# ELODIE_Z	float32	Redshift fit to ELODIE star
# # ELODIE_Z_ERR	float32	Error in redshift fit to ELODIE star
# ELODIE_Z_MODELERR	float32	Standard deviation in redshift among the 12 best-fit ELODIE stars
# ELODIE_RCHI2	float32	Reduced chi-squared of fit to best ELODIE star
# ELODIE_DOF	int32	Degrees of freedom in fit to best ELODIE star
# Z_NOQSO	float32	Best redshift when ignoring QSO fits, recommended for BOSS CMASS and LOWZ targets; calculated only for survey='boss' spectra, not for any SDSS spectrograph data
# Z_ERR_NOQSO	float32	Error in Z_NOQSO redshift
# ZWARNING_NOQSO	int32	For Z_NOQSO redshift, the bitmask of spectroscopic warning values; 0 means everything is OK
# CLASS_NOQSO	str	Spectroscopic classification for Z_NOQSO redshift
# SUBCLASS_NOQSO	str	Spectroscopic subclassification for Z_NOQSO redshift
# RCHI2DIFF_NOQSO	float32	Difference in reduced chi-squared between best and second best fit for Z_NOQSO redshift
# Z_PERSON	float32	Visual-inspection redshift
# CLASS_PERSON	int32	Visual-inspection classification (0=not inspected or unknown, 1=star, 2=narrow emission-line galaxy, 3=QSO, 4=galaxy)
# Z_CONF_PERSON	int32	Visual-inspection confidence (0=not inspected or no confidence, 1,2=low confidence, 3,4=high confidence)
# COMMENTS_PERSON	str	Visual-inspection comments
# CALIBFLUX	float32[5]	ugriz fluxes used for calibrations (nanomaggies)
# CALIBFLUX_IVAR	float32[5]	Inverse variances of ugriz fluxes used for calibrations (nanomaggies)
