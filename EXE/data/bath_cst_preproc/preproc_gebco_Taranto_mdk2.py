#!/usr/bin/python

# Application:
# to extract subsets from GEBCO 30'' files and transform it into MEDSLIK-II
# compatible .bath files

import netCDF4
import numpy as np
from scipy.interpolate import RectBivariateSpline
import os


# Functions
# open GEBCO
def load_gebco(filename):

	fGEBCO = netCDF4.Dataset(filename)
	zz=fGEBCO.variables["z"][:]
	cols, rows = fGEBCO.variables["dimension"]
	iDs = fGEBCO.variables["spacing"][0]

	# Coordinates for GEBCO corners - LAT,LON
	iNECorner=np.array([89.+(59./60)+45./(60*60), -179.-(59./60)-45./(60*60)])
	iSWCorner=np.array([-90.,180.])

	iLatitude=np.arange(iNECorner[0],iSWCorner[0],-iDs)
	iLongitude=np.arange(iNECorner[1],iSWCorner[1],iDs)

	# Reshape bathymetry into m x n matrix
	
	a=np.shape(iLongitude)[0]
	b=np.shape(iLatitude)[0]
	Z = zz.reshape(b, a)
	return iLatitude, iLongitude, Z

# crop GEBCO
def interp_gebco(iLatitude, iLongitude, Z, x_mod, y_mod):

	iLatitudeMin= np.min(y_mod)
	iLatitudeMax= np.max(y_mod)
	iLongitudeMin= np.min(x_mod)
	iLongitudeMax= np.max(x_mod)

	# Crop to area of interest
	iLonIndex=np.argwhere((iLongitude>=iLongitudeMin) & (iLongitude<=iLongitudeMax))
	iLatIndex=np.argwhere((iLatitude>=iLatitudeMin) & (iLatitude<=iLatitudeMax))

	x_crop, y_crop = np.meshgrid(iLongitude[iLonIndex],iLatitude[iLatIndex])
	z_crop = Z[np.min(iLatIndex):(np.max(iLatIndex)+1),np.min(iLonIndex):(np.max(iLonIndex)+1)]
    	#z_crop[z_crop>0] = 0.
    

	# Generate interpolator
	y_crop = np.flipud(y_crop)
	x_crop = np.flipud(x_crop)
	z_crop = np.flipud(z_crop)

	z_int = RectBivariateSpline(y_crop[:,1],x_crop[1,:],z_crop)

	# Interpolate
	z_proc = z_int(y_mod,x_mod)

     # Fix orientation 
	#y_mod = np.flipud(y_mod)
	#x_mod = np.flipud(x_mod)
	z_proc = np.flipud(z_proc)
    
	# Convert bathymetry to MDK-II
	mdk_z=[]
	mdk_x=[]
	mdk_y=[]


	r_,c_ = z_proc.shape
	for i in range(0,r_,1):
		for j in range(0,c_,1):
			mdk_z.append(z_proc[i,j])
			mdk_x.append(x_mod[j])
			mdk_y.append(y_mod[i])

	mdk_z = np.array(mdk_z)
	land_mask = np.where(mdk_z >= 0)
	mdk_z[land_mask]=-9999
	mdk_z = -1.*(mdk_z)
	return mdk_x, mdk_y, mdk_z, c_, r_

# extract grid corners from ocean fields
def oce_grid(filename):
	fOCM = netCDF4.Dataset(filename)
	x_mod = fOCM.variables['nav_lon'][:]
	y_mod = fOCM.variables['nav_lat'][:]
	return x_mod, y_mod

def open_mercator(sFileName):
	# load file and geo coordinates
	fM=netCDF4.Dataset(sFileName,"r")
	y=fM.variables["nav_lat"][:]
	x=fM.variables["nav_lon"][:]
	iT=fM.variables["vomecrtx"][:] # water temperature in deg C
	return x,y,iT

################################################################################
# USER INPUTS
################################################################################

# set path for your original gebco netCDF file
gebco_filename = '/work/tessa15/meglob/EXE_test/data/bath_cst_preproc/GEBCO_2014_1D.nc'

# define where outputs will be placed (MEDSLIK-adapted nc files)
output_dir= "/work/tessa15/meglob/EXE_test/data/bath_cst_preproc/"

# the bathymetry fields will be interpolated to your hydrodynamic grid
# therefore, let us know where you've placed your current files
oce_dir = "/work/tessa15/meglob/DATA/fcst_data/H3k/"

################################################################################
# USER INPUTS - OVER!
################################################################################
# From here onwards, the script should do everything pretty much automatic
# bugs/errors are expected and, in case you do find one,
# feel free to send us comments/corrections.

# load GEBCO file
iLatitude, iLongitude, Z = load_gebco(gebco_filename)


# open an ocean forecast file
###oce_filename = (oce_dir + "/MDK_ocean_190804")
oce_filename = (oce_dir + "/MDK_ocean_190804_T.nc")
x_mod, y_mod = oce_grid(oce_filename)


# extract bathymetry
mdk_x, mdk_y, mdk_z, c_, r_ = interp_gebco(iLatitude, iLongitude, Z, x_mod, y_mod)

# write .bath file
BathFile=open(output_dir + "medf_.bath.mercator.taranto", "w")
BathFile.write("MEDSLIK-II compatible bathymetry file. Degraded resolution based on GEBCO 30''\n")
BathFile.write("%-7.4f %-7.4f %-7.4f %-7.4f \n" % (np.min(mdk_x),np.max(mdk_x),np.min(mdk_y),np.max(mdk_y)))
BathFile.write("%d %d \n" % (c_,r_))
np.savetxt(BathFile,mdk_z,fmt="%04.0f")
BathFile.close()

