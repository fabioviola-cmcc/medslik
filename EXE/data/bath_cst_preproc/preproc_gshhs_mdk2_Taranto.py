#!/usr/bin/python
# pre_processor MEDSLIK - II

# Developed for coastline extraction from gshhs full resolution. 
# It will extract coastline files for whichever part of the globe you
# wish and make it compatible with MEDSLIK-II. 

import numpy as np
from netCDF4 import Dataset
import os

# extract grid corners from ocean fields
def oce_grid(filename):

	fOCM = Dataset(filename)
	x_mod = fOCM.variables['nav_lon'][:]
	y_mod = fOCM.variables['nav_lat'][:]
	return x_mod, y_mod

# Getting header
def getheader(fid):

	A=np.fromfile(fid,dtype='>i',count=8,sep="")

	cnt=A.size

	if cnt<8:  # This gets triggered by the EOF
		return

	else:
		# ver=bitand(bitshift(A(3),-8),255);
		a=A[2] >> 8
		b=255
		ver= a & b

	# Assuming that the GSHHS version used is recent (and comes with MEDSLIK-II Installation pack)

		level = A[2] & b
		a=A[2] >> 16
		greenwich= a & b

		a=A[2] >> 24
		source=a & b

		A[2]=level

		#A[8]=greenwich*65536

		if A.size<=8:
			A=np.hstack((A,greenwich*65536))
		else:
			A[8]=greenwich*65536

		A2=np.fromfile(fid,dtype='>i',count=3,sep="")


		return A, cnt

def extract_coastline(x_mod,y_mod,input_file,output_file):

	# Define limits
	iLatitudeMin= np.min(y_mod)
	iLatitudeMax= np.max(y_mod)
	iLongitudeMin= np.min(x_mod)
	iLongitudeMax= np.max(x_mod)

	lats=np.array([iLatitudeMin, iLatitudeMax])
	longs=np.array([iLongitudeMin, iLongitudeMax])

	llim=np.remainder(longs[0]+360,360)*1e6
	rlim=np.remainder(longs[1]+360,360)*1e6
	tlim=lats[1]*1e6
	blim=lats[0]*1e6


	mrlim=np.remainder(longs[1]+360+180,360)-180
	mllim=np.remainder(longs[0]+360+180,360)-180
	mtlim=lats[1]
	mblim=lats[0]


	# Preparing matrices
	ncst=np.nan+np.zeros((10561669,2))
	Area=np.zeros((188611,1))
	k=np.ones((188612,1))
	decfac=12500

	# Open GSHHS_f file

	fid = open(input_file, 'r')

	Area2=Area.copy()

	A,cnt=getheader(fid)


	l=-1

	while (cnt>0):

	# A: 1:ID, 2:num points, 3:land/lake etc., 4-7:w/e/s/n, 8:area, 9:greenwich crossed.

		iPointsInSegment=A[1]*2
		C=np.fromfile(fid,dtype='>i',count=iPointsInSegment,sep="")  # Read all points in the current segment

	#For Antarctica the lime limits are 0 to 360 (exactly), thus c==0 and the
	#line is not chosen for (e.g. a conic projection of part of Antarctica)
	#Fix 30may/02

		if A[4]==360e6: A[4]=A[4]-1

		a=rlim>llim # Map limits cross longitude jump? (a==1 is no)
		b=A[8]<65536 # Cross boundary? (b==1 if no).
		c=llim<np.remainder(A[4]+360e6,360e6)
		d=rlim>np.remainder(A[3]+360e6,360e6)
		e=tlim>A[5] and blim<A[6]

	#  This test checks whether the lat/long box containing the line overlaps that of
	# the map. There are various cases to consider, depending on whether map limits
	# and/or the line limits cross the longitude jump or not.


		if e and (  (a and ( (b and c and d) or (not b and (c or d)) )) or (not a and (not b or (b and (c or d))) ) ):

			l=l+1
			iLengthC=C.size

			x=C[np.arange(0,iLengthC,2)]*1e-6
			y=C[np.arange(1,iLengthC,2)]*1e-6
			dx=np.diff(x,n=1,axis=0)

			higher=dx>356
			lower=dx<-356
			eastern=x[0]>180
			x=x-360*np.cumsum([np.hstack((eastern.astype(int),higher.astype(int) - lower.astype(int)))])

	#Antarctic is a special case - extend contour to make nice closed polygon
	#that doesn't surround the pole.

			if np.abs(x[0])<1 and np.abs(y[0]+68.9)<1:

				y=np.concatenate((-89.9, -78.4, y[x<=-180], y[x>-180], -78.4, -89.9*np.ones((18,1))))
				iAntarctic=np.arange(-180,161,20)
				x=np.concatenate((180, 180,x[x<=-180]+360,x[x>-180],-180, iAntarctic.conj().transpose()))
				del iAntarctic


	# First and last point should be the same IF THIS IS A POLYGON
	# if the Area=0 then this is a line, and don't add points!

			if A[7]>0:

				if x[-1] != x[0] or y[-1] != y[0]:

					x=np.hstack((x,x[0]))
					y=np.hstack((y,y[0]))


	# get correct curve orientation for patch-fill algorithm.

				iCurveOrientation=np.diff(x,n=1,axis=0) * ( y[0:-1] + y[1:] ) / 2
				Area2[l]=iCurveOrientation.sum(axis=0),
				del iCurveOrientation
				Area[l]=A[7] / 10

				# ok for now

				if np.remainder(A[2],2)==0:
					Area[l]=-np.abs(Area[l]);

					if Area2[l]>0:

						x=x[::-1]
						y=y[::-1]
						#ok for now
				else:
					if Area2[l]<0:
						x=x[::-1]
						y=y[::-1]
						#ok for now

			else:

	# Later on 2 point lines are clipped so we want to avoid that

				if x.size==2:

					x=np.hstack((x[0],x.mean(axis=0), x[1]))
					y=np.hstack((y[0],y.mean(axis=0), y[1]))


	# Here we try to reduce the number of points.

			xflag=0

			if x.max(0)>180: # First, save original curve for later if we anticipate

				sx=x.copy()
				sy=y.copy()   # a 180-problem.
				xflag=1

	#Look for points outside the lat/long boundaries, and then decimate them
	 # by a factor of about 'decfac' (don't get rid of them completely because that
	   # can sometimes cause problems when polygon edges cross curved map edges).

			tol=0.2

			#ok for now

	#Do y limits, then x so we can keep corner points.

			nn=np.logical_or((y>mtlim+tol),(y<mblim-tol))

	# keep one extra point when crossing limits, also the beginning/end point.

			bb=np.concatenate(([0],np.diff(nn.astype(int),n=1,axis=0)))>0

			cc=np.concatenate((np.diff(nn.astype(int),n=1,axis=0),[0]))<0

			dd=bb.astype(int) + cc.astype(int)

			minlon=dd>dd.min(0)

			nn=nn-minlon.astype(int)

			nn[0]=0
			nn[-1]=0

			del bb
			del cc
			del dd
			del minlon

	# decimate vigorously

			bb=np.remainder(np.arange(1,np.size(nn)+1,1),decfac)
			nn=np.logical_and(nn,bb.astype(int))

			x=np.delete(x,nn.nonzero(),axis=0)
			y=np.delete(y,nn.nonzero(),axis=0)

			#ok for now

			if mrlim>mllim: # no wraparound

				aa=np.logical_or(x>mrlim+tol,x<mllim-tol)
				bb=np.logical_and(aa,y<mtlim)
				nn=np.logical_and(bb,y>mblim)

				#ok for now

			else:

				aa=np.logical_and(x>mrlim+tol,x<mllim-tol)
				bb=np.logical_and(aa,y<mtlim)
				nn=np.logical_and(bb,y>mblim)



	# ------------------------------------------------------------------------------------------------------
			del bb

			bb=np.concatenate(([0],np.diff(nn.astype(int),n=1,axis=0)))>0
			cc=np.concatenate((np.diff(nn.astype(int),n=1,axis=0),[0]))<0
			dd=bb.astype(int) + cc.astype(int)

			minlon=dd>dd.min(0)
			nn=nn-minlon.astype(int)

			nn[0]=0
			nn[-1]=0

			del bb
			del cc
			del dd
			del aa

			#ok for now

			# decimating vigorously again

			bb=np.remainder(np.arange(1,np.size(nn)+1,1),decfac)

			nn=np.logical_and(nn.astype(int),bb.astype(int))

			x=np.delete(x,nn.nonzero(),axis=0)
			y=np.delete(y,nn.nonzero(),axis=0)

			# Move all points "near" to map boundaries.
			# I'm not sure about the wisdom of this - it might be better to clip
			# to the boundaries instead of moving. Hmmm.


			y[y>(mtlim+tol)]=mtlim+tol

			y[y<(mblim-tol)]=mblim-tol

			# ok for now

			if mrlim>mllim: # Only clip long bdys if I can tell I'm on the right
			# or left (i.e. not in wraparound case
				x[x>(mrlim+tol)]=mrlim+tol
				x[x<(mllim-tol)]=mllim-tol

			# ok for now

			# plot(x,y);pause;

			k[l+1]=k[l]+x.size+1;

			ncst[np.arange((k[l].astype(int)),(k[l+1].astype(int)-1),1),0]=x

			ncst[np.arange((k[l].astype(int)),(k[l+1].astype(int)-1),1),1]=y

			ncst[k[l+1].astype(int),]=[np.nan,np.nan]

			# This is a little tricky...the filling algorithm expects data to be in the
			# range -180 to 180 deg long. However, there are some land parts that cut across
			# this divide so they appear at +190 but not -170. This causes problems later...
			# so as a kludge I replicate some of the problematic features at 190-360=-170.
			# Small islands are just duplicated, for the Eurasian landmass I just clip
			# off the eastern part.

			if xflag:

				l=l+1
				Area[l]=Area[l-1]

				#ok for now

				if abs(Area[l])>1e5:

					iIndexBase=np.arange(0,sx.size,1)

					nn=iIndexBase[np.where(sx>180)]
					nn=np.hstack((nn,nn[0]))

					k[l+1]=k[l]+nn.size+1

					ncst[np.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),0]=sx[nn]-360
					ncst[np.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),1]=sy[nn]

					# ok for now

				else:   # repeat the island at the other edge.

					k[l+1]=k[l]+sx.size+1

					ncst[np.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),0]=sx-360
					ncst[np.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),1]=sy


				ncst[k[l+1].astype(int),]=[np.nan,np.nan]


		try:
			A,cnt=getheader(fid)
		except:
			print 'Data extraction process has finished. Slicing...'
			break

	#Getting rid of unused rows

	ncst=np.delete(ncst, (np.arange(k[l+1].astype(int), ncst.size/2)), axis=0) # get rid of unused part of data matrices
	Area=np.delete(Area, (np.arange((l+1), Area.size+1)), axis=0)
	k=np.delete( k , ( np.arange ( (l+2) , k.size+1)), axis=0)

	print '...Done. Printing it on a txt file...'

	iSegmentLength=np.diff(k,n=1,axis=0)

	# Organizing segments according to their length -- descending order
	IX=iSegmentLength.argsort(axis=0)[::-1]

	# writes the first line of the .map file. It should contain the # of "isles"
	CoastFile=open(output_file,'w')
	CoastFile.write("%-4.0f \n" % (IX.size))
	iTargetSites=[]
        iTargetSite=np.array([0.,0.,0.,0.])

	# prints coastline data island by island


	for t in IX:

		# extracts the island
		iCoastalSection=ncst[np.arange(k[t].astype(int),k[t+1].astype(int)-1,1),]
		iCoastalLength=iCoastalSection.size/2
		#prints the length of the island
		CoastFile.write("%-4.0f %1.0f \n" % (iCoastalLength,0))
		#prints data related to the island
		for r in np.arange(0,iCoastalSection.size/2,1):
			CoastFile.write("%-8.5f %-6.5f \n" % (iCoastalSection[r,0], iCoastalSection[r,1]))

		for i in range(0,len(iCoastalSection)-1):
                        iTargetSite[0]=iCoastalSection[i,0]
                        iTargetSite[1]=iCoastalSection[i,1]
                        iTargetSite[2]=iCoastalSection[i+1,0]
                        iTargetSite[3]=iCoastalSection[i+1,1]
                        iTargetSites.append(list(iTargetSite))

	CoastFile.close()
	np.savetxt(output_file[:-3] + 'txt', np.array(iTargetSites))


################################################################################
# USER INPUTS
################################################################################
    
# set path for your original gebco netCDF file
gshhs_filename = '/work/tessa15/meglob/EXE_test/data/bath_cst_preproc/gshhs_f.b'

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

# open an ocean forecast file
###oce_filename = (oce_dir + "/" + os.listdir(oce_dir)[2])
oce_filename = (oce_dir + "/MDK_ocean_190804_T.nc")

x_mod, y_mod = oce_grid(oce_filename)

# extract coastline
output_file = output_dir + "/medf.map.taranto"
extract_coastline(x_mod,y_mod,gshhs_filename,output_file)
