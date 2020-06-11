#!/usr/bin/python

import os
import numpy
import math
import datetime
import netCDF4
import scipy.interpolate


def xp_environ(itosra_dir,xp_name):
	xp_files = itosra_dir + '/' + xp_name + '/xp_files'
	met_files = itosra_dir + '/' + xp_name + '/met_files'
	oce_files = itosra_dir + '/' + xp_name + '/oce_files'
	bnc_files = itosra_dir + '/' + xp_name + '/bnc_files'

	os.system('mkdir ' + itosra_dir + '/' + xp_name)
	os.system('mkdir ' + xp_files)
	os.system('mkdir ' + met_files)
	os.system('mkdir ' + oce_files)
	os.system('mkdir ' + bnc_files)

def runlist_gen(rps_file,itosra_dir, xp_name,grid_corners,oil_API,spill_duration,\
    spill_volume, members_id,start_time,end_time,time_interval_btw_spills):

	x1 = grid_corners[0]
	x2 = grid_corners[1]
	y1 = grid_corners[2]
	y2 = grid_corners[3]

	rps_info = numpy.loadtxt(rps_file)
	rp_id = rps_info[:,0]
	rp_x = rps_info[:,1]
	rp_y = rps_info[:,2]
	xx = numpy.logical_and(rp_x>x1, rp_x<=x2)
	yy = numpy.logical_and(rp_y>y1,rp_y<=y2)
	tt = numpy.logical_and(xx, yy)
	rp_id = rp_id[tt]
	rp_x = rp_x[tt]
	rp_y = rp_y[tt]

	runlist_file = open(itosra_dir + '/' + xp_name + '/xp_files/runlist.txt','w')
	numpy.savetxt(itosra_dir + '/' + xp_name + '/xp_files/rp_ids.txt',rps_info[tt,:])
	numpy.savetxt(itosra_dir + '/' + xp_name + '/xp_files/glamor_rp_ids.txt',rps_info[tt,:],fmt='%6.2f', delimiter=',', newline='\n', header='id, lon, lat',)

	model = 2 # for now only mercator

	# Generating runlist.dat
	for rp in range(0,len(rp_id)):
    		for spill_time_n in numpy.arange(start_time,end_time,time_interval_btw_spills):
        		member_id=1
        		for member in members_id:
            			# MDK setup info
            			oil_type=oil_API[member]
            			duration=spill_duration[member]
            			volume=spill_volume[member]
            			rate=volume/duration
            			spill_time=datetime.date.fromordinal(spill_time_n)

            			# Coast check
            			position_x = math.modf(rp_x[rp])
            			position_y = math.modf(rp_y[rp])

            			simulation_name_string=('%05d' % (rp_id[rp]) + '_' + spill_time.strftime('%d') + spill_time.strftime('%m') + spill_time.strftime('%y')  + '_' + '%1d' % model + '_' + '%02d' % member_id)
            			spill_time_string= (spill_time.strftime('%d') + ' ' + spill_time.strftime('%m') + ' ' + spill_time.strftime('%y') + ' ' + spill_time.strftime('%H') + ' ' +spill_time.strftime('%M'))
            			spill_position_string=('%02d' % position_y[1] +  ' ' + '%02d' % ((position_y[0])*60) +  ' ' + '%02d' % position_x[1] +  ' ' + '%02d' % (position_x[0]*60))
            			characteristics_string=('%04d' % duration + ' ' + '%08.2f' % rate + ' ' +  '%02d' % oil_type)
            			runlist_file.write("%s %s %s %s\n" % (simulation_name_string, spill_time_string, spill_position_string, characteristics_string))
        			member_id=member_id+1
	runlist_file.close()

def open_mercator_subset(sFileName,grid_corners):
	fM=netCDF4.Dataset(sFileName,"r")
	iLatitude=fM.variables["latitude"][:]
	iLongitude=fM.variables["longitude"][:]

	# cropping before extraction
	xmin = grid_corners[0]
	xmax = grid_corners[1]
	ymin = grid_corners[2]
	ymax = grid_corners[3]

	xwindow = numpy.logical_and(iLongitude>=xmin, iLongitude<=xmax)
	ywindow = numpy.logical_and(iLatitude>=ymin, iLatitude<=ymax)

	x = iLongitude[xwindow]
	y = iLatitude[ywindow]


	# Extraction and conversion of variables
	iU=fM.variables["uo"][:,:,ywindow,xwindow]#*scale_velocities + offset_velocities # u component of water velocity in m/s
	iV=fM.variables["vo"][:,:,ywindow,xwindow]#*scale_velocities + offset_velocities  # v component
	iT=fM.variables["thetao"][:,:,ywindow,xwindow]#*scale_temperature + offset_temperature # water temperature in deg C
	iTime=fM.variables["time"][:]  # hours since 1950-01-01 00:00:00
	iTime=datetime.date(1950,01,1).toordinal() + iTime/24 # time converted to standard pythonic time
	iDepth=fM.variables["depth"][:]          # in m
	fM.close()
	return x,y,iU,iV,iT,iTime,iDepth

def interp_z(iU,iDepth):

	iShapeO=numpy.shape(iU)
	# daily resolution and four depths - 0,10,30,120m
	iOutputUd=numpy.zeros((iShapeO[0],4,iShapeO[2],iShapeO[3]))	# u(time,depth,y,x)
	# surface
	iOutputUd[:,0,:,:]=iU[:,0,:,:]
	# 10 m
	iOutputUd[:,1,:,:]=iU[:,7,:,:] + (iU[:,8,:,:]-iU[:,7,:,:])*((10. - iDepth[7])/(iDepth[8] - iDepth[7]))
	# 30 m
	iOutputUd[:,2,:,:]=iU[:,14,:,:] + (iU[:,15,:,:]-iU[:,14,:,:])*((30. - iDepth[14])/(iDepth[15] - iDepth[14]))
	# 120 m
	iOutputUd[:,3,:,:]=iU[:,22,:,:] + (iU[:,23,:,:]-iU[:,22,:,:])*((120. - iDepth[22])/(iDepth[23] - iDepth[22]))
	return numpy.squeeze(iOutputUd)

def interp_t(iOutputUd, iHRTimeLine, iLRTimeLine):
		# prepare storing matrixes
		iShapeO = numpy.shape(iOutputUd)
		iOutputU=numpy.zeros((24,4,iShapeO[2],iShapeO[3]))-32767.0	# u(time,depth,y,x)
		#interpolating
		for ii in range(0,iShapeO[2]):
			for jj in range(0,iShapeO[3]):
				for dd in range(0,4):
					if iOutputUd[0,dd,ii,jj]>-20:
						iOutputU[:,dd,ii,jj]=numpy.interp(iHRTimeLine,iLRTimeLine,iOutputUd[:,dd,ii,jj])
						#Extrapolating data
						if (dd!=0 and any(iOutputU[:,dd,ii,jj]<-2) and all(iOutputU[:,dd-1,ii,jj]>-2)):
							iOutputU[:,dd,ii,jj]=iOutputU[:,dd-1,ii,jj]
		return iOutputU

def curr_gen(rp_info, start_date,end_date,rawdata_folder,\
	itosra_dir, xp_name, expected_displacement, ds, input_ref_hours, output_ref_hours):

	rp_id = rp_info[0]
	rp_x = rp_info[1]
	rp_y = rp_info[2]
	output_folder = itosra_dir + '/' + xp_name + '/oce_files'

	# create output dirs
	os.system('mkdir ' + output_folder + '/H3k_' + "%05i" % (rp_id,))
	pt_folder= (output_folder + "/H3k_" + "%05i" % (rp_id,) + "/")

	# set the limits for the extraction of the current fields
	xmin = rp_x- expected_displacement*ds
	xmax = rp_x+ expected_displacement*ds
	ymin = rp_y- expected_displacement*ds
	ymax = rp_y+ expected_displacement*ds
	grid_corners = [xmin,xmax,ymin,ymax]


	counter = 0
	for i in range(start_date,end_date,1):

		iDateY=datetime.date.fromordinal(int(i-1))
		iDate=datetime.date.fromordinal(int(i))
		iDateT=datetime.date.fromordinal(int(i+1))

		if os.path.isfile(pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_U.nc'):
			print('... file found. Next time step...')
			continue

		# Open MERCATOR files and interpolate in z for preset dephts 0,10,30,120m
		if counter == 0:
			iLongitude,iLatitude,iU_Y,iV_Y,iT_Y,iTime_Y,iDepth = open_mercator_subset(rawdata_folder + 'mercator_paria_20' + iDateY.strftime('%y') + iDateY.strftime('%m') + iDateY.strftime('%d') + '.nc',grid_corners)
			iU_Y = (interp_z(iU_Y,iDepth))
			iV_Y = (interp_z(iV_Y,iDepth))
			iT_Y = (interp_z(iT_Y,iDepth))
			# day of interest
			iLongitude,iLatitude,iU_P,iV_P,iT_P,iTime,iDepth = open_mercator_subset(rawdata_folder + 'mercator_paria_20' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '.nc',grid_corners)
			iU_P = (interp_z(iU_P,iDepth))
			iV_P = (interp_z(iV_P,iDepth))
			iT_P = (interp_z(iT_P,iDepth))
			#  day after
			iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = open_mercator_subset(rawdata_folder + 'mercator_paria_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '.nc',grid_corners)
			iU_T = (interp_z(iU_T,iDepth))
			iV_T = (interp_z(iV_T,iDepth))
			iT_T = (interp_z(iT_T,iDepth))

			counter = counter + 1
		else:
			# Present becomes yesterday -- oh life...
			iU_Y = iU_P
			iV_Y = iV_P
			iT_Y = iT_P

			# Tomorrow becomes today -- oh dear life... time flies.
			iU_P = iU_T
			iV_P = iV_T
			iT_P = iT_T

			iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = open_mercator_subset(rawdata_folder + 'mercator_paria_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '.nc',grid_corners)
			iU_T = (interp_z(iU_T,iDepth))
			iV_T = (interp_z(iV_T,iDepth))
			iT_T = (interp_z(iT_T,iDepth))

		# Put the three fields (t - 1, t, t + 1) into a single matrix to interpolate
		iOutputUd = numpy.stack((iU_Y,iU_P,iU_T),axis=0)
		iOutputVd = numpy.stack((iV_Y,iV_P,iV_T),axis=0)
		iOutputTd = numpy.stack((iT_Y,iT_P,iT_T),axis=0)

		# Temporal interpolation:
		iOutputU = interp_t(iOutputUd, output_ref_hours, input_ref_hours)
		iOutputV = interp_t(iOutputVd, output_ref_hours, input_ref_hours)
		iOutputT = interp_t(iOutputTd, output_ref_hours, input_ref_hours)

		# Masking
		iOutputU=numpy.where(iOutputU < -3,9999,iOutputU)
		iOutputV=numpy.where(iOutputV < -3,9999,iOutputV)
		iOutputT=numpy.where(iOutputT < -4,9999,iOutputT)

		iOutputU=numpy.where(iOutputU > 3,9999,iOutputU)
		iOutputV=numpy.where(iOutputV > 3,9999,iOutputV)
		iOutputT=numpy.where(iOutputT > 40,9999,iOutputT)

		#Generates NC files
		# U component
		sConvertedFilename=pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_U.nc'
		f = netCDF4.Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
		f.history = 'MERCATOR Forecast - MyOcean -- Adapted to Medslik-II'
		f.createDimension('time_counter',None)
		f.createDimension('y', len(iLatitude))
		f.createDimension('x', len(iLongitude))
		f.createDimension('depthu',4)
		time = f.createVariable('time_counter', 'd', ('time_counter',))
		# Geo coordinates
		yCoord= f.createVariable('nav_lat', 'd', ('y'))
		xCoord= f.createVariable('nav_lon', 'd', ('x'))
		# Zonal velocity component
		uCoord= f.createVariable('vozocrtx', 'd',('time_counter','depthu','y','x'))
		time[:] = i + numpy.arange(1./24,25./24,1./24)
		#Filling the sausage
		yCoord[:]=iLatitude
		xCoord[:]=iLongitude
		#Filling the sausage chap. II
		uCoord[:]=iOutputU
		#Adding units
		yCoord.units = 'degrees_north'
		xCoord.units = 'degrees_east'
		time.units = 'julian days'
		uCoord.units='m*s-1'
		uCoord.missing_value=9999
		f.close()

		# V component
		sConvertedFilename=pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_V.nc'
		f = netCDF4.Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
		f.history = 'MERCATOR Forecast - MyOcean -- Adapted to Medslik-II'
		f.createDimension('time_counter',None)
		f.createDimension('y', len(iLatitude))
		f.createDimension('x', len(iLongitude))
		f.createDimension('depthv',4)
		time = f.createVariable('time_counter', 'd', ('time_counter',))
		# Geo coordinates
		yCoord= f.createVariable('nav_lat', 'd', ('y'))
		xCoord= f.createVariable('nav_lon', 'd', ('x'))
		# Zonal velocity component
		vCoord= f.createVariable('vomecrty', 'd',('time_counter','depthv','y','x'))
		time[:] = i + numpy.arange(1./24,25./24,1./24)
		#Filling the sausage
		yCoord[:]=iLatitude
		xCoord[:]=iLongitude
		vCoord[:]=iOutputV
		#Adding units
		yCoord.units = 'degrees_north'
		xCoord.units = 'degrees_east'
		time.units = 'julian days'
		vCoord.units='m*s-1'
		vCoord.missing_value=9999
		f.close()

		# Temperature
		sConvertedFilename=pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_T.nc'
		f = netCDF4.Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
		f.history = 'MERCATOR Forecast - MyOcean -- Adapted to Medslik-II v2.0'
		f.createDimension('time_counter',None)
		f.createDimension('y', len(iLatitude))
		f.createDimension('x', len(iLongitude))
		f.createDimension('deptht',4)
		time = f.createVariable('time_counter', 'd', ('time_counter',))
		# Geo coordinates
		yCoord= f.createVariable('nav_lat', 'd', ('y'))
		xCoord= f.createVariable('nav_lon', 'd', ('x'))
		# Zonal velocity component
		tCoord= f.createVariable('votemper', 'd',('time_counter','deptht','y','x'))
		time[:] = i + numpy.arange(1./24,25./24,1./24)
		#Filling the sausage
		yCoord[:]=iLatitude
		xCoord[:]=iLongitude
		tCoord[:]=iOutputT
		#Adding units
		yCoord.units = 'degrees_north'
		xCoord.units = 'degrees_east'
		time.units = 'julian days'
		tCoord.units='degree Celsius'
		tCoord.missing_value=9999
		f.close()

# extract grid corners from ocean fields
def oce_grid(filename):

	fOCM = netCDF4.Dataset(filename)
	x_mod = fOCM.variables['nav_lon'][:]
	y_mod = fOCM.variables['nav_lat'][:]
	fOCM.close()
	return x_mod, y_mod

# Getting header
def getheader(fid):

	A=numpy.fromfile(fid,dtype='>i',count=8,sep="")

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
			A=numpy.hstack((A,greenwich*65536))
		else:
			A[8]=greenwich*65536

		A2=numpy.fromfile(fid,dtype='>i',count=3,sep="")


		return A, cnt

def extract_coastline(x_mod,y_mod,sOutputFile):

	# Define limits
	iLatitudeMin= numpy.min(y_mod)
	iLatitudeMax= numpy.max(y_mod)
	iLongitudeMin= numpy.min(x_mod)
	iLongitudeMax= numpy.max(x_mod)

	# Path for your full resolution GSHHS file

	file='/Users/asepp/work/gshhs_i.b'

	lats=numpy.array([iLatitudeMin, iLatitudeMax])
	longs=numpy.array([iLongitudeMin, iLongitudeMax])

	llim=numpy.remainder(longs[0]+360,360)*1e6
	rlim=numpy.remainder(longs[1]+360,360)*1e6
	tlim=lats[1]*1e6
	blim=lats[0]*1e6


	mrlim=numpy.remainder(longs[1]+360+180,360)-180
	mllim=numpy.remainder(longs[0]+360+180,360)-180
	mtlim=lats[1]
	mblim=lats[0]


	# Preparing matrices
	ncst=numpy.nan+numpy.zeros((10561669,2))
	Area=numpy.zeros((188611,1))
	k=numpy.ones((188612,1))
	decfac=12500

	# Open GSHHS_f file

	fid = open(file, 'r')

	Area2=Area.copy()

	A,cnt=getheader(fid)


	l=-1

	while (cnt>0):

	# A: 1:ID, 2:num points, 3:land/lake etc., 4-7:w/e/s/n, 8:area, 9:greenwich crossed.

		iPointsInSegment=A[1]*2
		C=numpy.fromfile(fid,dtype='>i',count=iPointsInSegment,sep="")  # Read all points in the current segment

	#For Antarctica the lime limits are 0 to 360 (exactly), thus c==0 and the
	#line is not chosen for (e.g. a conic projection of part of Antarctica)
	#Fix 30may/02

		if A[4]==360e6: A[4]=A[4]-1

		a=rlim>llim # Map limits cross longitude jump? (a==1 is no)
		b=A[8]<65536 # Cross boundary? (b==1 if no).
		c=llim<numpy.remainder(A[4]+360e6,360e6)
		d=rlim>numpy.remainder(A[3]+360e6,360e6)
		e=tlim>A[5] and blim<A[6]

	#  This test checks whether the lat/long box containing the line overlaps that of
	# the map. There are various cases to consider, depending on whether map limits
	# and/or the line limits cross the longitude jump or not.


		if e and (  (a and ( (b and c and d) or (not b and (c or d)) )) or (not a and (not b or (b and (c or d))) ) ):

			l=l+1
			iLengthC=C.size

			x=C[numpy.arange(0,iLengthC,2)]*1e-6
			y=C[numpy.arange(1,iLengthC,2)]*1e-6
			dx=numpy.diff(x,n=1,axis=0)

			higher=dx>356
			lower=dx<-356
			eastern=x[0]>180
			x=x-360*numpy.cumsum([numpy.hstack((eastern.astype(int),higher.astype(int) - lower.astype(int)))])

	#Antarctic is a special case - extend contour to make nice closed polygon
	#that doesn't surround the pole.

			if numpy.abs(x[0])<1 and numpy.abs(y[0]+68.9)<1:

				y=numpy.concatenate((-89.9, -78.4, y[x<=-180], y[x>-180], -78.4, -89.9*numpy.ones((18,1))))
				iAntarctic=numpy.arange(-180,161,20)
				x=numpy.concatenate((180, 180,x[x<=-180]+360,x[x>-180],-180, iAntarctic.conj().transpose()))
				del iAntarctic


	# First and last point should be the same IF THIS IS A POLYGON
	# if the Area=0 then this is a line, and don't add points!

			if A[7]>0:

				if x[-1] != x[0] or y[-1] != y[0]:

					x=numpy.hstack((x,x[0]))
					y=numpy.hstack((y,y[0]))


	# get correct curve orientation for patch-fill algorithm.

				iCurveOrientation=numpy.diff(x,n=1,axis=0) * ( y[0:-1] + y[1:] ) / 2
				Area2[l]=iCurveOrientation.sum(axis=0),
				del iCurveOrientation
				Area[l]=A[7] / 10

				# ok for now

				if numpy.remainder(A[2],2)==0:
					Area[l]=-numpy.abs(Area[l]);

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

					x=numpy.hstack((x[0],x.mean(axis=0), x[1]))
					y=numpy.hstack((y[0],y.mean(axis=0), y[1]))


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

			nn=numpy.logical_or((y>mtlim+tol),(y<mblim-tol))

	# keep one extra point when crossing limits, also the beginning/end point.

			bb=numpy.concatenate(([0],numpy.diff(nn.astype(int),n=1,axis=0)))>0

			cc=numpy.concatenate((numpy.diff(nn.astype(int),n=1,axis=0),[0]))<0

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

			bb=numpy.remainder(numpy.arange(1,numpy.size(nn)+1,1),decfac)
			nn=numpy.logical_and(nn,bb.astype(int))

			x=numpy.delete(x,nn.nonzero(),axis=0)
			y=numpy.delete(y,nn.nonzero(),axis=0)

			#ok for now

			if mrlim>mllim: # no wraparound

				aa=numpy.logical_or(x>mrlim+tol,x<mllim-tol)
				bb=numpy.logical_and(aa,y<mtlim)
				nn=numpy.logical_and(bb,y>mblim)

				#ok for now

			else:

				aa=numpy.logical_and(x>mrlim+tol,x<mllim-tol)
				bb=numpy.logical_and(aa,y<mtlim)
				nn=numpy.logical_and(bb,y>mblim)



	# ------------------------------------------------------------------------------------------------------
			del bb

			bb=numpy.concatenate(([0],numpy.diff(nn.astype(int),n=1,axis=0)))>0
			cc=numpy.concatenate((numpy.diff(nn.astype(int),n=1,axis=0),[0]))<0
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

			bb=numpy.remainder(numpy.arange(1,numpy.size(nn)+1,1),decfac)

			nn=numpy.logical_and(nn.astype(int),bb.astype(int))

			x=numpy.delete(x,nn.nonzero(),axis=0)
			y=numpy.delete(y,nn.nonzero(),axis=0)

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

			ncst[numpy.arange((k[l].astype(int)),(k[l+1].astype(int)-1),1),0]=x

			ncst[numpy.arange((k[l].astype(int)),(k[l+1].astype(int)-1),1),1]=y

			ncst[k[l+1].astype(int),]=[numpy.nan,numpy.nan]

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

					iIndexBase=numpy.arange(0,sx.size,1)

					nn=iIndexBase[numpy.where(sx>180)]
					nn=numpy.hstack((nn,nn[0]))

					k[l+1]=k[l]+nn.size+1

					ncst[numpy.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),0]=sx[nn]-360
					ncst[numpy.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),1]=sy[nn]

					# ok for now

				else:   # repeat the island at the other edge.

					k[l+1]=k[l]+sx.size+1

					ncst[numpy.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),0]=sx-360
					ncst[numpy.arange((k[l].astype(int)),(k[l+1].astype(int))-1,1),1]=sy


				ncst[k[l+1].astype(int),]=[numpy.nan,numpy.nan]


		try:
			A,cnt=getheader(fid)
		except:
			print 'Data extraction process has finished. Slicing...'
			break

	#Getting rid of unused rows

	ncst=numpy.delete(ncst, (numpy.arange(k[l+1].astype(int), ncst.size/2)), axis=0) # get rid of unused part of data matrices
	Area=numpy.delete(Area, (numpy.arange((l+1), Area.size+1)), axis=0)
	k=numpy.delete( k , ( numpy.arange ( (l+2) , k.size+1)), axis=0)

	print '...Done. Printing it on a txt file...'

	iSegmentLength=numpy.diff(k,n=1,axis=0)

	# Organizing segments according to their length -- descending order
	IX=iSegmentLength.argsort(axis=0)[::-1]

	# writes the first line of the .map file. It should contain the # of "isles"
	CoastFile=open(sOutputFile,'w')
	CoastFile.write("%-4.0f \n" % (IX.size))
	iTargetSites=[]
        iTargetSite=numpy.array([0.,0.,0.,0.])

	# prints coastline data island by island


	for t in IX:

		# extracts the island
		iCoastalSection=ncst[numpy.arange(k[t].astype(int),k[t+1].astype(int)-1,1),]
		iCoastalLength=iCoastalSection.size/2
		#prints the length of the island
		CoastFile.write("%-4.0f %1.0f \n" % (iCoastalLength,0))
		#prints data related to the island
		for r in numpy.arange(0,iCoastalSection.size/2,1):
			CoastFile.write("%-8.5f %-6.5f \n" % (iCoastalSection[r,0], iCoastalSection[r,1]))

		for i in range(0,len(iCoastalSection)-1):
                        iTargetSite[0]=iCoastalSection[i,0]
                        iTargetSite[1]=iCoastalSection[i,1]
                        iTargetSite[2]=iCoastalSection[i+1,0]
                        iTargetSite[3]=iCoastalSection[i+1,1]
                        iTargetSites.append(list(iTargetSite))

	CoastFile.close()
	numpy.savetxt(sOutputFile[:-3] + 'txt', numpy.array(iTargetSites))

def coast_gen(rp_id,itosra_dir,xp_name):

	oce_files = itosra_dir + '/' + xp_name + '/oce_files'
	bnc_files = itosra_dir + '/' + xp_name + '/bnc_files'

	# get the geo limits for your area
	oce_filename = os.popen('ls ' + oce_files + "/H3k_" + "%05i" % (rp_id[0],) + ' | sort -n | head -1').read()
	x_mod, y_mod = oce_grid(oce_files + "/H3k_" + "%05i" % (rp_id[0],) + "/" + oce_filename[:-1])

	# create output dir
	os.system("mkdir " + bnc_files + "/BnC_" + "%05i" % (rp_id[0],))

	# extract coastline
	sOutputFile = bnc_files + "/BnC_" + "%05i" % (rp_id[0],) + "/medf.map"
	extract_coastline(x_mod,y_mod,sOutputFile)

def load_gebco(filename):

	fGEBCO = netCDF4.Dataset(filename)
	zz=fGEBCO.variables["z"][:]
	cols, rows = fGEBCO.variables["dimension"]
	iDs = fGEBCO.variables["spacing"][0]

	# Coordinates for GEBCO corners - LAT,LON
	iNECorner=numpy.array([89.+(59./60)+45./(60*60), -179.-(59./60)-45./(60*60)])
	iSWCorner=numpy.array([-90.,180.])

	iLatitude=numpy.arange(iNECorner[0],iSWCorner[0],-iDs)
	iLongitude=numpy.arange(iNECorner[1],iSWCorner[1],iDs)

	# Reshape bathymetry into m x n matrix
	a=numpy.shape(iLongitude)[0]
	b=numpy.shape(iLatitude)[0]
	Z = zz.reshape(b, a)
	fGEBCO.close()
	return iLatitude, iLongitude, Z

# crop GEBCO
def interp_gebco(iLatitude, iLongitude, Z, x_mod, y_mod):

	iLatitudeMin= numpy.min(y_mod)
	iLatitudeMax= numpy.max(y_mod)
	iLongitudeMin= numpy.min(x_mod)
	iLongitudeMax= numpy.max(x_mod)

	print 'Cropping bat file for:'
	print str(iLongitudeMin), str(iLongitudeMax), ' Longitude'
	print str(iLatitudeMin), str(iLatitudeMax), ' Latitude'
	# Crop to area of interest
	iLonIndex=numpy.argwhere((iLongitude>=iLongitudeMin) & (iLongitude<=iLongitudeMax))
	iLatIndex=numpy.argwhere((iLatitude>=iLatitudeMin) & (iLatitude<=iLatitudeMax))

	x_crop, y_crop = numpy.meshgrid(iLongitude[iLonIndex],iLatitude[iLatIndex])
	z_crop = Z[numpy.min(iLatIndex):(numpy.max(iLatIndex)+1),numpy.min(iLonIndex):(numpy.max(iLonIndex)+1)]

	# Generate interpolator
	y_crop = numpy.flipud(y_crop)
	x_crop = numpy.flipud(x_crop)
	z_crop = numpy.flipud(z_crop)

	z_int = scipy.interpolate.RectBivariateSpline(y_crop[:,1],x_crop[1,:],z_crop)

	# Interpolate
	#z_proc = z_int(y_mod[:,1],x_mod[1,:])
	z_proc = z_int(y_mod,x_mod)

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

	mdk_z = numpy.array(mdk_z)
	land_mask = numpy.where(mdk_z > 0)
	mdk_z[land_mask]=-9999
	mdk_z = -1.*(mdk_z)
	return mdk_x, mdk_y, mdk_z, c_, r_

def bath_gen(rp_id, itosra_dir, xp_name, iLatitude, iLongitude, Z):

	oce_files = itosra_dir + '/' + xp_name + '/oce_files'
	bnc_files = itosra_dir + '/' + xp_name + '/bnc_files'

	# get the geo limits for your area
	oce_filename = os.popen('ls ' + oce_files + "/H3k_" + "%05i" % (rp_id[0],) + ' | sort -n | head -1').read()
	x_mod, y_mod = oce_grid(oce_files + "/H3k_" + "%05i" % (rp_id[0],) + "/" + oce_filename[:-1])
	
	print "Using as reference: "
        print (oce_files + "/H3k_" + "%05i" % (rp_id[0],) + "/" + oce_filename[:-1])
	
	# create output dir
	os.system("mkdir " + bnc_files + "/BnC_" + "%05i" % (rp_id[0],))

	# extract bathymetry
	mdk_x, mdk_y, mdk_z, c_, r_ = interp_gebco(iLatitude, iLongitude, Z, x_mod, y_mod)

	# Write .bath file
	BathFile=open(bnc_files + "/BnC_" + "%05i" % (rp_id[0],) + "/medf_.bath", "w")
	BathFile.write("Bathymetry for global relocatable model. Degraded resolution based on GEBCO 30''\n")
	BathFile.write("%-7.4f %-7.4f %-7.4f %-7.4f \n" % (numpy.min(mdk_x),numpy.max(mdk_x),numpy.min(mdk_y),numpy.max(mdk_y)))
	BathFile.write("%d %d \n" % (c_,r_))
	numpy.savetxt(BathFile,mdk_z,fmt="%04.0f")
	BathFile.close()


def open_era(sFileName):
	fERA=netCDF4.Dataset(sFileName,"r",format="NETCDF4")
	y=fERA.variables["latitude"][:]
	x=fERA.variables["longitude"][:]#-360.

	# Extraction and conversion of variables
	iU=fERA.variables["u10"][:]#*iScale_factor_U + iOffset_U # u component of water velocity in m/s
	iV=fERA.variables["v10"][:]#*iScale_factor_V + iOffset_V # v component
	iTime=fERA.variables["time"][:]/24.            # hours since 1900-01-01 00:00:0.0
	iTime = datetime.date(1900,01,1).toordinal() + iTime
	fERA.close()
	return x,y,iU,iV,iTime

def subset_era(iLongitude,iLatitude,iU,iV,grid_corners):
        # cropping before extraction
        xmin = grid_corners[0]
        xmax = grid_corners[1]
        ymin = grid_corners[2]
        ymax = grid_corners[3]

        xwindow = numpy.where((iLongitude>=xmin) & (iLongitude<=xmax))
        ywindow = numpy.where((iLatitude>=ymin) & (iLatitude<=ymax))

        xs = iLongitude[xwindow]
        ys = iLatitude[ywindow]

        iUs = iU[:,numpy.transpose(ywindow),xwindow]
        iVs = iV[:,numpy.transpose(ywindow),xwindow]

        return xs,ys,iUs,iVs

def open_era_subset(sFileName,grid_corners):
        fM=netCDF4.Dataset(sFileName,"r",format="NETCDF4")
        iLatitude=fM.variables["latitude"][:]
        iLongitude=fM.variables["longitude"][:]

        # cropping before extraction
        xmin = grid_corners[0]
        xmax = grid_corners[1]
        ymin = grid_corners[2]
        ymax = grid_corners[3]

        xwindow = numpy.logical_and(iLongitude>=xmin, iLongitude<=xmax)
        ywindow = numpy.logical_and(iLatitude>=ymin, iLatitude<=ymax)

        x = iLongitude[xwindow]
        y = iLatitude[ywindow]

        # Extraction and conversion of variables
        iU=fM.variables["u10"][:,ywindow,xwindow]
        iV=fM.variables["v10"][:,ywindow,xwindow]
        iTime=fM.variables["time"][:]  # hours since 1900-01-01 00:00:00
        iTime=datetime.date(1900,01,1).toordinal() + iTime/24. # time converted to standard pythonic time
        fM.close()
	return x,y,iU,iV,iTime

def interp_t_wind(iOutputUd, iOutputRefHours, iInputRefHours):
	# prepare storing matrixes
	iShapeO = numpy.shape(iOutputUd)
	iOutputU = numpy.zeros((len(iOutputRefHours),iShapeO[1],iShapeO[2]))

	#interpolating
	for ii in range(0,iShapeO[1]):
		for jj in range(0,iShapeO[2]):
				iOutputU[:,ii,jj]=numpy.interp(iOutputRefHours,iInputRefHours,iOutputUd[:,ii,jj])

	return iOutputU

def wind_gen(rp_info, start_date, end_date, rawdata_folder, itosra_dir, \
    xp_name, expected_displacement, ds, input_ref_hours, output_ref_hours):

	rp_id = rp_info[0]
	rp_x = rp_info[1]
	rp_y = rp_info[2]
	output_folder = itosra_dir + '/' + xp_name + '/met_files'

	# create output dirs
	os.system('mkdir ' + output_folder + '/SK1_' + "%05i" % (rp_id,))
	pt_folder= (output_folder + "/SK1_" + "%05i" % (rp_id,) + "/")

	# set the limits for the extraction of the current fields
	xmin = rp_x- expected_displacement*ds
	xmax = rp_x+ expected_displacement*ds
	ymin = rp_y- expected_displacement*ds
	ymax = rp_y+ expected_displacement*ds
	grid_corners = [xmin,xmax,ymin,ymax]

	for i in range(int(start_date),int(end_date)):

		iDate1=datetime.date.fromordinal(i)
        	iDate2=datetime.date.fromordinal(i+1)

		if os.path.isfile(pt_folder + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc'):
			print('... file found. Next time step...')
			continue
		
		# open ERA files
      		x1,y1,u1,v1,t1 = open_era_subset(rawdata_folder + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc',grid_corners)
	        x2,y2,u2,v2,t2 = open_era_subset(rawdata_folder + iDate2.strftime('%Y') + iDate2.strftime('%m') + iDate2.strftime('%d') + '.nc',grid_corners)
 
		# Temporal interpolation
		iOutputUd = numpy.concatenate((u1,u2),axis=0)[0:5,:,:]
        	iOutputVd = numpy.concatenate((v1,v2),axis=0)[0:5,:,:]
		
		# ERA interim has a 6h temporal resolution
		iOutputU = interp_t_wind(iOutputUd, output_ref_hours, input_ref_hours)
		iOutputV = interp_t_wind(iOutputVd, output_ref_hours, input_ref_hours)

		# Masking
		iOutputU=numpy.where(iOutputU < -100,9999,iOutputU)
		iOutputV=numpy.where(iOutputV < -100,9999,iOutputV)

		#Generates NC file
		sConvertedFilename= pt_folder + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc'
	        f = netCDF4.Dataset(sConvertedFilename,"w","NETCDF3_CLASSIC")
       		f.history = 'ERA Interim 10m winds -- Adapted to Medslik-II'
	        f.createDimension('time',None)
  		f.createDimension('lat', len(y1))
        	f.createDimension('lon', len(x1))
        	time = f.createVariable('time', 'd', ('time',))
        	# Geo coordinates
        	yCoord= f.createVariable('lat', 'd', ('lat',))
        	xCoord= f.createVariable('lon', 'd', ('lon',))
        	# Zonal velocity component
        	uCoord= f.createVariable('U10M', 'd',('time','lat','lon'))
        	vCoord= f.createVariable('V10M', 'd',('time','lat','lon'))
        	time[:] = i + numpy.arange(1./24,24./24,1./24)
        	#Filling the sausage
        	yCoord[:]=y1
        	xCoord[:]=x1
        	#Filling the sausage chap. II
        	uCoord[:]=iOutputU
        	vCoord[:]=iOutputV
        	#Adding units
        	yCoord.units = 'degrees_north'
        	xCoord.units = 'degrees_east'
        	xCoord.lon_min = numpy.min(x1)
        	yCoord.lat_min = numpy.min(y1)
        	time.units = 'pythonic days'
        	uCoord.units='m*s-1'
        	uCoord.missing_value=9999
        	vCoord.units='m*s-1'
        	vCoord.missing_value=9999
        	f.close()


def wind_gen_parallel(f1,f2,rp_info, start_date, end_date, rawdata_folder, itosra_dir, \
    xp_name, expected_displacement, ds, input_ref_hours, output_ref_hours):

        rp_id = rp_info[0]
        rp_x = rp_info[1]
        rp_y = rp_info[2]
        output_folder = itosra_dir + '/' + xp_name + '/met_files'

        # create output dirs
        os.system('mkdir ' + output_folder + '/SK1_' + "%05i" % (rp_id,))
        pt_folder= (output_folder + "/SK1_" + "%05i" % (rp_id,) + "/")

        # set the limits for the extraction of the current fields
        xmin = rp_x- expected_displacement*ds
        xmax = rp_x+ expected_displacement*ds
        ymin = rp_y- expected_displacement*ds
        ymax = rp_y+ expected_displacement*ds
        grid_corners = [xmin,xmax,ymin,ymax]

        for i in range(int(start_date),int(end_date)):

                iDate1=datetime.date.fromordinal(i)
                iDate2=datetime.date.fromordinal(i+1)

                if os.path.isfile(pt_folder + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc'):
                        print('... wind file found. Next time step...')
                        continue

                # open ERA files
                x1,y1,u1,v1,t1 = f1(rawdata_folder + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc',grid_corners)
                x2,y2,u2,v2,t2 = f1(rawdata_folder + iDate2.strftime('%Y') + iDate2.strftime('%m') + iDate2.strftime('%d') + '.nc',grid_corners)

                # Temporal interpolation
                iOutputUd = numpy.concatenate((u1,u2),axis=0)[0:5,:,:]
                iOutputVd = numpy.concatenate((v1,v2),axis=0)[0:5,:,:]

                # ERA interim has a 6h temporal resolution
                iOutputU = f2(iOutputUd, output_ref_hours, input_ref_hours)
                iOutputV = f2(iOutputVd, output_ref_hours, input_ref_hours)

                # Masking
                iOutputU=numpy.where(iOutputU < -100,9999,iOutputU)
                iOutputV=numpy.where(iOutputV < -100,9999,iOutputV)

                #Generates NC file
                sConvertedFilename= pt_folder + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc'
                f = netCDF4.Dataset(sConvertedFilename,"w","NETCDF3_CLASSIC")
                f.history = 'ERA Interim 10m winds -- Adapted to Medslik-II'
                f.createDimension('time',None)
                f.createDimension('lat', len(y1))
                f.createDimension('lon', len(x1))
                time = f.createVariable('time', 'd', ('time',))
                # Geo coordinates
                yCoord= f.createVariable('lat', 'd', ('lat',))
                xCoord= f.createVariable('lon', 'd', ('lon',))
                # Zonal velocity component
                uCoord= f.createVariable('U10M', 'd',('time','lat','lon'))
                vCoord= f.createVariable('V10M', 'd',('time','lat','lon'))
                time[:] = i + numpy.arange(1./24,24./24,1./24)
                #Filling the sausage
                yCoord[:]=y1
                xCoord[:]=x1
                #Filling the sausage chap. II
                uCoord[:]=iOutputU
                vCoord[:]=iOutputV
                #Adding units
                yCoord.units = 'degrees_north'
                xCoord.units = 'degrees_east'
                xCoord.lon_min = numpy.min(x1)
                yCoord.lat_min = numpy.min(y1)
                time.units = 'pythonic days'
                uCoord.units='m*s-1'
                uCoord.missing_value=9999
                vCoord.units='m*s-1'
                vCoord.missing_value=9999
                f.close()

def curr_gen_parallel(f1,f2,f3,rp_info, start_date,end_date,rawdata_folder,\
	itosra_dir, xp_name, expected_displacement, ds, input_ref_hours, output_ref_hours):

	rp_id = rp_info[0]
	rp_x = rp_info[1]
	rp_y = rp_info[2]
	output_folder = itosra_dir + '/' + xp_name + '/oce_files'

	# create output dirs
	os.system('mkdir ' + output_folder + '/H3k_' + "%05i" % (rp_id,))
	pt_folder= (output_folder + "/H3k_" + "%05i" % (rp_id,) + "/")

	# set the limits for the extraction of the current fields
	xmin = rp_x- expected_displacement*ds
	xmax = rp_x+ expected_displacement*ds
	ymin = rp_y- expected_displacement*ds
	ymax = rp_y+ expected_displacement*ds
	grid_corners = [xmin,xmax,ymin,ymax]


	counter = 0
	for i in range(start_date,end_date,1):

		iDateY=datetime.date.fromordinal(int(i-1))
		iDate=datetime.date.fromordinal(int(i))
		iDateT=datetime.date.fromordinal(int(i+1))

		if os.path.isfile(pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_U.nc'):
			print('... current file found. Next time step...')
			continue

		# Open MERCATOR files and interpolate in z for preset dephts 0,10,30,120m
		if counter == 0:
			print rawdata_folder + 'mercator_paria_20' + iDateY.strftime('%y') + iDateY.strftime('%m') + iDateY.strftime('%d') + '.nc'
			iLongitude,iLatitude,iU_Y,iV_Y,iT_Y,iTime_Y,iDepth = f1(rawdata_folder + 'mercator_paria_20' + iDateY.strftime('%y') + iDateY.strftime('%m') + iDateY.strftime('%d') + '.nc',grid_corners)
			iU_Y = (f2(iU_Y,iDepth))
			iV_Y = (f2(iV_Y,iDepth))
			iT_Y = (f2(iT_Y,iDepth))
			# day of interest
			iLongitude,iLatitude,iU_P,iV_P,iT_P,iTime,iDepth = f1(rawdata_folder + 'mercator_paria_20' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '.nc',grid_corners)
			iU_P = (f2(iU_P,iDepth))
			iV_P = (f2(iV_P,iDepth))
			iT_P = (f2(iT_P,iDepth))
			#  day after
			iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = f1(rawdata_folder + 'mercator_paria_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '.nc',grid_corners)
			iU_T = (f2(iU_T,iDepth))
			iV_T = (f2(iV_T,iDepth))
			iT_T = (f2(iT_T,iDepth))

			counter = counter + 1
		else:
			# Present becomes yesterday -- oh life...
			iU_Y = iU_P
			iV_Y = iV_P
			iT_Y = iT_P

			# Tomorrow becomes today -- oh dear life... time flies.
			iU_P = iU_T
			iV_P = iV_T
			iT_P = iT_T

			iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = f1(rawdata_folder + 'mercator_paria_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '.nc',grid_corners)
			iU_T = (f2(iU_T,iDepth))
			iV_T = (f2(iV_T,iDepth))
			iT_T = (f2(iT_T,iDepth))

		# Put the three fields (t - 1, t, t + 1) into a single matrix to interpolate
		iOutputUd = numpy.stack((iU_Y,iU_P,iU_T),axis=0)
		iOutputVd = numpy.stack((iV_Y,iV_P,iV_T),axis=0)
		iOutputTd = numpy.stack((iT_Y,iT_P,iT_T),axis=0)

		# Temporal interpolation:
		iOutputU = f3(iOutputUd, output_ref_hours, input_ref_hours)
		iOutputV = f3(iOutputVd, output_ref_hours, input_ref_hours)
		iOutputT = f3(iOutputTd, output_ref_hours, input_ref_hours)

		# Masking
		iOutputU=numpy.where(iOutputU < -3,9999,iOutputU)
		iOutputV=numpy.where(iOutputV < -3,9999,iOutputV)
		iOutputT=numpy.where(iOutputT < -4,9999,iOutputT)

		iOutputU=numpy.where(iOutputU > 3,9999,iOutputU)
		iOutputV=numpy.where(iOutputV > 3,9999,iOutputV)
		iOutputT=numpy.where(iOutputT > 40,9999,iOutputT)

		#Generates NC files
		# U component
		sConvertedFilename=pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_U.nc'
		f = netCDF4.Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
		f.history = 'MERCATOR Forecast - MyOcean -- Adapted to Medslik-II'
		f.createDimension('time_counter',None)
		f.createDimension('y', len(iLatitude))
		f.createDimension('x', len(iLongitude))
		f.createDimension('depthu',4)
		time = f.createVariable('time_counter', 'd', ('time_counter',))
		# Geo coordinates
		yCoord= f.createVariable('nav_lat', 'd', ('y'))
		xCoord= f.createVariable('nav_lon', 'd', ('x'))
		# Zonal velocity component
		uCoord= f.createVariable('vozocrtx', 'd',('time_counter','depthu','y','x'))
		time[:] = i + numpy.arange(1./24,25./24,1./24)
		#Filling the sausage
		yCoord[:]=iLatitude
		xCoord[:]=iLongitude
		#Filling the sausage chap. II
		uCoord[:]=iOutputU
		#Adding units
		yCoord.units = 'degrees_north'
		xCoord.units = 'degrees_east'
		time.units = 'julian days'
		uCoord.units='m*s-1'
		uCoord.missing_value=9999
		f.close()

		# V component
		sConvertedFilename=pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_V.nc'
		f = netCDF4.Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
		f.history = 'MERCATOR Forecast - MyOcean -- Adapted to Medslik-II'
		f.createDimension('time_counter',None)
		f.createDimension('y', len(iLatitude))
		f.createDimension('x', len(iLongitude))
		f.createDimension('depthv',4)
		time = f.createVariable('time_counter', 'd', ('time_counter',))
		# Geo coordinates
		yCoord= f.createVariable('nav_lat', 'd', ('y'))
		xCoord= f.createVariable('nav_lon', 'd', ('x'))
		# Zonal velocity component
		vCoord= f.createVariable('vomecrty', 'd',('time_counter','depthv','y','x'))
		time[:] = i + numpy.arange(1./24,25./24,1./24)
		#Filling the sausage
		yCoord[:]=iLatitude
		xCoord[:]=iLongitude
		vCoord[:]=iOutputV
		#Adding units
		yCoord.units = 'degrees_north'
		xCoord.units = 'degrees_east'
		time.units = 'julian days'
		vCoord.units='m*s-1'
		vCoord.missing_value=9999
		f.close()

		# Temperature
		sConvertedFilename=pt_folder + 'HOPS2km_H_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_T.nc'
		f = netCDF4.Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
		f.history = 'MERCATOR Forecast - MyOcean -- Adapted to Medslik-II v2.0'
		f.createDimension('time_counter',None)
		f.createDimension('y', len(iLatitude))
		f.createDimension('x', len(iLongitude))
		f.createDimension('deptht',4)
		time = f.createVariable('time_counter', 'd', ('time_counter',))
		# Geo coordinates
		yCoord= f.createVariable('nav_lat', 'd', ('y'))
		xCoord= f.createVariable('nav_lon', 'd', ('x'))
		# Zonal velocity component
		tCoord= f.createVariable('votemper', 'd',('time_counter','deptht','y','x'))
		time[:] = i + numpy.arange(1./24,25./24,1./24)
		#Filling the sausage
		yCoord[:]=iLatitude
		xCoord[:]=iLongitude
		tCoord[:]=iOutputT
		#Adding units
		yCoord.units = 'degrees_north'
		xCoord.units = 'degrees_east'
		time.units = 'julian days'
		tCoord.units='degree Celsius'
		tCoord.missing_value=9999
		f.close()

# POST PROCESSING SCRIPTS

def get_beached_oil(sNcFile,time_index,iTargetSites,iSegmentLengths):
	# EXTRACT VALUES OF INTEREST
	ncfile= netCDF4.Dataset(sNcFile,'r')
	oil_density = ncfile.variables['non_evaporative_volume'].oil_density
	parcel_volume = ncfile.variables['non_evaporative_volume'].volume_of_parcel
	rbm3=0.158987
	barrel2tonnes=1/(rbm3*(oil_density/1000))

	# EXTRACT VALUES FROM MEDSLIK OUTPUT NETCDF4 FILE
	lats = ncfile.variables['latitude'][time_index,:]
	lons = ncfile.variables['longitude'][time_index,:]
	particle_status = ncfile.variables['particle_status'][time_index,:]

	# PICKING BEACHED PARTICLES ONLY
	iNoise=numpy.argwhere(particle_status >= 0)
	lats = numpy.delete(lats, (iNoise), axis=1)
	lons = numpy.delete(lons, (iNoise), axis=1)
	particle_status = numpy.delete(particle_status, (iNoise), axis=1)
	iBeaching=(numpy.transpose(lons),numpy.transpose(lats),numpy.transpose(particle_status))

	iAssignedSegment=numpy.zeros(numpy.shape(lons)[1])
	iConcentrationsParcels=numpy.zeros(len(iTargetSites))
	iCP=numpy.zeros(len(iTargetSites))

	# ASSIGNING A COASTAL SEGMENT TO EACH BEACHED PARTICLE
	# for each beached parcel
	for jj in range(0,numpy.shape(lons)[1]):

		iParcelPosition=(iBeaching[0][jj],iBeaching[1][jj]) # x,y of the parcel
		iParcelDistance=ens_tools.dist(iParcelPosition,iTargetSites) # dist between parcel and all coastal segments
		iParcelDistance[numpy.isnan(iParcelDistance)]=9999 # border segmets are removed
		iClosestSegmentDist=numpy.min(iParcelDistance) # pick the closest segment
		# in case the distance is acceptable... assing a segment to the parcel
		if iClosestSegmentDist<.3/110:

			if len(numpy.argwhere(iParcelDistance==iClosestSegmentDist))>1:
				iAssignedSegment[jj]=numpy.argwhere(iParcelDistance==iClosestSegmentDist)[0]
			else:
				iAssignedSegment[jj]=numpy.argwhere(iParcelDistance==iClosestSegmentDist)

		# assing to the assigned segment a concentration
		iObservedOil=((-iBeaching[2][jj])/iSegmentLengths[int(iAssignedSegment[jj])])/barrel2tonnes
		iCP[int(iAssignedSegment[jj])]=iObservedOil
	ncfile.close()
	return iCP

def post_proc_file_gen(iReleasePoint,iModels,iMembers,time_index,\
	iStartDate,iEndDate,sModelFolder,iSegmentLengths,sOutfolder):

	sOutputFilename=[]
	iOctoputs=[]

	for iModel in iModels:

		for iMember in iMembers:

			jobs=[]

			for iDate in range(iStartDate,iEndDate,1):

				dDate=date.fromordinal(iDate)

				# DEFINE FOLDER IN WHICH THE SIMULATION OUTPUTS ARE
				sFolder=('RELOCATABLE_20' + dDate.strftime('%y') + '_' + dDate.strftime('%m') + \
					'_' + dDate.strftime('%d')  + '_0000_' + "%05d" % (iReleasePoint,) + '_' + \
					dDate.strftime('%d') + dDate.strftime('%m') + dDate.strftime('%y') + '_' + \
					str(iModel) + '_' + "%02d" % (iMember,))

				sNcFile = (sModelFolder + sFolder + '/spill_properties.nc')
				# CHECKING THE EXISTENCE OF THE FILE
				if os.path.isfile(sModelFolder + sFolder + '/medslik.fte')>0:

					try:
						print (sModelFolder + sFolder + '/spill_properties.nc')
						ncfile= Dataset(sModelFolder + sFolder + '/spill_properties.nc','r')
						oil_density = ncfile.variables['non_evaporative_volume'].oil_density

					except:
						print "file missing"

					jobs.append(job_server.submit(get_beached_oil,(sNcFile,time_index,iTargetSites,iSegmentLengths),(),("numpy","netCDF4","ens_tools","sys")))
					iCheck=iCheck+1

				else:
					continue
			if len(jobs)>0:
				sOutputFilename.append((sOutfolder + '/' + "%03d" % (time_index[0]+1,) + 'h_' + str(iModel) + '_' + dDate.strftime('%m') + dDate.strftime('%y') + '_'  + "%05d" % (iReleasePoint,) + '_' + "%02d" % (iMember,)))
				iOctoputs.append(jobs)


	print str(iCheck),' jobs related to release point', str(iReleasePoint), ' have been submitted. Harvesting results. This may take a while...'

	# Harvesting results

	jobs_index=0

	for jobs in iOctoputs:

		num_of_simulations=0

		for job in jobs:

			print sOutputFilename[jobs_index],num_of_simulations

			if num_of_simulations==0:
				iEnsembleCP=job()
			else:
				iEnsembleCP=np.vstack((iEnsembleCP,job()))
			num_of_simulations=num_of_simulations+1

		cPickle.dump([iEnsembleCP], open((sOutputFilename[jobs_index] + '.p'), 'wb'))
		jobs_index=jobs_index+1
	return



