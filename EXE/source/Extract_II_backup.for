c-----------------------------------------------------------------------------------
c  MEDSLIK-II_1.01 
c  oil spill fate and transport model 
c-----------------------------------------------------------------------------------
c  Extract_II.for
c  This routine reads winds and currents from 
c  meteo-oceanogrpahic model output (NetCDF files)
c-----------------------------------------------------------------------------------
c  Copyright (C) <2012>
c  This program was originally written
c  by Robin Lardner and George Zodiatis.
c  Subsequent additions and modifications
c  have been made by Michela De Dominicis. 
c----------------------------------------------------------------------------------
c  The development of the MEDSLIK-II model is supported by a formal agreement
c  Memorandum of Agreement for the Operation and Continued Development of MEDSLIK-II
c  signed by the following institutions:
c  INGV - Istituto Nazionale di Geofisica e Vulcanologia
c  OC-UCY - Oceanography Center at the University of Cyprus
c  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per 
c  lo Studio dell’Ambiente Marino Costiero
c  CMCC - Centro Euro-Mediterraneo sui Cambiamenti Climatici
c 
c  This program is free software: you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation, either version 3 of the License, or
c  any later version.
c
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c  You should have received a copy of the GNU General Public License
c  along with this program.  If not, see <http://www.gnu.org/licenses/>.
c-----------------------------------------------------------------------------------

      character regn*4, indate(30)*8,indate_wind(30)*8, fc_dir*120
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       numfiles_wind,indate_wind,iviod, icurrents,fc_dir
      integer len_dir
c--------------------------------------------------------------------
c     read subregion limits & file dates. Adjust limits to lie on OPA grid
c--------------------------------------------------------------------
      call getarg(1,fc_dir)
      write(*,*) 'Directory = ', fc_dir 
      len_dir=120
      do while(fc_dir(len_dir:len_dir).eq.' ')
      len_dir=len_dir-1
      enddo


      open(1,file='medslik.tmp')
	read(1,*) regn, icurrents, iwind
	read(1,*) alon1,alon2
	read(1,*) alat1,alat2
	read(1,*) numfiles
	write(*,*) 'Wind type = ', iwind, numfiles
	write(*,*) indate_wind(1)

	do n=1,numfiles+1
	  read(1,'(a8)') indate(n)
	enddo

	read(1,*) numfiles_wind
	do n=1,numfiles_wind
	  read(1,'(a8)') indate_wind(n)

	enddo
	read(1,*) iviod
	close(1)

      open(99,file='Extract.log')
      if(icurrents.eq.74) call ExtractRELO(fc_dir,len_dir)
      if(iwind.eq.5)      call ExtractSKIRON(fc_dir,len_dir)

      stop
	end
c********************************************************************
c     Extract medslik files from Relocatable model output (1hr-3km)
c********************************************************************
      subroutine ExtractRELO(fc_dir,len_dir)

      parameter(ktmx=24, imx=121, jmx=121, kmx=4) 

c netcdf variables

      integer ncid, lenfile, lon_varid,lat_varid
      integer u_varid,v_varid,t_varid

      character*(*) LAT_NAME, LON_NAME
      character*(*) U_NAME, V_NAME, T_NAME
      
      parameter (LAT_NAME='nav_lat', LON_NAME='nav_lon')
      parameter (U_NAME='vozocrtx', V_NAME='vomecrty')
      parameter (T_NAME='votemper')
      
      integer lat_dimid, lon_dimid
      integer retval     
      
c other variables      

      real fmis !netcdf
      parameter(fmis=0.) !netcdf
      integer start(4), count(4) !netcdf  
      integer id, idU, idV, idT !netcdf
      integer Status !netcdf
      dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          ts(imx,jmx,kmx), u(imx,jmx,kmx), v(imx,jmx,kmx),
     &          ts_tmp(imx,jmx,kmx), u_tmp(imx,jmx,kmx),
     &          v_tmp(imx,jmx,kmx),              
     &          ts_24(imx,jmx,kmx,ktmx), u_24(imx,jmx,kmx,ktmx),
     &          v_24(imx,jmx,kmx,ktmx)
      character indate(30)*8, prdate*16, outfile*40,infile*120,                
     &          heads*150, empty*80, regn*4, ora*2, ore*2, fc_dir*120,
     &          start_lat*19
      logical ex
      integer t,len_dir
      common regn, alon1, alon2, alat1, alat2, numfiles, indate, 
     &       iviod, icurrents
      data udef /9999/,      rhoa /1.19/

c--------------------------------------------------------------------
c     main program
c--------------------------------------------------------------------
      do 60 n=1,numfiles
    
      Status = 0
      infile=fc_dir(1:len_dir)//'/fcst_data/H3k/'
     &                                //indate(n)(1:6)//'_U.nc'
      
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo
      
c++++++++++++++++++++++++      
c added by Augusto Neves between "+"
c simplifying reading the nc files
      if (n.eq.1) then

c open netcdf file      
      retval = nf_open(infile,0, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c inquire lat and lon variables
      retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

c load lat lons
      retval = nf_get_var_real(ncid, lon_varid, oplon)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_real(ncid, lat_varid, oplat)
      if (retval .ne. nf_noerr) call handle_err(retval)  
      
c generate vector lat lons
      
c      DO j=1,jmx
c	oplat(j) = lat2D(1,j)
c      ENDDO
                
c      DO i=1,imx
c        oplon(i) = lon2D(i,1)
c      ENDDO      
      
c estimate data resolution and min lat lon

      op_dlon = oplon(2) - oplon(1)
      op_dlat = oplat(2) - oplat(1)
      oplon0 = oplon(1)
      oplat0 = oplat(1)
      
      i_first = int( (alon1 - oplon0) / op_dlon ) + 1
      i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
      j_first = int( (alat1 - oplat0) / op_dlat ) + 1
      j_last  = int( (alat2 - oplat0) / op_dlat ) + 2

c cropping ur data

      IF (i_first.lt.1) i_first = 1
      IF (i_last.gt.imx) i_last = imx
      IF (j_first.lt.1) j_first = 1
      IF (j_last.gt.jmx) j_last = jmx

      alon1 = oplon0 + (i_first - 1) * op_dlon
      alon2 = oplon0 + (i_last  - 1) * op_dlon
      alat1 = oplat0 + (j_first - 1) * op_dlat
      alat2 = oplat0 + (j_last  - 1) * op_dlat

      imax = i_last - i_first + 1
      jmax = j_last - j_first + 1
      
      
      write(99,*) 'i-limits   = ',i_first,i_last,imax
      write(99,*) 'j-limits   = ',j_first,j_last,jmax
      write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
      write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16
      

      endif
c++++++++++++++++++++++++

c     set *_U.nc filename 
      infile=fc_dir(1:len_dir)//'/fcst_data/H3k/'
     &                                //indate(n)(1:6)//'_U.nc'

c     open netCDF file containing "u"
      retval = nf_open(infile,0, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, U_NAME, u_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     

c     extract data from u file
       do t = 1,ktmx
       start(4)=t

       retval = nf_get_var_real(ncid, u_varid,u_tmp)
       if (retval .ne. nf_noerr) call handle_err(retval)

       u_24(1:imx,1:jmx,1:kmx,t) = u_tmp(1:imx,1:jmx,1:kmx)

       enddo

       retval = nf_close (ncid)
       if (retval .ne. nf_noerr) call handle_err(retval)

c     set *_V.nc filename 
      infile=fc_dir(1:len_dir)//'/fcst_data/H3k/'
     &                                //indate(n)(1:6)//'_V.nc'

c     open netCDF file containing "v"
      retval = nf_open(infile,0, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, V_NAME, v_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     

c     extract data from v file
       do t = 1,ktmx
       start(4)=t

       retval = nf_get_var_real(ncid, v_varid,v_tmp)
       if (retval .ne. nf_noerr) call handle_err(retval)

       v_24(1:imx,1:jmx,1:kmx,t) = v_tmp(1:imx,1:jmx,1:kmx)

       enddo

       retval = nf_close (ncid)
       if (retval .ne. nf_noerr) call handle_err(retval)


c     set *_T.nc filename 
      infile=fc_dir(1:len_dir)//'/fcst_data/H3k/'
     &                                //indate(n)(1:6)//'_T.nc'

c     open netCDF file containing "temperature"
      retval = nf_open(infile,0, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, T_NAME, t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = kmx
      count(4) = 1     

c     extract data from t file
       do t = 1,ktmx
       start(4)=t

       retval = nf_get_var_real(ncid, t_varid,ts_tmp)
       if (retval .ne. nf_noerr) call handle_err(retval)

       ts_24(1:imx,1:jmx,1:kmx,t) = ts_tmp(1:imx,1:jmx,1:kmx)

       enddo

       retval = nf_close (ncid)
       if (retval .ne. nf_noerr) call handle_err(retval)



          
         do t=1,ktmx
          if(t.lt.10) then
           write(ore,'(i2)') t
	   ora='0'//ore(2:2)
	 else
	 write(ora,'(i2)') t
	 endif 
       
      
	   
	  prdate = indate(n)(5:6)//'/'//indate(n)(3:4)//'/20'//
     &                  indate(n)(1:2)//' '//ora//':00'
        write(6,*) 'Writing medslik file for date '//prdate
        write(99,*) 'Writing medslik file for date '//prdate
        outfile =
     & 'fcst_data/H3k//'//'relo'//indate(n)(1:6)//ora(1:2)//'.rel'

      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
c--------------------------------------------------------------------
c     read Relocatable data files
c--------------------------------------------------------------------
        u(1:imx,1:jmx,1:kmx) = u_24(1:imx,1:jmx,1:kmx,t)
	v(1:imx,1:jmx,1:kmx) = v_24(1:imx,1:jmx,1:kmx,t)
	ts(1:imx,1:jmx,1:kmx) = ts_24(1:imx,1:jmx,1:kmx,t)

c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(ts(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
       
	

        call extrap3d(ts, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(u, i1, i2, j1, j2, imx, jmx, kmx)
        call extrap3d(v, i1, i2, j1, j2, imx, jmx, kmx)
       
c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
c        do i=i_first,i_last
c        do j=j_first,j_last
c          if(msk(i,j).eq.1) then
c            do k=1,kmx
c              if(i.lt.imx) u(i,j,k) = (u(i,j,k) + u(i+1,j,k)) / 2.
c              if(j.lt.jmx) v(i,j,k) = (v(i,j,k) + v(i,j+1,k)) / 2.
c            enddo   
c          endif
c        enddo
c        enddo
	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/H3k//'//'relo'//indate(n)(1:6)//ora(1:2)//'.rel'
      open(20,file = outfile)
        write(20,*) 'MERCATOR model 9 km forecast data for '//prdate 
        write(20,*) 'Subregion of the Global Ocean:' 
c (AUGUSTO NEVES) one extra digit has been added
c to encompass southern hemisphere lat and longs
c     ORIGINAL
c        write(20,'(4f9.5,2i5,''   Geog. limits'')')
c     NEW
	write(20,'(4f10.5,2i5,''   Geog. limits'')')
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
        heads = '    lat        lon        SST        '//
     &  'u_srf      v_srf      u_10m      v_10m'//
     &  '       u_30m      v_30m      u_120m     v_120m'
        write(20,'(a150)') heads
        
	
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then

            blon = oplon(i)
            blat = oplat(j)
            sst = ts(i,j,1)
            us = u(i,j,1)
            vs = v(i,j,1)
            u10 = u(i,j,2)
            v10 = v(i,j,2)
            u10 = u(i,j,2)
            v10 = v(i,j,2)
            u30 = u(i,j,3)
            v30 = v(i,j,3) 
            u120 = u(i,j,4)
            v120 = v(i,j,4) 
c	      km = 1
c            do k=2,kmx
c              if(u(i,j,k).lt.udef.and.v(i,j,k).lt.udef) km=k
c            enddo
c
c            if(km.ge.5) then     
c              u10 = (u(i,j,4) * 1.559 + u(i,j,5) * 2.056 ) / 3.615 
c              v10 = (v(i,j,4) * 1.559 + v(i,j,5) * 2.056 ) / 3.615 
c            else
c              u10 = u(i,j,km)
c              v10 = v(i,j,km) 
c            endif
c            if(km.ge.10) then     
c              u30 = (u(i,j,9) * 4.164 + u(i,j,10) * 1.032 ) / 5.196 
c              v30 = (v(i,j,9) * 4.164 + v(i,j,10) * 1.032 ) / 5.196 
c            else
c              u30 = u(i,j,km)
c              v30 = v(i,j,km) 
c            endif
c            if(km.ge.21) then   
c              u120 = (u(i,j,20) * 3.459 + u(i,j,21) * 7.753 ) / 11.212 
c              v120 = (v(i,j,20) * 3.459 + v(i,j,21) * 7.753 ) / 11.212
c            else
c              u120 = u(i,j,km)
c              v120 = v(i,j,km) 
c            endif

            
            write(20,'(13f11.4)') blat,blon,sst,us,vs,
     &                                   u10,v10,u30,v30,u120,v120
            write (1,*) i_first,i_last,j_first,j_last
  
           endif
        enddo
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   
                  
c********************************************************************
c     Extract medslik files from SKIRON MFSTEP
c********************************************************************
      subroutine ExtractSKIRON(fc_dir,len_dir)
      

       parameter(ktmx=24, imx=257, jmx=113)

c netcdf variables

      integer ncid, lenfile, lon_varid,lat_varid
      integer u_varid,v_varid

      character*(*) LAT_NAME, LON_NAME
      character*(*) U_NAME, V_NAME
      
      parameter (LAT_NAME='lat', LON_NAME='lon')
      parameter (U_NAME='U10M', V_NAME='V10M')
        
      integer lat_dimid, lon_dimid
      integer retval     

c other variables
       real fmis !netcdf
       parameter(fmis=0.) !netcdf
       integer start(3), count(3) !netcdf  
       integer id, idWX, idWY !netcdf
       integer Status !netcdf
       dimension oplon(imx), oplat(jmx), msk(imx,jmx),
     &          wx(imx,jmx,ktmx),wy(imx,jmx,ktmx),
     &          wx_tmp(imx,jmx),wy_tmp(imx,jmx)

       character indate(30)*8, prdate*16, outfile*40, infile*120,
     &          heads*150, empty*80, regn*4, ora*2, ore*2,fc_dir*120,
     &          indate_wind(30)*8
       logical ex
       integer t,len_dir
       common regn, alon1, alon2, alat1, alat2, numfiles,indate,  
     &       numfiles_wind,indate_wind, iviod, icurrents
       data udef /9999./,      rhoa /1.19/

      INCLUDE 'netcdf.inc'
      PRINT *, NF_INQ_LIBVERS()
      
c--------------------------------------------------------------------
c  main program

      numfiles=numfiles_wind
      indate=indate_wind
            
      Status = 0
      
      do 60 n=1,numfiles   
      
      infile=fc_dir(1:len_dir)//'/fcst_data/SK1/'
     &                                //indate_wind(n)(1:8)//'.nc'
      
          
      len_file=120
      do while(infile(len_file:len_file).eq.' ')
      len_file=len_file-1
      enddo      
      
c++++++++++++++++++++++++      
c added by Augusto Neves between "+"
c simplifying reading the nc files
      if (n.eq.1) then

c open netcdf file      
      retval = nf_open(infile,0, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c inquire lat and lon variables
      retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

c load lat lons
      retval = nf_get_var_real(ncid, lon_varid, oplon)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_real(ncid, lat_varid, oplat)
      if (retval .ne. nf_noerr) call handle_err(retval)  
            
c estimate data resolution and min lat lon
      op_dlon = oplon(2) - oplon(1)
      op_dlat = oplat(1) - oplat(2) ! N to S

      oplon0 = oplon(1)
      oplat0 = oplat(1)
      
      i_first = int( (alon1 - oplon0) / op_dlon ) + 1
      i_last  = int( (alon2 - oplon0) / op_dlon ) + 2
      j_first = int( (oplat0 - alat2) / op_dlat ) + 1
      j_last  = int( (oplat0 - alat1) / op_dlat ) + 2
	
c cropping ur data
      IF (i_first.lt.1) i_first = 1
      IF (i_last.gt.imx) i_last = imx
      IF (j_first.lt.1) j_first = 1
      IF (j_last.gt.jmx) j_last = jmx
      
      alon1 = oplon0 + (i_first - 1) * op_dlon
      alon2 = oplon0 + (i_last  - 1) * op_dlon
      alat1 = oplat0 - (j_first - 1) * op_dlat
      alat2 = oplat0 - (j_last  - 1) * op_dlat

      imax = i_last - i_first + 1
      jmax = j_last - j_first + 1
	
      
      write(99,*) 'i-limits   = ',i_first,i_last,imax
      write(99,*) 'j-limits   = ',j_first,j_last,jmax
      write(99,*) 'lon-limits = ',alon1,alon2,(alon2-alon1)*16
      write(99,*) 'lat-limits = ',alat1,alat2,(alat2-alat1)*16
      

      endif
c++++++++++++++++++++++++      
c--------------------------------------------------------------------
c     read SKIRON data files
c--------------------------------------------------------------------

      ! Open *.nc File 

      Status = nf_open(infile(1:len_file),nf_nowrite,id)
 
      call handle_err(Status)
  
      Status = nf_inq_varid (id, 'U10M', idWX)
 
      call handle_err(Status)
  
      start(1) = 1
      start(2) = 1
      start(3) = 1
     
      count(1) = imx
      count(2) = jmx
      count(3) = 1
      
       do t = 1,ktmx
       start(3)=t
       Status = nf_get_vara_real (id, idWX, start, count, wx_tmp)
       call handle_err(Status)
       wx(1:imx,1:jmx,t) = wx_tmp(1:imx,1:jmx)
       enddo
             
      Status = nf_close (id)
       
       Status = nf_open(infile(1:len_file),nf_nowrite,id)
       call handle_err(Status)
       Status = nf_inq_varid (id, 'V10M', idWY)
       call handle_err(Status)
       start(1) = 1
       start(2) = 1
       start(3) = 1
        
       count(1) = imx
       count(2) = jmx
       count(3) = 1
      
       do t = 1,ktmx
       start(3)=t
       Status = nf_get_vara_real (id, idWY, start, count, wy_tmp)
       call handle_err(Status)
       wy(1:imx,1:jmx,t) = wy_tmp(1:imx,1:jmx)
       enddo

       Status = nf_close (id)
    
  	   
	  prdate = indate_wind(n)(7:8)//'/'//indate(n)(5:6)//'/20'//
     &                  indate(n)(3:4)
        write(6,*) 'Writing medslik SKIRON file for date '//prdate
        write(99,*) 'Writing medslik SKIRON file for date '//prdate
        outfile =
     & 'fcst_data/E25/'//'ecm_'//indate(n)(3:8)//'.ecm'
  
      
	  inquire(file = outfile, EXIST = ex)
	  if(ex) then
	    open(20,file = outfile)
          read(20,*) empty 
          read(20,*) empty 
c (AUGUSTO NEVES) one extra digit has been added
c to deal with southern hemisphere lat lon
c ORIGINAL
c          read(20,'(4f9.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
c NEW
          read(20,'(4f10.5,2i5)') blon1,blon2,blat1,blat2,imax1,jmax1
          if(blon1.eq.alon1.and.blon2.eq.alon2.and.blat1.eq.alat1.and.
     &       blat2.eq.alat2.and.imax1.eq.imax.and.jmax1.eq.jmax) then
            write(6,*) outfile//' already exists for this subregion'
            go to 60
          endif
          close(20)
        endif            
	
c--------------------------------------------------------------------
c     mask & nwp
c--------------------------------------------------------------------
        do i=1,imx
	  do j=1,jmx
	    msk(i,j) = 0
	    if(wx(i,j,1).lt.udef) msk(i,j) = 1
        enddo
        enddo

        nwp = 0
        do i=i_first,i_last
        do j=j_first,j_last
          if(msk(i,j).eq.1) then
            nwp = nwp + 1
          endif   
        enddo
        enddo
      
        write(99,*) 'nwp = ',nwp
        write(99,*) 'mask: '
        do j=j_last,j_first,-1
          write(99,'(300i1)') (msk(i,j),i=i_first,i_last)
        enddo
c        write(99,*) 'ts: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo

        i1 = i_first-2
        i2 = i_last+2
        j1 = j_first-2
        j2 = j_last+2
        if(i1.lt.1) i1 = 1 
        if(i2.gt.imx) i2 = imx 
        if(j1.lt.1) j1 = 1 
        if(j2.gt.jmx) j2 = jmx 
     
	
       call extrap2d(wx, i1, i2, j1, j2, imx, jmx)
       call extrap2d(wy, i1, i2, j1, j2, imx, jmx)


c        write(99,*) 'ts after exprapolation: '
c        do j=j_last,j_first,-1
c          write(99,*) (ts(i,j),i=i_first,i_last)
c        enddo
        
  	
c--------------------------------------------------------------------
c     write medslik files
c--------------------------------------------------------------------

        outfile =
     & 'fcst_data/SK1/'//'sk1_'//indate(n)(3:8)//'.sk1'
      open(20,file = outfile)
        write(20,*) 'SKIRON forecast data for '//prdate 
        write(20,*) 'Subregion of the Global Ocean with limits:' 
c (AUGUSTO NEVES) one extra digit has been added
c to deal with southern hemisphere lat lon
c ORIGINAL
c        write(20,'(4f9.5,2i5,''   Geog. limits'')') 
c NEW	
        write(20,'(4f10.5,2i5,''   Geog. limits'')') 
     &                                alon1,alon2,alat1,alat2,imax,jmax
        write(20,'(i6,''   0.0'')') nwp 
	heads = '    lat        lon        u00   '//
     &  '    v00      u01      v01      u02      v02'//
     &  '      u03     v03      u04      v04'//
     &  '      u05     v05      u06      v06'//
     &  '      u07     v07      u08      v08'//
     &  '      u09     v09      u10      v10'//
     &  '      u11     v11      u12      v12'//
     &  '      u13     v13      u14      v14'//
     &  '      u15     v15      u16      v16'//
     &  '      u17     v17      u18      v18'//
     &  '      u19     v19      u20      v20'//
     &  '      u21     v21      u22      v22'//
     &  '      u23     v23'  
        write(20,'(a900)') heads
c        print *, i_first, i_last, j_first, j_last
	
        do i=i_first,i_last
        do j=j_first,j_last
c         do j=j_last,j_first,-1
c           if(msk(i,j).eq.1) then

	    blon = oplon(i)
            blat = oplat(j)
            
c	    wx0 = wx(i,j,1)
c	    wy0 = wy(i,j,1)
c            endif

	     
              wx00 = wx(i,j,1)
              wy00 = wy(i,j,1)
              wx01 = wx(i,j,2)
              wy01 = wy(i,j,2)
              wx02 = wx(i,j,3)
              wy02 = wy(i,j,3)
              wx03 = wx(i,j,4)
              wy03 = wy(i,j,4)
              wx04 = wx(i,j,5)
              wy04 = wy(i,j,5)
              wx05 = wx(i,j,6)
              wy05 = wy(i,j,6)
              wx06 = wx(i,j,7)
              wy06 = wy(i,j,7)
              wx07 = wx(i,j,8)
              wy07 = wy(i,j,8)
              wx08 = wx(i,j,9)
              wy08 = wy(i,j,9)
              wx09 = wx(i,j,10)
              wy09 = wy(i,j,10)
              wx10 = wx(i,j,11)
              wy10 = wy(i,j,11)
              wx11 = wx(i,j,12)
              wy11 = wy(i,j,12)
              wx12 = wx(i,j,13)
              wy12 = wy(i,j,13)
              wx13 = wx(i,j,14)
              wy13 = wy(i,j,14)
              wx14 = wx(i,j,15)
              wy14 = wy(i,j,15)
              wx15 = wx(i,j,16)
              wy15 = wy(i,j,16)
              wx16 = wx(i,j,17)
              wy16 = wy(i,j,17)
              wx17 = wx(i,j,18)
              wy17 = wy(i,j,18)
              wx18 = wx(i,j,19)
              wy18 = wy(i,j,19)
              wx19 = wx(i,j,20)
              wy19 = wy(i,j,20)
              wx20 = wx(i,j,21)
              wy20 = wy(i,j,21)
              wx21 = wx(i,j,22)
              wy21 = wy(i,j,22)
              wx22 = wx(i,j,23)
              wy22 = wy(i,j,23)
              wx23 = wx(i,j,24)
              wy23 = wy(i,j,24)
           
             
            write(20,'(50f11.4)') blat,blon,wx00,wy00,wx01,wy01
     &      ,wx02,wy02,wx03,wy03,wx04,wy04,wx05,wy05,wx06,wy06
     &      ,wx07,wy07,wx08,wy08,wx09,wy09,wx10,wy10,wx11,wy11  
     &      ,wx12,wy12,wx13,wy13,wx14,wy14,wx15,wy15,wx16,wy16 
     &      ,wx17,wy17,wx18,wy18,wx19,wy19,wx20,wy20,wx21,wy21
     &      ,wx22,wy22,wx23,wy23
        
        enddo
        enddo
        close(20)
   60 continue
        
      return
      end   

c************************************************************************
c      Extrapolation of 2-D fields over land points
c------------------------------------------------------------------------

      subroutine extrap2d(data, i_first, i_last, j_first, j_last, 
     &                          imx, jmx)

      real data(imx,jmx),carpet(imx,jmx)
    
      data udef /9999/

	ngridpts = (i_last - i_first + 1) * (j_last - j_first + 1)

      do iter=1,20
        knt=0
        do j = j_first, j_last
        do i = i_first, i_last

          if(data(i,j). ge. udef) then
             knt = knt + 1
               im1=i-1
               ip1=i+1
               jm1=j-1
               jp1=j+1
               if(im1.lt.i_first) im1 = i_first
               if(ip1.gt.i_last)  ip1 = i_last
               if(jm1.lt.j_first) jm1 = j_first
               if(jp1.gt.j_last)  jp1 = j_last

               datan=0.
               jcn=0
               do jj=jm1,jp1
               do ii=im1,ip1

                    if(data(ii,jj). lt. udef) then
                      datan = datan + data(ii,jj)
                      jcn = jcn + 1
                   end if
               enddo 
               enddo 

               if(jcn. gt. 2) then
                 data(i,j) = datan / float(jcn)
               end if

          end if

        enddo 
        enddo
c        write(99,*) iter,knt
        if(knt.eq.0.or.knt.eq.ngridpts) return 
        
      enddo
                    
      return
      end      

c************************************************************************
c      Extrapolation of 3-D fields over land points, 7 levels only
c------------------------------------------------------------------------

      subroutine extrap3d (data, i_first, i_last, j_first, j_last, 
     &                           imx, jmx, kmx)


      real data(imx,jmx,kmx), tmp(imx,jmx)
      integer kd(7)

      data udef /9999/
      data kd /1,3,4,8,9,19,20/
c
      do k=1,kmx
c        k = kd(n)

        do i=1,imx
        do j=1,jmx
          tmp(i,j) = data(i,j,k)
        enddo
        enddo
        
        call extrap2d(tmp, i_first, i_last, j_first, j_last, imx, jmx)

        do i=1,imx
        do j=1,jmx
          data(i,j,k) = tmp(i,j)
        enddo
        enddo
        
      enddo   !k

      return
      end      

c****************************************************************************
c****************************************************************************
c****************************************************************************
      SUBROUTINE HANDLE_ERR(Status)
      include 'netcdf.inc'
c      include '/usr/local/include/netcdf.inc'
      INTEGER Status
      IF (Status .NE. NF_NOERR) THEN
        PRINT *, NF_STRERROR(Status)
        STOP 'Stopped'
      ENDIF
      END
