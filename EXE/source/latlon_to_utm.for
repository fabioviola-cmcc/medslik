c----------------------------------------------------------------------------
c  latlon_to_utm.for
c  ------------------
c  FORTRAN subroutine for coordinate conversion from lat,lon (WGS84) to UTM
-----------------------------------------------------------------------------
c  Copyright (C) 2013 Achilleas G. Samaras
c----------------------------------------------------------------------------
c  This program is free software: you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation, either version 3 of the License, or
c  any later version.
c  ------------------
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c  See the GNU General Public License for more details.
c  You should have received a copy of the GNU General Public License
c  along with this program.  If not, see <http://www.gnu.org/licenses/>.
c----------------------------------------------------------------------------
     
      subroutine latlon_to_utm(x,y,lat,lon,izone)

c----------------------------------------------------------------------------
c     IMPORTANT NOTES
c     ---------------
c     - Latitude and longitude must be given in decimal degrees
c     - Northern hemisphere: lat>0 / Southern hemisphere: lat<0
c     - Eastern hemisphere: lon>0 / Western hemisphere: lon<0
c     - If latitude is poleward of 84 degrees or dlam is not 
c       within 0.16 radians of the central meridian,
c       calculation will continue but output values may not be valid
c----------------------------------------------------------------------------

      implicit none

      real*8 lat,lon,ccdr,phi,lam,x,y
      real*8 a,b,f,fm,rm,k0,ecc,ecct2,n,rho,nu
      real*8 a0,b0,c0,d0,e0,s2,s4,s6,s8,scap
      real*8 dlam,s1
      real*8 cutm1,cutm2,cutm3,cutm4,cutm5

      integer izone,izcm


c     Zone and ZoneCM definition
c     ---------------------------

        if(lon.lt.0.d0) then
			izone=int((abs(lon)+180)/6)+1
		else
			izone=int((abs(lon))/6)+31
		endif
		
		izcm=6.d0*izone-183


c     Conversion factor from decimal degrees to radians
c     -------------------------------------------------
      data ccdr/57.29577951308232087d0/
      
c     Coordinates in radians
c     -------------------------------------------------
      phi = lat/ccdr
      lam = lon/ccdr


c     Datum Constants
c     ---------------

c     Equatorial radius
      data a/6378137.d0/
c     Polar radius
      data b/6356752.314245d0/
c     flattening
	  f=(a-b)/a	
c     inverse flatening
	  fm=1/f
c     mean radius
	  rm=dsqrt(a*b)
c     scale factor
      data k0/0.9996d0/
c     eccentricity squared
      ecc= dsqrt(1-(b/a)**2)
c     eccentricity primed squared                  
      ecct2 = ecc*ecc/(1.0d0-ecc*ecc)
c     n
      n=(a-b)/(a+b)
c     rho
      rho=(a*(1-ecc*ecc))/(1-(ecc*dsin(phi))**2)**(1.5)
c     nu      
      nu=a/dsqrt(1-(ecc*dsin(phi))**2)


c     Calculate Meridional Arc Length      
c     -------------------------------

      a0=a*(1-n+(5*n*n/4)*(1-n)+(81*n**4/64)*(1-n))
      b0=(3*a*n/2)*(1-n-(7*n*n/8)*(1-n)+55*n**4/64)
      c0=(15*a*n*n/16)*(1-n+(3*n*n/4)*(1-n))
      d0=(35*a*n**3/48)*(1-n+11*n*n/16)
      e0=(315*a*n**4/51)*(1-n)
      s2=dsin(2.0d0*phi)
      s4=dsin(4.0d0*phi)
      s6=dsin(6.0d0*phi)
      s8=dsin(8.0d0*phi)
      scap=a0*phi-b0*s2+c0*s4-d0*s6+e0*s8
      
      
c     Calculation constants  
c     ---------------------

 	  dlam=(lon-izcm)/ccdr
 	  s1=0
 	  
 	  
c     Coefficients for UTM coordinates  
c     --------------------------------

      cutm1=scap*k0
      cutm2=nu*dsin(phi)*dcos(phi)*k0/2  
      cutm3=((nu*dsin(phi)*dcos(phi)**3)/24) 
     &        *(5-dtan(phi)**2+9*ecct2*dcos(phi)**2
     &        +4*ecct2**2*dcos(phi)**4)*k0
      cutm4=nu*dcos(phi)*k0
      cutm5=(dcos(phi))**3*(nu/6)*(1-dtan(phi)**2
     &        +ecct2*dcos(phi)**2)*k0
	  
	  
c     UTM Coordinates
c     ---------------

c     Northing
      y=(cutm1+cutm2*dlam*dlam+cutm3*dlam**4)
c     Easting
      x=500000+(cutm4*dlam+cutm5*dlam**3)	  
	  
       
      return
      end 