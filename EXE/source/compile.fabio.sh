# load modules
module load intel19.5/19.5.281
module load intel19.5/netcdf/C_4.7.2-F_4.5.2_CXX_4.3.1
module load intel19.5/hdf5/1.10.5
module load intel19.5/szip/2.1.1
module load curl/7.66.0

#export CC=icc
#export F9X=ifort
#export CXX=icpc

#export LC_ALL=en_US

#gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for

# set folders
wdir=/work/opa/witoil-dev/medslik-sanifs/meglob-fabio/
DIR_EXE=${wdir}/EXE
DIR_SRC=${wdir}/EXE/source
NETCDF=/zeus/opt/intel19.5/netcdf

################sergio vecchie
ifort -O3 -o $DIR_EXE/Extract_II.exe $DIR_SRC/Extract_II.for -I$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L/$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/lib -lnetcdff -lnetcdf -fPIC -shared-intel -mcmodel=large
#gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f
#ifort -O3 -o $DIR_EXE/medslik_II.exe $DIR_SRC/medslik_II.for

# compile
# gfortran -I$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/lib  $DIR_SRC/Extract_II.for -lnetcdf -lnetcdff -o $DIR_EXE/Extract_II.exe -fPIC -mcmodel=large -g -fcheck=all -fbacktrace

# This line is commented because I already have another file to generate this
#gfortran  -I/usr/local/netcdf/current/include -L/usr/local/netcdf/current/lib  $DIR_SRC/Extract_II.for -lnetcdf -lnetcdff -o $DIR_EXE/Extract_II.exe

#gfortran -I$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/lib  $DIR_SRC/Extract_II.for -lnetcdf -lnetcdff -o $DIR_EXE/Extract_II.exe

gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f 

# gfortran  -I/usr/local/netcdf/current/include -L/usr/local/netcdf/current/lib  $DIR_SRC/medslik_II.for -lnetcdf -lnetcdff -o $DIR_EXE/medslik_II.exe
gfortran -I$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/include -L$NETCDF/C_4.7.2-F_4.5.2_CXX_4.3.1/lib  $DIR_SRC/medslik_II.for -lnetcdf -lnetcdff -o $DIR_EXE/medslik_II.exe


gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for
