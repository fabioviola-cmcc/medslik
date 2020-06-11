module load INTEL/intel_xe_2013
module load NETCDF/netcdf-4.2.1.1
module load HDF5/hdf5-1.8.10-patch1

echo $NETCDF

gfortran -o  ../Extract_II.exe Extract_II.for -I$NETCDF/include -L$NETCDF/lib  -lnetcdf -lnetcdff
