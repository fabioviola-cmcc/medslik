module load INTEL/intel_xe_2013
module load NETCDF/netcdf-4.3
module load HDF5/hdf5-1.8.10-patch1
gfortran -o medslik_II.exe medslik_II.for -I$NETCDF/include -L$NETCDF/lib  -lnetcdf -lnetcdff

