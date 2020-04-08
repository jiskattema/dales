# Juwels GCC with parastation  - v. 2019
# use with: 
# load modules:
#  module use /gpfs/software/juwels/otherstages
#  module load Stages/2018b
# module load GCC/8.2.0 # module load Intel/2019.0.117-GCC-7.3.0
# module load ParaStationMPI/5.2.1-1  # module load ParaStationMPI
# module load netCDF-Fortran/4.4.4-serial  # module load HDF
# module load HDF    # module load netCDF
# module load netCDF # module load netCDF-Fortran/4.4.4
# module load CMake # module load CMake
#
#
set(CMAKE_Fortran_COMPILER "gfortran" ) # "ifort") #"gfortran")
set(Fortran_COMPILER_WRAPPER mpif90 ) # mpiifort)

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

# set(CMAKE_Fortran_COMPILER "ifort" ) # "gfortran" ) # "ifort") #"gfortran")
# set(Fortran_COMPILER_WRAPPER mpif90 ) # mpiifort)
# 
# set(USER_Fortran_FLAGS "-autodouble -unroll-loops")
# set(USER_Fortran_FLAGS_RELEASE "-autodouble  -unroll-loops -O3")
# set(USER_Fortran_FLAGS_DEBUG "-autodouble -unroll-loops -O3")

set(NETCDF_INCLUDE_DIR "/gpfs/software/juwels/stages/2018b/software/netCDF-Fortran/4.4.4-GCC-8.2.0-serial/include/")
set(NETCDF_LIB_1      "/gpfs/software/juwels/stages/2018b/software/netCDF/4.6.1-GCC-8.2.0-serial/lib64/libnetcdf.so")
set(NETCDF_LIB_2      "/gpfs/software/juwels/stages/2018b/software/netCDF-Fortran/4.4.4-GCC-8.2.0-serial/lib/libnetcdff.so")
set(HDF5_LIB_1        "/gpfs/software/juwels/stages/2018b/software/HDF/4.2.13-GCC-8.2.0/lib/libdf.a")
set(HDF5_LIB_2        "/gpfs/software/juwels/stages/2018b/software/HDF/4.2.13-GCC-8.2.0/lib/libmfhdf.a")
 
# set(NETCDF_INCLUDE_DIR "/gpfs/software/juwels/stages/2018b/software/netCDF-Fortran/4.4.4-iccifort-2019.0.117-GCC-7.3.0-serial/include") # /rrzk/lib/netcdf/4.1.3-gcc-4.8.2/include") # set(NETCDF_INCLUDE_DIR "/opt/rrzk/lib/netcdf/4.1.3/include")
# set(NETCDF_LIB_1       "/gpfs/software/juwels/stages/2018b/software/netCDF/4.6.1-iccifort-2019.0.117-GCC-7.3.0-serial/lib/libnetcdf.so") # "/opt/rrzk/lib/netcdf/4.1.3-gcc-4.8.2/lib/libnetcdff.so") # set(NETCDF_LIB_1       "/opt/rrzk/lib/netcdf/4.1.3/lib/libnetcdff.so")
# set(NETCDF_LIB_2       "/gpfs/software/juwels/stages/2018b/software/netCDF-Fortran/4.4.4-iccifort-2019.0.117-GCC-7.3.0-serial/lib/libnetcdff.so") # "/opt/rrzk/lib/netcdf/4.1.3-gcc-4.8.2/lib/libnetcdf.so") # set(NETCDF_LIB_2       "/opt/rrzk/lib/netcdf/4.1.3/lib/libnetcdf.so")
# set(NETCDF_LIB_3       "/opt/rrzk/lib/netcdf/netcdf-4.1.3/f90/netcdf.mod")
# set(HDF5_LIB_1         "/gpfs/software/juwels/stages/2018b/software/HDF/4.2.13-iccifort-2019.0.117-GCC-7.3.0/lib/libdf.a") #  "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5_hl.so")  # set(HDF5_LIB_1         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5_hl.so")
# set(HDF5_LIB_2         "/gpfs/software/juwels/stages/2018b/software/HDF/4.2.13-iccifort-2019.0.117-GCC-7.3.0/") #  "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5.so") # set(HDF5_LIB_2         "/opt/rrzk/lib/hdf5/1.8.11/lib/libhdf5.so")

set(SZIP_LIB           "") #   "/opt/rrzk/lib/szip/szip-2.1/lib/libsz.so")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2}  ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
# set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${NETCDF_LIB_3} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)  
#



