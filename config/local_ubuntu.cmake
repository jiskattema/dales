# Ubuntu Linux
# cmake setting for computer: raou
# --------------------------------------
#
# - check whether these libraries are installed:
#       libnetcdf-dev
#       libnetcdff-dev
#       libhdf5-dev or  libhdf4-dev
#       curl
#       z
#             
#        inf not ,l install them and adjust paths accordingly
#
set(CMAKE_Fortran_COMPILER "gfortran")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

set(NETCDF_INCLUDE_DIR "/usr/include")
# changing from /lib/ to /lib/x86_64-linux-gnu/
set(NETCDF_LIB_1       "/usr/lib/x86_64-linux-gnu/libnetcdff.a")
# replacing set(NETCDF_LIB_2       "/usr/lib/libnetcdf.a")
set(NETCDF_LIB_2       "/usr/lib/x86_64-linux-gnu/libnetcdf.so")
#
# set(HDF5_LIB_1         "/usr/lib/libhdf5_hl.a")
# replacing it with
# ./lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5_hl.a
# OR 
# ./lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a
set(HDF5_LIB_1         "/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a")
#
# similarly for:
# set(HDF5_LIB_2         "/usr/lib/libhdf5.a")
set(HDF5_LIB_2         "/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a")
#
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
