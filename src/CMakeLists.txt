set(target_name dales4)

# Add all Fortran90 files from this directory
# This will also add some unwated files, which we will remove again below
FILE(GLOB sourcefiles "*.f90")

if(USE_HYPRE)
  # Build the iterative solver
  list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/modnohypre.f90)
  add_executable(${target_name} ${sourcefiles})
  target_link_libraries(${target_name} ${NETCDF_LIBS} ${HYPRE_LIB})
else(USE_HYPRE)
  # Do not build the iterative solver, but use (empty) fallback functions
  list(REMOVE_ITEM sourcefiles ${CMAKE_CURRENT_SOURCE_DIR}/modhypre.f90)
  add_executable(${target_name} ${sourcefiles})
  target_link_libraries(${target_name} ${NETCDF_LIBS})
endif(USE_HYPRE)

install(TARGETS dales4 DESTINATION ${CMAKE_BINARY_DIR})
