########################################
# Set up how to compile the source files
########################################

# Define the executable name
SET(1DPMEEXE 1dpme.exe)

# Add the source files
SET(1DPME_Sources ${SRC}/main.f90
                  ${SRC}/linear_solvers.f90
                  ${SRC}/special_functions.f90
 )

# Define the executable in terms of the source files
ADD_EXECUTABLE(${1DPMEEXE} ${1DPME_Sources})

# This links foo to the bar library
#TARGET_LINK_LIBRARIES(${1DPMEEXE} ${SOME_LIBRARY})

#####################################
# Tell how to install this executable
#####################################

set_target_properties( ${1DPMEEXE}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
