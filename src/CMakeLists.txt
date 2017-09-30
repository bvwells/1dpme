########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(1DPME_src ${SRC}/main.f90
              ${SRC}/linear_solvers.f90
              ${SRC}/special_functions.f90
 )

# Define the executable in terms of the source files
ADD_EXECUTABLE(${1DPMEEXE} ${1DPME_src})

# This links foo to the bar library
#TARGET_LINK_LIBRARIES(${1DPMEEXE} ${SOME_LIBRARY})

#####################################
# Tell how to install this executable
#####################################

INSTALL(TARGETS ${1DPMEEXE} RUNTIME DESTINATION ${BIN})

