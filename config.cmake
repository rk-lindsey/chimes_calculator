#
# User adjustable config options
#

# Debug level
set(DEBUG "0" CACHE STRING "Debug level")

# Enable if you want to build a Fortran interface
option(WITH_FORTRAN_API "Whether the Fortran API should be built" TRUE)

# Enable if you want to build a Fortran 2008 style object oriented interface
option(WITH_FORTRAN08_API "Whether the Fortran 2008 style API should be built" TRUE)

# Turn this on, if the libraries should be built as shared libraries
option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" TRUE)
