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


# C++ compiler dependent config options
if("GNU" STREQUAL "${CMAKE_CXX_COMPILER_ID}")

    set(CXX_FLAGS "${CMAKE_CXX_FLAGS}"
        CACHE STRING "Extra flags for the C compiler")

    set(CXX_FLAGS_RELEASE "-O3 -funroll-all-loops"
        CACHE STRING "C compiler flags for Release build")

    set(CXX_FLAGS_DEBUG "-g -Wall -pedantic -fbounds-check"
        CACHE STRING "C compiler flags for Debug build")

elseif("Intel" STREQUAL "${CMAKE_CXX_COMPILER_ID}")

    set(CXX_FLAGS "${CMAKE_CXX_FLAGS}"
        CACHE STRING "Extra flags for the C compiler")

    set(CXX_FLAGS_RELEASE "-O3 -ip"
        CACHE STRING "C compiler flags for Release build")

    set(CXX_FLAGS_DEBUG "-g -Wall -traceback"
        CACHE STRING "C compiler flags for Debug build")

endif()


# C compiler dependent config options
if("GNU" STREQUAL "${CMAKE_C_COMPILER_ID}")

    set(C_FLAGS "${CMAKE_C_FLAGS}"
        CACHE STRING "Extra flags for the C compiler")

    set(C_FLAGS_RELEASE "-O3 -funroll-all-loops"
        CACHE STRING "C compiler flags for Release build")

    set(C_FLAGS_DEBUG "-g -Wall -pedantic -fbounds-check"
        CACHE STRING "C compiler flags for Debug build")

elseif("Intel" STREQUAL "${CMAKE_C_COMPILER_ID}")

    set(C_FLAGS "${CMAKE_C_FLAGS}"
        CACHE STRING "Extra flags for the C compiler")

    set(C_FLAGS_RELEASE "-O3 -ip"
        CACHE STRING "C compiler flags for Release build")

    set(C_FLAGS_DEBUG "-g -Wall -traceback"
        CACHE STRING "C compiler flags for Debug build")

endif()

# Compiler detection, so that CMAKE_Fortran_COMPILER_ID is defined below
if(WITH_FORTRAN_API OR WITH_FORTRAN08_API)
  enable_language(Fortran)
endif()

# Fortran compiler dependent config options
if("GNU" STREQUAL "${CMAKE_Fortran_COMPILER_ID}")

    set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008"
        CACHE STRING "Extra flags for the Fortran compiler")

    set(Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops"
        CACHE STRING "Fortran compiler flags for Release build")

    set(Fortran_FLAGS_DEBUG "-g -Wall -pedantic -fbounds-check"
        CACHE STRING "Fortran compiler flags for Debug build")

elseif("Intel" STREQUAL "${CMAKE_Fortran_COMPILER_ID}")

    set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f08"
        CACHE STRING "Extra flags for the Fortran compiler")

    set(Fortran_FLAGS_RELEASE "-O3 -ip"
        CACHE STRING "Fortran compiler flags for Release build")

    set(Fortran_FLAGS_DEBUG "-g -warn all -check -traceback"
        CACHE STRING "Fortran compiler flags for Debug build")

endif()
