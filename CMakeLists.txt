cmake_minimum_required(VERSION 3.13)

############################################################
# General settings
############################################################

project(ChimesCalc
    VERSION 1.0
    DESCRIPTION "ChIMES calculator: Utilities to compute ChIMES interactions"
	LANGUAGES CXX C
)

string(TOLOWER "${PROJECT_NAME}" PROJECT_NAME_LOWER)

include(GNUInstallDirs)
set(INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}")
set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/modfiles")

include(config.cmake)

if(WITH_FORTRAN_API OR WITH_FORTRAN08_API)
    enable_language(Fortran)
endif()


# Check whether libm is needed to provide the round() function in C
include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS("m" "round" "" NEEDS_LIB_M)
if(NEEDS_LIB_M)
  set(LIB_M "m")
endif()


###############################################################################
# C++/C library
###############################################################################

set(c-sources
    "chimesFF/src/chimesFF.cpp"
    "chimesFF/api/chimescalc_C.cpp"
    "serial_interface/src/serial_chimes_interface.cpp"
    "serial_interface/api/chimescalc_serial_C.cpp"
)
set(c-includes
    "chimesFF/src/chimesFF.h"
    "chimesFF/api/chimescalc_C.h"
    "serial_interface/src/serial_chimes_interface.h"
    "serial_interface/api/chimescalc_serial_C.h"
)

add_library(ChimesCalc "${c-sources}")

set_target_properties(ChimesCalc
    PROPERTIES
    OUTPUT_NAME "chimescalc"
    PUBLIC_HEADER "${c-includes}"
)
target_include_directories(ChimesCalc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>
)
target_compile_definitions(ChimesCalc PRIVATE "DEBUG=${DEBUG}")
target_compile_features(ChimesCalc PRIVATE cxx_std_11)

install(TARGETS ChimesCalc
    EXPORT ${PROJECT_NAME_LOWER}-targets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_INCLUDEDIR}"
)


############################################################
# Executables for C/C++
############################################################

add_executable(CPP-interface serial_interface/examples/cpp/main.cpp)
target_link_libraries(CPP-interface ChimesCalc)
target_compile_features   (CPP-interface PUBLIC cxx_std_11)
target_compile_options    (CPP-interface PUBLIC -O3)
target_compile_definitions(CPP-interface PRIVATE "DEBUG=${DEBUG}")


IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	install(TARGETS CPP-interface DESTINATION ${CMAKE_INSTALL_PREFIX}/serial_interface/examples/cpp/ OPTIONAL)
ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

add_executable(C_wrapper-direct_interface chimesFF/examples/c/main.c)
target_link_libraries     (C_wrapper-direct_interface ChimesCalc ${LIB_M})
target_compile_features   (C_wrapper-direct_interface PRIVATE cxx_std_11)
target_compile_options    (C_wrapper-direct_interface PRIVATE -O3)
target_compile_definitions(C_wrapper-direct_interface PRIVATE "DEBUG=${DEBUG}")


IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	install(TARGETS C_wrapper-direct_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/chimesFF/examples/c/ OPTIONAL)
ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


add_executable(C_wrapper-serial_interface serial_interface/examples/c/main.c)
target_link_libraries     (C_wrapper-serial_interface ChimesCalc)
target_include_directories(C_wrapper-serial_interface PRIVATE serial_interface/api/)
target_compile_features   (C_wrapper-serial_interface PRIVATE cxx_std_11)
target_compile_options    (C_wrapper-serial_interface PRIVATE -O3)
target_compile_definitions(C_wrapper-serial_interface PRIVATE "DEBUG=${DEBUG}")

IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	install(TARGETS C_wrapper-serial_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/serial_interface/examples/c/ OPTIONAL)
ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


############################################################
# Executables/libraries for Python
############################################################


# Shared: .dylib
# Static: .a
# Module: .so

add_library(lib-C_wrapper-direct_interface MODULE
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/chimescalc_C.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/chimesFF.cpp
)

set_target_properties     (lib-C_wrapper-direct_interface PROPERTIES PREFIX "") # Prevent CMake from prepending "lib" to libwrapper-C.so
set_target_properties     (lib-C_wrapper-direct_interface PROPERTIES LINKER_LANGUAGE CXX)  # Tell it which language to use
target_include_directories(lib-C_wrapper-direct_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/)
target_include_directories(lib-C_wrapper-direct_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/)
target_compile_features   (lib-C_wrapper-direct_interface PUBLIC cxx_std_11)
target_compile_options    (lib-C_wrapper-direct_interface PUBLIC -O3)
target_compile_options    (lib-C_wrapper-direct_interface PUBLIC -fPIC)

IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	install(TARGETS lib-C_wrapper-direct_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/chimesFF/examples/python/ OPTIONAL)
ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


add_library(lib-C_wrapper-serial_interface MODULE
    ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/api/chimescalc_serial_C.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src/serial_chimes_interface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/chimesFF.cpp
)

set_target_properties     (lib-C_wrapper-serial_interface PROPERTIES PREFIX "") # Prevent CMake from prepending "lib" to libwrapper-C.so
set_target_properties     (lib-C_wrapper-serial_interface PROPERTIES LINKER_LANGUAGE CXX)  # Tell it which language to use
target_include_directories(lib-C_wrapper-serial_interface PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/api/)
target_include_directories(lib-C_wrapper-serial_interface PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src/)
target_include_directories(lib-C_wrapper-serial_interface PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/)
target_compile_features   (lib-C_wrapper-serial_interface PUBLIC cxx_std_11)
target_compile_options    (lib-C_wrapper-serial_interface PUBLIC -O3)
target_compile_options    (lib-C_wrapper-serial_interface PUBLIC -fPIC)

#install(TARGETS lib-C_wrapper-serial_interface DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/examples/python/)
IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	install(TARGETS lib-C_wrapper-serial_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/serial_interface/examples/python/ OPTIONAL)
ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)



###############################################################################
# Fortran library (separate library as it needs Fortran runtime to be linked)
###############################################################################

if(WITH_FORTRAN_API OR WITH_FORTRAN08_API)

    set(f90safe-fortran-sources
        "chimesFF/api/chimescalc_F.f90"
        "serial_interface/api/chimescalc_serial_F.f90"
    )

    set(fortran-sources
        ${f90safe-fortran-sources}
    )

    # Optional, as ancient Fortran compilers may have difficulties with it
    if(WITH_FORTRAN08_API)
        list(APPEND fortran-sources
            "serial_interface/api/chimescalc_serial_F08.f90"
        )
    endif()

    add_library(ChimesCalc_Fortran "${fortran-sources}")
    set_target_properties(ChimesCalc_Fortran PROPERTIES OUTPUT_NAME "chimescalc_fortran")
    target_link_libraries(ChimesCalc_Fortran PUBLIC ChimesCalc)

    set(moddir "${CMAKE_CURRENT_BINARY_DIR}/modfiles")
    set_target_properties(ChimesCalc_Fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${moddir}"
    )
    target_include_directories(ChimesCalc_Fortran PUBLIC
        $<BUILD_INTERFACE:${moddir}>
        $<INSTALL_INTERFACE:${INSTALL_MODULEDIR}>
    )

    install(TARGETS ChimesCalc_Fortran
        EXPORT ${PROJECT_NAME_LOWER}-targets
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    )
    install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MODULEDIR}")

    ############################################################
   # Executables for Fortran (90/08)
   ############################################################


    add_executable(fortran_wrapper-direct_interface chimesFF/examples/fortran/main.F90)
    target_link_libraries     (fortran_wrapper-direct_interface ChimesCalc_Fortran)
    target_compile_options    (fortran_wrapper-direct_interface PRIVATE -std=f2003)
    target_compile_definitions(fortran_wrapper-direct_interface PRIVATE "DEBUG=${DEBUG}")

    IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    	    install(TARGETS fortran_wrapper-direct_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/chimesFF/examples/fortran/ OPTIONAL)
    ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


    add_executable(fortran_wrapper-serial_interface serial_interface/examples/fortran/main.F90)
    target_link_libraries     (fortran_wrapper-serial_interface ChimesCalc_Fortran)
    target_compile_options    (fortran_wrapper-serial_interface PRIVATE -std=f2003)
    target_compile_definitions(fortran_wrapper-serial_interface PRIVATE "DEBUG=${DEBUG}")

    IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    	    install(TARGETS fortran_wrapper-direct_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/serial_interface/examples/fortran/ OPTIONAL)
    ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

    if(WITH_FORTRAN08_API) # Optional, as ancient Fortran compilers may have difficulties with it
        add_executable(fortran08_wrapper-serial_interface serial_interface/examples/fortran08/main.F90)
        target_link_libraries	  (fortran08_wrapper-serial_interface ChimesCalc_Fortran)
	target_compile_options    (fortran08_wrapper-serial_interface PRIVATE -std=f2008)
	target_compile_definitions(fortran08_wrapper-serial_interface PRIVATE "DEBUG=${DEBUG}")

	IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	       install(TARGETS fortran08_wrapper-serial_interface DESTINATION ${CMAKE_INSTALL_PREFIX}/serial_interface/examples/fortran08/ OPTIONAL)
	ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    endif()

endif()

###############################################################################
# Create CMake export file
###############################################################################

include(CMakePackageConfigHelpers)

install(
    EXPORT ${PROJECT_NAME_LOWER}-targets
    FILE "${PROJECT_NAME_LOWER}-targets.cmake"
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}"
)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/export/${PROJECT_NAME_LOWER}-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}"
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config-version.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}"
)
