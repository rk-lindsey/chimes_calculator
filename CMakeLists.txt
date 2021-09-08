cmake_minimum_required(VERSION 3.10)

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


###############################################################################
# C++/C library
###############################################################################

set(c-sources
    "chimesFF/src/chimescalc.cpp"
    "chimesFF/api/chimescalc_C.cpp"
    "serial_interface/src/chimescalc_serial.cpp"
    "serial_interface/api/chimescalc_serial_C.cpp"
)
set(c-includes
    "chimesFF/src/chimescalc.h"
    "chimesFF/api/chimescalc_C.h"
    "serial_interface/src/chimescalc_serial.h"
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


###############################################################################
# Fortran library (separate library as it needs Fortran runtime to be linked)
###############################################################################

if(WITH_FORTRAN_API OR WITH_FORTRAN08_API)

    set(fortran-sources
        "chimesFF/api/chimescalc_F.f90"
        "serial_interface/api/chimescalc_serial_F.f90"
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
