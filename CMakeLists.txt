cmake_minimum_required(VERSION 3.13)

####################################################################################################
# General settings
####################################################################################################

project(ChimesCalc
    VERSION 1.0
    DESCRIPTION "ChIMES calculator: Utilities to compute ChIMES interactions"
  	LANGUAGES CXX C)

string(TOLOWER "${PROJECT_NAME}" PROJECT_NAME_LOWER)

include(GNUInstallDirs)
set(INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME_LOWER}")
set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/modfiles")

include(cmake/ChimesCalc.cmake)
include(config.cmake)

chimescalc_setup_build_type()
chimescalc_setup_compiler_flags()

# Check whether libm is needed to provide the round() function in C
include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS("m" "round" "" NEEDS_LIB_M)
if(NEEDS_LIB_M)
    set(LIB_M "m")
endif()


####################################################################################################
# Static/shared library
####################################################################################################

set(c-sources
    "chimesFF/src/chimesFF.cpp"
    "chimesFF/api/chimescalc_C.cpp"
    "serial_interface/src/serial_chimes_interface.cpp"
    "serial_interface/api/chimescalc_serial_C.cpp")

set(c-includes
    "chimesFF/src/chimesFF.h"
    "chimesFF/api/chimescalc_C.h"
    "serial_interface/src/serial_chimes_interface.h"
    "serial_interface/api/chimescalc_serial_C.h")

add_library(ChimesCalc "${c-sources}")

set_target_properties(ChimesCalc
    PROPERTIES
    OUTPUT_NAME "chimescalc"
    PUBLIC_HEADER "${c-includes}")

target_include_directories(ChimesCalc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/api>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>)

target_compile_definitions(ChimesCalc PRIVATE "DEBUG=${DEBUG}")
target_compile_features(ChimesCalc PRIVATE cxx_std_11)

install(TARGETS ChimesCalc
    EXPORT ${PROJECT_NAME_LOWER}-targets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDEDIR}")


####################################################################################################
# Executables for C/C++
####################################################################################################

# Main executable installed as chimescalc
add_executable(ChimesCalcExe serial_interface/examples/cpp/main.cpp)
set_target_properties(ChimesCalcExe PROPERTIES OUTPUT_NAME "chimescalc")
target_link_libraries(ChimesCalcExe ChimesCalc)
target_compile_features   (ChimesCalcExe PUBLIC cxx_std_11)
target_compile_definitions(ChimesCalcExe PRIVATE "DEBUG=${DEBUG}")

install(TARGETS ChimesCalcExe DESTINATION ${CMAKE_INSTALL_BINDIR})

# Test executables
add_executable(chimescalc-test_direct-C chimesFF/examples/c/main.c)
target_link_libraries     (chimescalc-test_direct-C ChimesCalc ${LIB_M})
target_compile_features   (chimescalc-test_direct-C PRIVATE cxx_std_11)
target_compile_definitions(chimescalc-test_direct-C PRIVATE "DEBUG=${DEBUG}")

add_executable(chimescalc-test_serial-C serial_interface/examples/c/main.c)
target_link_libraries     (chimescalc-test_serial-C ChimesCalc)
target_compile_features   (chimescalc-test_serial-C PRIVATE cxx_std_11)
target_compile_definitions(chimescalc-test_serial-C PRIVATE "DEBUG=${DEBUG}")


####################################################################################################
# Dynamically loadable library (e.g for Python)
####################################################################################################

add_library(ChimesCalc_dynamic MODULE "${c-sources}")

set_target_properties(ChimesCalc_dynamic
    PROPERTIES OUTPUT_NAME "chimescalc_dl")

target_include_directories(ChimesCalc_dynamic PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/api>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>)

target_compile_definitions(ChimesCalc_dynamic PRIVATE "DEBUG=${DEBUG}")
target_compile_features(ChimesCalc_dynamic PRIVATE cxx_std_11)

install(TARGETS ChimesCalc_dynamic
    EXPORT ${PROJECT_NAME_LOWER}-targets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDEDIR}")


####################################################################################################
# Fortran wrappers (not contained in base library as it requires Fortran runtime library)
####################################################################################################

if(WITH_FORTRAN_API OR WITH_FORTRAN08_API)

    set(f90safe-fortran-sources
        "chimesFF/api/chimescalc_F.f90"
        "serial_interface/api/chimescalc_serial_F.f90")

    set(fortran-sources
        ${f90safe-fortran-sources})

    # Optional, as ancient Fortran compilers may have difficulties with it
    if(WITH_FORTRAN08_API)
        list(APPEND fortran-sources
            "serial_interface/api/chimescalc_serial_F08.f90")
    endif()

    add_library(ChimesCalc_Fortran "${fortran-sources}")
    set_target_properties(ChimesCalc_Fortran PROPERTIES OUTPUT_NAME "chimescalc_fortran")
    target_link_libraries(ChimesCalc_Fortran PUBLIC ChimesCalc)

    set(moddir "${CMAKE_CURRENT_BINARY_DIR}/modfiles")
    set_target_properties(ChimesCalc_Fortran PROPERTIES
        Fortran_MODULE_DIRECTORY "${moddir}")

    target_include_directories(ChimesCalc_Fortran PUBLIC
        $<BUILD_INTERFACE:${moddir}>
        $<INSTALL_INTERFACE:${INSTALL_MODULEDIR}>)

    install(TARGETS ChimesCalc_Fortran
        EXPORT ${PROJECT_NAME_LOWER}-targets
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

    install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MODULEDIR}")


    ################################################################################################
    # Executables for Fortran (90/08)
    ################################################################################################

    add_executable(chimescalc-test_serial-F serial_interface/examples/fortran/main.F90)
    target_link_libraries(chimescalc-test_serial-F ChimesCalc_Fortran)
    target_compile_definitions(chimescalc-test_serial-F PRIVATE "DEBUG=${DEBUG}")

    # Optional, as ancient Fortran compilers may have difficulties with it
    if(WITH_FORTRAN08_API)
        add_executable(chimescalc-test_serial-F08 serial_interface/examples/fortran08/main.F90)
        target_link_libraries(chimescalc-test_serial-F08 ChimesCalc_Fortran)
	      target_compile_definitions(chimescalc-test_serial-F08 PRIVATE "DEBUG=${DEBUG}")
    endif()

endif()

####################################################################################################
# Create CMake export file
####################################################################################################

include(CMakePackageConfigHelpers)

install(
    EXPORT ${PROJECT_NAME_LOWER}-targets
    FILE "${PROJECT_NAME_LOWER}-targets.cmake"
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/export/${PROJECT_NAME_LOWER}-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config-version.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME_LOWER}-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME_LOWER}")


####################################################################################################
# Add regression tests for the various APIs
####################################################################################################

enable_testing()

set(apis
    "cpp;${CMAKE_CURRENT_BINARY_DIR}/chimescalc"
    "c;${CMAKE_CURRENT_BINARY_DIR}/chimescalc-test_serial-C"
    "fortran;${CMAKE_CURRENT_BINARY_DIR}/chimescalc-test_serial-F"
    "fortran08;${CMAKE_CURRENT_BINARY_DIR}/chimescalc-test_serial-F08")

set(_testdir "${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/tests")
file(STRINGS ${_testdir}/test_list.dat _testcases)
foreach(_testcase IN LISTS _testcases)
    if("${_testcase}" STREQUAL "")
        continue()
    endif()
    list(GET _testcase 0 _paramfile)
    string(STRIP "${_paramfile}" _paramfile)
    list(GET _testcase 1 _geometry)
    string(STRIP "${_geometry}" _geometry)
    string(REPLACE "#" "" _geometry_escaped "${_geometry}")
    list(GET _testcase 2 _configopt)
    string(STRIP "${_configopt}" _configopt)
    list(GET _testcase 3 _testlabels)
    string(STRIP "${_testlabels}" _testlabels)

    set(_add_test False)
    foreach(_label IN LISTS TEST_LABELS)
        if("${_testlabels}" MATCHES "${_label}")
            set(_add_test True)
            break()
        endif()
    endforeach()
    if(NOT _add_test)
        continue()
    endif()

    set(_apis "${apis}")
    list(LENGTH _apis _apis_length )
    while(_apis_length GREATER 0)
        list(POP_FRONT _apis _api_abbrev)
        list(POP_FRONT _apis _api_executable)
        #message(STATUS "Test added: ${_api_abbrev}:${_paramfile}:${_geometry_escaped}")
        list(LENGTH _apis _apis_length)

        add_test(
            NAME "${_api_abbrev}/${_paramfile}:${_geometry_escaped}"
            COMMAND
                ${_testdir}/run_single_test.sh
                ${_api_executable}
                ${_paramfile}
                ${_geometry}
                ${_configopt}
                ${CMAKE_CURRENT_BINARY_DIR}/_test/${_api_abbrev}/${_paramfile}:${_geometry_escaped})
        if("${_api_abbrev}" STREQUAL "py")
            set_tests_properties("${_api_abbrev}/${_paramfile}:${_geometry_escaped}"
                PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
                ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}")
        endif()
    endwhile()
endforeach()
