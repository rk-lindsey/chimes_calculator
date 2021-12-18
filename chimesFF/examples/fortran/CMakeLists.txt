cmake_minimum_required(VERSION 3.10)

############################################################
# General settings
############################################################

project(TestChimes
    DESCRIPTION "Tests the Fortran interface of the ChIMES calculator"
  	LANGUAGES Fortran CXX)

find_package(ChimesCalc REQUIRED)
add_executable(test_chimescalc "main.F90")
target_link_libraries(test_chimescalc ChimesCalc::ChimesCalc_Fortran)
