cmake_minimum_required(VERSION 3.10)

############################################################
# set the project name
############################################################

project(ChIMES_Calculator
	VERSION 1.0
	DESCRIPTION "Utilities to compute ChIMES interactions"
	LANGUAGES CXX C Fortran)
    
set(CMAKE_C_COMPILER)
set(CMAKE_CXX_COMPILER)
set(CMAKE_F_COMPILER)

############################################################
# Executables/libraries for C++ via the serial wrapper
############################################################

# Define an executable target

add_executable(test-CPP
    ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/examples/cpp/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src/serial_chimes_interface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/chimesFF.cpp
)

# install(TARGETS myapp DESTINATION bin)    


# Define properties for the executable target

target_include_directories(test-CPP PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/examples/cpp)
target_include_directories(test-CPP PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src/)
target_include_directories(test-CPP PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/)

target_compile_features   (test-CPP PRIVATE cxx_std_11)
target_compile_definitions(test-CPP PRIVATE DEBUG=0)
target_compile_options    (test-CPP PRIVATE -O3)




############################################################
# Executables/libraries for C (direct)
############################################################

# Define an executable target

add_executable(test_wrapper-C-direct
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/examples/c/main.c
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/wrapper-C.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/chimesFF.cpp
)

# Define properties for the executable target

target_include_directories(test_wrapper-C-direct PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/examples/c/)
target_include_directories(test_wrapper-C-direct PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/)
target_include_directories(test_wrapper-C-direct PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/)

target_compile_features   (test_wrapper-C-direct PRIVATE cxx_std_11)
target_compile_options    (test_wrapper-C-direct PRIVATE -O3)



############################################################
# Executables/libraries for C via the serial wrapper
############################################################

# Define an executable target

add_executable(test_wrapper-C-serial_wrapper
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/examples/c/main.c
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/wrapper-C.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/chimesFF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src/serial_chimes_interface.cpp
)

# Define properties for the executable target

target_include_directories(test_wrapper-C-serial_wrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/examples/c/)
target_include_directories(test_wrapper-C-serial_wrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/)
target_include_directories(test_wrapper-C-serial_wrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/)
target_include_directories(test_wrapper-C-serial_wrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/src/)

target_compile_features   (test_wrapper-C-serial_wrapper PRIVATE cxx_std_11)
target_compile_options    (test_wrapper-C-serial_wrapper PRIVATE -O3)



############################################################
# Executables/libraries for Fortran (direct)
############################################################

# Define an executable target

add_executable(test_wrapper-F-direct
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/examples/fortran/main.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/wrapper-F.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/wrapper-C.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/chimesFF.cpp
)

# Define properties for the executable target

target_include_directories(test_wrapper-F-direct PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/examples/fortran/)
target_include_directories(test_wrapper-F-direct PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/)
target_include_directories(test_wrapper-F-direct PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/)

target_compile_features   (test_wrapper-F-direct PRIVATE cxx_std_11)
#target_compile_features   (test_wrapper-F-direct PRIVATE fortran_std_03)
target_compile_options    (test_wrapper-F-direct PRIVATE -O3)
target_compile_options    (test_wrapper-F-direct PRIVATE -fPIC)
target_compile_options    (test_wrapper-F-direct PRIVATE -std=f2003)



############################################################
# Executables/libraries for Fortran via the serial wrapper
############################################################

# Define an executable target

add_executable(test_wrapper-F-serial_wrapper
    ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/examples/fortran/main.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/wrapper-F.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/wrapper-C.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/chimesFF.cpp
)

# Define properties for the executable target

target_include_directories(test_wrapper-F-serial_wrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/serial_interface/examples/fortran/)
target_include_directories(test_wrapper-F-serial_wrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/api/)
target_include_directories(test_wrapper-F-serial_wrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/chimesFF/src/)

target_compile_features   (test_wrapper-F-serial_wrapper PRIVATE cxx_std_11)
target_compile_options    (test_wrapper-F-serial_wrapper PRIVATE -O3)
target_compile_options    (test_wrapper-F-serial_wrapper PRIVATE -fPIC)
target_compile_options    (test_wrapper-F-serial_wrapper PRIVATE -std=f2003)
target_compile_options    (test_wrapper-F-serial_wrapper PRIVATE -cpp)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_options(test_wrapper-F-serial_wrapper PRIVATE -lc++)
else()
    target_compile_options(test_wrapper-F-serial_wrapper PRIVATE -stdc++)    
endif()






















































#add_library(chimesFF-obj OBJECT ${SOURCES})


#add_executable(chimesFF $<TARGET_OBJECTS:chimesFF-obj>)

#target_compile_features   (chimesFF-obj PRIVATE cxx_std_11)








#cmake_minimum_required(VERSION 3.11)

#project(foo C)

#set(SOURCES foo.c)

#add_library(foo_obj OBJECT ${SOURCES})

#add_executable(foo $<TARGET_OBJECTS:foo_obj>)

#set(OBJS $<TARGET_OBJECTS:foo_obj>)
#message(STATUS "these are my objects: ${OBJS}") # generator expression printed here, not evaluated yet

#add_custom_target(print_foo_objs
#  COMMAND ${CMAKE_COMMAND} -E echo $<TARGET_OBJECTS:foo_obj>)

#set(OBJ_ROOT_DIR $<TARGET_PROPERTY:foo,BINARY_DIR>) 

#add_custom_target(copy_foo_objs
#  COMMAND ${CMAKE_COMMAND} -E make_directory ${OBJ_ROOT_DIR}/myobjects/
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_OBJECTS:foo_obj> ${OBJ_ROOT_DIR}/myobjects/)