cmake_minimum_required(VERSION 3.4.0)

#set(CMAKE_C_COMPILER "icc")
#set(CMAKE_Fortran_COMPILER "ifort")
set(CMAKE_Fortran_COMPILER 
        /opt/oneapi/compiler/2022.1.0/linux/bin/intel64/ifort
        CACHE FILEPATH "Fortran compiler")

set(CMAKE_Fortran_FLAGS "-g -O0 -ffree-line-length-none")

set(CMAKE_BUILD_TYPE release)
project(1DPIC LANGUAGES C Fortran)
#file(REMOVE_RECURSE  ${PROJECT_SOURCE_DIR}/run)

file(GLOB_RECURSE SRCS ./mc/*.f90 ./mc/*.F90 ./base/*.f90  ./base/*.F90  ./Interface_IFE/*.f90  ./Interface_IFE/*.F90
     ./AxialSymmetric/*.f90
    ./Data/*.f90 ./Debug_zwz/*.f90 ./DSMC/*.f90 ./Function/*.f90 ./IFE-error/*.f90 
    ./IFE-Rectangular/*.f90 ./IFE-Solver/*.f90 ./IMPIC/*.f90 ./In-Output/*.f90 
     ./PIC/*.f90 ./Setup/*.f90 ./PPR/*.f90  ./diagnostics/*.f90 ./SIDG/*.f90 )
add_executable(${PROJECT_NAME} ${SRCS})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

