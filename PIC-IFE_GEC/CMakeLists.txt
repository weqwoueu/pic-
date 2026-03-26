cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE release)
project(1DPIC LANGUAGES C Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -O0 -extend-source 132")
else()
  set(CMAKE_Fortran_FLAGS "-g -O0 -ffree-line-length-none")
endif()

# 源码在 code/ 下（与仓库布局一致）；登录超算后请用:
#   cmake -S . -B build -DCMAKE_Fortran_COMPILER=ifort
file(GLOB_RECURSE SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/code/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/code/*.F90"
)
add_executable(${PROJECT_NAME} ${SRCS})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

