cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE release)
project(1DPIC LANGUAGES C Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -O0 -extend-source 132")
else()
  set(CMAKE_Fortran_FLAGS "-g -O0 -ffree-line-length-none")
endif()

# 主程序在 code/PIC/Main_IFE_Test_2.f90；另需 MCC_jw（JW 的 Monte Carlo 碰撞接口 ModuleMCCInterface 等）
# 以及工程根目录下与主程序链接的若干子程序（勿与 MCC_jw/code/Main.F90 同时编译，后者为独立可执行）
file(GLOB_RECURSE SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/code/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/code/*.F90"
)
file(GLOB_RECURSE MCC_JW_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/MCC_jw/code/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/MCC_jw/code/*.F90"
)
list(REMOVE_ITEM MCC_JW_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/MCC_jw/code/Main.F90")

set(ROOT_FORTRAN_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/OUTPUT_velocity.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/Output_Energy.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/generate_Elementmap.f90"
)

add_executable(${PROJECT_NAME} ${SRCS} ${MCC_JW_SRCS} ${ROOT_FORTRAN_SRCS})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

