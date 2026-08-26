cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE release)
project(1DPIC LANGUAGES C Fortran)

set(PIC_CPU_TARGET "host" CACHE STRING "CPU target: host, avx2, or portable")
set_property(CACHE PIC_CPU_TARGET PROPERTY STRINGS host avx2 portable)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
  if(PIC_CPU_TARGET STREQUAL "host")
    set(PIC_ARCH_FLAGS "-xHost")
  elseif(PIC_CPU_TARGET STREQUAL "avx2")
    set(PIC_ARCH_FLAGS "-xCORE-AVX2")
  elseif(PIC_CPU_TARGET STREQUAL "portable")
    set(PIC_ARCH_FLAGS "")
  else()
    message(FATAL_ERROR "PIC_CPU_TARGET must be host, avx2, or portable")
  endif()
  set(CMAKE_Fortran_FLAGS "-O3 ${PIC_ARCH_FLAGS} -ipo -extend-source 132")
else()
  if(PIC_CPU_TARGET STREQUAL "host")
    set(PIC_ARCH_FLAGS "-march=native")
  elseif(PIC_CPU_TARGET STREQUAL "avx2")
    set(PIC_ARCH_FLAGS "-march=haswell")
  elseif(PIC_CPU_TARGET STREQUAL "portable")
    set(PIC_ARCH_FLAGS "")
  else()
    message(FATAL_ERROR "PIC_CPU_TARGET must be host, avx2, or portable")
  endif()
  set(CMAKE_Fortran_FLAGS "-O3 ${PIC_ARCH_FLAGS} -ffree-line-length-none")
endif()

message(STATUS "PIC CPU target: ${PIC_CPU_TARGET} (${PIC_ARCH_FLAGS})")

# 主程序在 code/PIC/Main_IFE_Test_2.f90；另需 MCC_jw（JW 的 Monte Carlo 碰撞接口 ModuleMCCInterface 等）
# 以及工程根目录下与主程序链接的若干子程序（勿与 MCC_jw/code/Main.F90 同时编译，后者为独立可执行）
file(GLOB_RECURSE SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/code/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/code/*.F90"
)
list(REMOVE_ITEM SRCS "${CMAKE_CURRENT_SOURCE_DIR}/code/PIC/DRandom.f90")
file(GLOB_RECURSE LEGACY_MCC_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/code/MCC/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/code/MCC/*.F90"
)
list(REMOVE_ITEM SRCS ${LEGACY_MCC_SRCS})
file(GLOB_RECURSE MCC_JW_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/MCC_jw/code/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/MCC_jw/code/*.F90"
)
list(REMOVE_ITEM MCC_JW_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/MCC_jw/code/Main.F90")

set(ROOT_FORTRAN_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/DiagCounters.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/PerformanceTimers.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/OUTPUT_velocity.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/Output_Energy.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/GlobalDiagnostics.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/generate_Elementmap.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/InjectPositioning.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/ParticlePositioning.f90"
)

add_executable(${PROJECT_NAME} ${SRCS} ${MCC_JW_SRCS} ${ROOT_FORTRAN_SRCS})

find_package(OpenMP REQUIRED COMPONENTS Fortran)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

