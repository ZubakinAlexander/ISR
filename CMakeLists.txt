cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(XSR)

# Options
option(BOOST_PREFIX "The path to the Boost installation")
option(GSL_PREFIX "The path to the GSL installation")
option(ROOT_PREFIX "The path to the CERN ROOT installation")

# ROOT
if (NOT ROOT_PREFIX)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
else()
  list(APPEND CMAKE_PREFIX_PATH ${ROOT_PREFIX})
endif()
find_package(ROOT REQUIRED COMPONENTS RIO Net Core)
if (ROOT_FOUND)
  message("[v] CERN ROOT")
endif()
include(${ROOT_USE_FILE})

# GSL
if (NOT GSL_PREFIX)
  execute_process(COMMAND gsl-config --prefix OUTPUT_VARIABLE GSL_PREFIX)
  string(REGEX REPLACE "\n$" "" GSL_PREFIX "${GSL_PREFIX}")
  message("GSL_PREFIX: " ${GSL_PREFIX})
endif()
set(GSL_INCLUDE_DIR ${GSL_PREFIX}/include)
set(GSL_LIBRARY ${GSL_PREFIX}/lib/libgsl.so)
set(GSL_CBLAS_LIBRARY ${GSL_PREFIX}/lib/libgslcblas.so)
find_package(GSL REQUIRED)
if(${GSL_FOUND})
  message("[v] GSL")
endif()

# BOOST

if (BOOST_PREFIX)
  set (BOOST_ROOT BOOST_PREFIX)
  set (BOOST_INCLUDEDIR ${BOOST_PREFIX}/include)
  set (BOOST_LIBRARYDIR ${BOOST_PREFIX}/lib)
  set (Boost_NO_SYSTEM_PATHS ON)
  set (Boost_NO_BOOST_CMAKE ON)
endif()
find_package(Boost COMPONENTS program_options REQUIRED)
if(${Boost_FOUND})
  message("[v] Boost")
endif()

# Define variables

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB SOURCES_NOMAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM SOURCES_NOMAIN main.cpp)
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.hpp)
set(PUBLIC_HEADERS ${HEADERS})

# Create library

add_library(ISR SHARED ${SOURCES_NOMAIN})
target_link_libraries(ISR ${ROOT_LIBRARIES}
  GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_include_directories(ISR
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
  $<INSTALL_INTERFACE:inc>
  PRIVATE ${GSL_INCLUDE_DIR}
  PRIVATE ${Boost_INCLUDE_DIRS})

set_target_properties(ISR
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")

install(TARGETS ISR EXPORT ISRTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION inc)

install(EXPORT ISRTargets
  FILE XSRConfig.cmake
  DESTINATION cmake)
