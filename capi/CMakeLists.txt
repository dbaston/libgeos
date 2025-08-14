################################################################################
# Part of CMake configuration for GEOS
#
# Copyright (C) 2018 Mateusz Loskot <mateusz@loskot.net>
#
# This is free software; you can redistribute and/or modify it under
# the terms of the GNU Lesser General Public Licence as published
# by the Free Software Foundation.
# See the COPYING file for more information.
################################################################################

target_sources(geos_c PUBLIC
                      FILE_SET HEADERS 
                      BASE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                      FILES ${CMAKE_CURRENT_BINARY_DIR}/geos_c.h)
target_sources(geos_c PRIVATE geos_c.cpp geos_ts_c.cpp)

target_include_directories(geos_c
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/geos>)

target_include_directories(geos_c
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/geos>)

# Copy these over so they match the @VARIABLES@ used by autoconf
# in geos_c.h.in
set(VERSION ${GEOS_VERSION})
set(VERSION_MAJOR ${GEOS_VERSION_MAJOR})
set(VERSION_MINOR ${GEOS_VERSION_MINOR})
set(VERSION_PATCH ${GEOS_VERSION_PATCH})

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/geos_c.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/geos_c.h
  @ONLY)

unset(VERSION)
unset(VERSION_MAJOR)
unset(VERSION_MINOR)
unset(VERSION_PATCH)
