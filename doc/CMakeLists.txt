################################################################################
# Part of CMake configuration for GEOS
#
# Copyright (C) 2018 Daniel Baston <dbaston@gmail.com>
#
# This is free software; you can redistribute and/or modify it under
# the terms of the GNU Lesser General Public Licence as published
# by the Free Software Foundation.
# See the COPYING file for more information.
################################################################################

find_package(Doxygen)

if (DOXYGEN_FOUND)
  set(DOXYGEN_IN ${CMAKE_CURRENT_LIST_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  geos_configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT})

  add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
  )

  unset(DOXYGEN_IN)
  unset(DOXYGEN_OUT)
else()
  message(ERROR "Could not find doxygen.")
endif()
