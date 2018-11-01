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

include(geos-get-version-source)

get_version_source(
  HASH_LENGTH 7
  HASH  GEOS_GIT_HASH
  DIRTY GEOS_GIT_IS_DIRTY)

set(REV ${GEOS_GIT_HASH})
if (${GEOS_GIT_IS_DIRTY})
  set(REV "${REV} (dirty)")
endif()

message(STATUS "Setting GEOS_REVISION to ${REV}")

file(WRITE ${DST} "#define GEOS_REVISION \"${REV}\"
")
