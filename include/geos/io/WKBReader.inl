/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: io/WKBReader.java rev. 1.1 (JTS-1.7)
 *
 **********************************************************************/

#ifndef GEOS_IO_WKBREADER_INL
#define GEOS_IO_WKBREADER_INL

#include <geos/io/WKBReader.h>
#include <geos/geom/GeometryFactory.h>

#if GEOS_DEBUG
# include <iostream>
#endif

namespace geos {
namespace io {

INLINE geom::Coordinate
WKBReader::readCoordinate()
{
    const geom::PrecisionModel& pm = *factory.getPrecisionModel();
    geom::Coordinate ret{dis.readDouble(), dis.readDouble()};

    if (inputDimension == 3) {
        ret.z = dis.readDouble();
    }

    pm.makePrecise(ret);

    return ret;
}

} // namespace io
} // namespace geos

#endif // #ifndef GEOS_IO_WKTREADER_INL
