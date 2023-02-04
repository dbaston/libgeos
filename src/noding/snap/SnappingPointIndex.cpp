/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2020 Paul Ramsey <pramsey@cleverelephant.ca>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/noding/snap/SnappingPointIndex.h>

using namespace geos::geom;

namespace geos {
namespace noding { // geos.noding
namespace snap { // geos.noding.snap

SnappingPointIndex::SnappingPointIndex(double p_snapTolerance) :
    // snapTolerance(p_snapTolerance),
    snapPointIndex(new index::kdtree::KdTree(p_snapTolerance)) {}


const CoordinateXY&
SnappingPointIndex::snap(const CoordinateXY& p)
{
    /**
    * Inserting the coordinate snaps it to any existing
    * one within tolerance, or adds it if not.
    */
    const index::kdtree::KdNode* node = snapPointIndex->insert(p);
    return node->getCoordinate();
}



} // namespace geos.noding.snap
} // namespace geos.noding
} // namespace geos
