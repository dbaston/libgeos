/**********************************************************************
*
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2025 ISciences, LLC
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/noding/ArcString.h>

namespace geos::noding {

std::unique_ptr<geom::CoordinateSequence>
ArcString::getCoordinates() const
{
 // FIXME Z/M?

 if (m_arcs.empty()) {
  return std::make_unique<geom::CoordinateSequence>();
 }

 const std::size_t seqSize = m_arcs.size() * 2 + 1;

 // fixme
 constexpr bool hasZ = false;
 constexpr bool hasM = false;
 constexpr bool initialize = false;

 auto seq = std::make_unique<geom::CoordinateSequence>(seqSize, hasZ, hasM, initialize);
 std::size_t i = 0;
 for (const auto& arc : m_arcs) {
  if (i == 0) {
   seq->setAt(arc.p0, i++);
  }
  seq->setAt(arc.p1, i++);
  seq->setAt(arc.p2, i++);
 }

 return seq;
}

}