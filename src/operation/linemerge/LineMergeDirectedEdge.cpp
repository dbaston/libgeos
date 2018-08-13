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
 * Last port: operation/linemerge/LineMergeDirectedEdge.java r378 (JTS-1.12)
 *
 **********************************************************************/

#include <geos/operation/linemerge/LineMergeDirectedEdge.h>
#include <geos/planargraph/DirectedEdge.h>
#include <geos/planargraph/Node.h>
#include <geos/planargraph/detail.hpp>

#include <cassert>

//using namespace geos::planargraph;
using namespace geos::geom;
using geos::planargraph::Node;
using geos::planargraph::DirectedEdge;

namespace geos {
namespace operation { // geos.operation
namespace linemerge { // geos.operation.linemerge

LineMergeDirectedEdge::LineMergeDirectedEdge(
		NodePtr newFrom,
		NodePtr newTo,
		const Coordinate& newDirectionPt,
		bool nEdgeDirection)
	:
	planargraph::DirectedEdge(newFrom, newTo,
			newDirectionPt, nEdgeDirection)
{}

/**
 * Returns the directed edge that starts at this directed edge's end point,
 * or null if there are zero or multiple directed edges starting there.
 * @return
 */
DirectedEdge::DirectedEdgePtr
LineMergeDirectedEdge::getNext() const
{
	if (getToNode()->getDegree() != 2) {
		return nullptr;
	}
	if (getToNode()->getOutEdges().getEdges()[0] == getSym()) {
	  return getToNode()->getOutEdges().getEdges()[1];
	}
	assert(getToNode()->getOutEdges().getEdges()[1] == getSym());

	auto nextedge = getToNode()->getOutEdges().getEdges()[0];
	assert(nextedge);

	return nextedge;
}

} // namespace geos.operation.linemerge
} // namespace geos.operation
} //namespace geos
