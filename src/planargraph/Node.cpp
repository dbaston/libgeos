/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/planargraph/Node.h>
#include <geos/planargraph/DirectedEdge.h>

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

namespace geos {
namespace planargraph {

#if 0
/* static public */
/* UNUSED */
vector<Edge*>*
Node::getEdgesBetween(Node *node0, Node *node1)
{
	auto edges0 = DirectedEdge::toEdges(node0->getOutEdges().getEdges());
	auto edges1 = DirectedEdge::toEdges(node1->getOutEdges().getEdges());

	// Sort edge lists (needed for set_intersection below
	std::sort( edges0.begin(), edges0.end() );
	std::sort( edges1.begin(), edges1.end() );

	std::vector<Edge*>* commonEdges = new std::vector<Edge*>();

	// Intersect the two sets
	std::set_intersection(
		edges0.begin(), edges0.end(),
		edges1.begin(), edges1.end(),
		commonEdges->end()
		);

	return commonEdges;

}
#endif

std::ostream& operator<<(std::ostream& os, const Node& n) {
	os << "Node " << n.pt << " with degree " << n.getDegree();
	if ( n.isMarked() ) os << " Marked ";
	if ( n.isVisited() ) os << " Visited ";
	return os;
}

bool Node::hasDegree(size_t degree) const {
	return deStar.getDegree() == degree;
}


} // namespace planargraph
} // namespace geos
