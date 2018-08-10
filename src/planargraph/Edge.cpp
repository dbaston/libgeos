/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/planargraph/Edge.h>
#include <geos/planargraph/DirectedEdge.h>
#include <geos/planargraph/Node.h>

#include <iostream>

namespace geos {
namespace planargraph {

/*
 * Initializes this Edge's two DirectedEdges, and for each DirectedEdge:
 *  sets the Edge,
 *  sets the symmetric DirectedEdge, and
 *  adds this Edge to its from-Node.
 */
void
Edge::setDirectedEdges(DirectedEdgePtr de0, DirectedEdgePtr de1)
{
	dirEdge.push_back(de0);
	dirEdge.push_back(de1);
	de0->set_parentEdge(this);
	de1->set_parentEdge(this);
	de0->setSym(de1);
	de1->setSym(de0);
	de0->getFromNode()->addOutEdge(de0);
	de1->getFromNode()->addOutEdge(de1);
}

/**
* Returns one of the DirectedEdges associated with this Edge.
* @param i 0 or 1
*/
Edge::DirectedEdgePtr
Edge::getDirEdge(int i)
{
	return dirEdge[i];
}

/*
 * Returns the DirectedEdge that starts from the given node, or null if the
 * node is not one of the two nodes associated with this Edge.
 */
Edge::DirectedEdgePtr
Edge::getDirEdge(NodePtr fromNode)
{
	if (dirEdge[0]->getFromNode()==fromNode) return dirEdge[0];
	if (dirEdge[1]->getFromNode()==fromNode) return dirEdge[1];
	// node not found
	// possibly should throw an exception here?
	return nullptr;
}

/**
* If <code>node</code> is one of the two nodes associated with this Edge,
* returns the other node; otherwise returns null.
*/
Edge::NodePtr
Edge::getOppositeNode(NodePtr node)
{
	if (dirEdge[0]->getFromNode()==node) return dirEdge[0]->getToNode();
	if (dirEdge[1]->getFromNode()==node) return dirEdge[1]->getToNode();
	// node not found
	// possibly should throw an exception here?
	return nullptr;
}

std::ostream& operator<<(std::ostream& os, const Edge& n) {
	os << "Edge ";
	if ( n.isMarked() ) os << " Marked ";
	if ( n.isVisited() ) os << " Visited ";
	return os;
}

} // namespace planargraph
} // namespace geos

