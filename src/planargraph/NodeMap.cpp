/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 * Copyright (C) 2005 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/planargraph/NodeMap.h>

/* for std::pair */
#include <utility>

#include <geos/planargraph/Node.h>

namespace geos {
namespace planargraph {

NodeMap::NodeContainer&
NodeMap::getNodeMap()
{
	return nodeMap;
}

/**
 * Adds a node to the map, replacing any that is already at that location.
 * @return the added node
 */
NodeMap::NodePtr
NodeMap::add(NodePtr n)
{
	nodeMap.insert(std::pair<geom::Coordinate, NodePtr> (n->getCoordinate(), n));
	return n;
}

/**
 * Removes the Node at the given location, and returns it
 * (or null if no Node was there).
 */
NodeMap::NodePtr
NodeMap::remove(geom::Coordinate& pt)
{
	auto n = find(pt);
	nodeMap.erase(pt);
	return n;
}

/* public */
NodeMap::NodeVector
NodeMap::getNodes() const
{
	NodeVector values;

	for (const auto e : nodeMap) {
		values.push_back(e.second);
	}

	return values;
}

// [[deprecated]]
void
NodeMap::getNodes(NodeVector& values) const
{
	values = getNodes();
}

/**
 * Returns the Node ptr at the given location, or null if no Node was there.
 */
NodeMap::NodePtr
NodeMap::find(const geom::Coordinate& coord) const
{
	auto found = nodeMap.find(coord);
	return (found == nodeMap.end()) ?
		nullptr :
		found->second;
}

}  // namespace planargraph
}  // namespace geos

