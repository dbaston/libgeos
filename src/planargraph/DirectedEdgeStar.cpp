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

#include <geos/planargraph/DirectedEdgeStar.h>
#include <geos/planargraph/DirectedEdge.h>
#include <geos/planargraph/detail.hpp>

#include <vector>
#include <algorithm>

using namespace geos::geom;

namespace geos {
namespace planargraph {


/*
 * Adds a new member to this DirectedEdgeStar.
 */
void
DirectedEdgeStar::add(DirectedEdgePtr de)
{
	outEdges.push_back(de);
	sorted = false;
}

/*
 * Drops a member of this DirectedEdgeStar.
 */
void
DirectedEdgeStar::remove(DirectedEdgePtr de)
{
	find_and_erase(de, outEdges);
}

DirectedEdgeStar::DirectedEdges::iterator
DirectedEdgeStar::begin() {
	sortEdges();
	return outEdges.begin();
}

DirectedEdgeStar::DirectedEdges::iterator
DirectedEdgeStar::end() {
	sortEdges();
	return outEdges.end();
}

DirectedEdgeStar::DirectedEdges::const_iterator
DirectedEdgeStar::begin() const {
	sortEdges();
	return outEdges.begin();
}

std::vector<DirectedEdge*>::const_iterator
DirectedEdgeStar::end() const {
	sortEdges();
	return outEdges.end();
}

/*
 * Returns the coordinate for the node at wich this star is based
 */
Coordinate&
DirectedEdgeStar::getCoordinate() const
{
	sortEdges();
	if (outEdges.empty())
		return Coordinate::getNull();
	return outEdges[0]->getCoordinate();
}

/*
 * Returns the DirectedEdges, in ascending order by angle with
 * the positive x-axis.
 */
std::vector<DirectedEdge*>&
DirectedEdgeStar::getEdges()
{
	sortEdges();
	return outEdges;
}

bool
pdeLessThan(DirectedEdge *first, DirectedEdge * second)
{
	return first->compareTo(*second) < 0;
}

/*private*/
void
DirectedEdgeStar::sortEdges() const
{
	if (!sorted) {
		std::sort(outEdges.begin(), outEdges.end(), pdeLessThan);
		sorted = true;
	}
}

#ifdef GEOS_USEDEPRECATED
/*
 * Returns the zero-based index of the given Edge, after sorting in
 * ascending order by angle with the positive x-axis.
 */
int
DirectedEdgeStar::getIndex(const Edge *edge)
{
	sortEdges();
	for (unsigned int i = 0; i<outEdges.size(); ++i)
	{
		DirectedEdge *de =outEdges[i];
		if (de->parentEdge() == edge)
		return i;
	}
	return -1;
}

/*
 * Returns the zero-based index of the given DirectedEdge, after sorting
 * in ascending order by angle with the positive x-axis.
 */
int
DirectedEdgeStar::getIndex(const DirectedEdge *dirEdge)
{
	sortEdges();
	for (unsigned int i = 0; i <outEdges.size(); ++i)
	{
		DirectedEdge *de =outEdges[i];
		if (de == dirEdge)
		return i;
	}
	return -1;
}

/*
 * Returns the remainder when i is divided by the number of edges in this
 * DirectedEdgeStar.
 */
int
DirectedEdgeStar::getIndex(int i) const
{
	int modi = i % (int)outEdges.size();
	//I don't think modi can be 0 (assuming i is positive) [Jon Aquino 10/28/2003]
	if (modi < 0) modi += (int)outEdges.size();
	return modi;
}
#endif

std::vector<geos::planargraph::DirectedEdge*>::iterator
DirectedEdgeStar::findEdge(const DirectedEdge *dirEdge) const
{
	sortEdges();
	return std::find(outEdges.begin(), outEdges.end(), dirEdge);
}

/*
 * Returns the DirectedEdge on the left-hand side of the given
 * DirectedEdge (which must be a member of this DirectedEdgeStar).
 */
DirectedEdge*
DirectedEdgeStar::getNextEdge(DirectedEdge *dirEdge)
{
	auto it = findEdge(dirEdge);
	return *(it + 1);
}

} // namespace planargraph
} // namespace geos

