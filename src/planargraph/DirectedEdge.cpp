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

#include <geos/planargraph/DirectedEdge.h>
#include <geos/planargraph/Node.h>
#include <geos/geomgraph/Quadrant.h>
#include <geos/algorithm/CGAlgorithms.h>

#include <cmath>
#include <sstream>
#include <vector>
#include <typeinfo>

using namespace std;
using namespace geos::geom;

namespace geos {
namespace planargraph {

/*public deprecated*/
void
DirectedEdge::toEdges(vector<DirectedEdge*>& dirEdges, vector<Edge*>& edges)
{
	for (auto e : dirEdges)
	{
		edges.push_back(e->m_parentEdge);
	}
}

/*public*/
vector<Edge*>
DirectedEdge::toEdges(vector<DirectedEdge*>& dirEdges)
{
	vector<Edge*> edges;
	for (auto e : dirEdges) edges.push_back(e->m_parentEdge);
	return edges;
}

/*public*/
DirectedEdge::DirectedEdge(Node* newFrom, Node* newTo,
	const Coordinate &directionPt, bool newEdgeDirection) :
	m_from(newFrom),
	m_to(newTo),
	m_p0(m_from->getCoordinate()),
	m_p1(directionPt),
	m_edgeDirection(newEdgeDirection)
{
	double dx = m_p1.x - m_p0.x;
	double dy = m_p1.y - m_p0.y;
	m_quadrant = geomgraph::Quadrant::quadrant(dx, dy);
	m_angle = atan2(dy, dx);
	//Assert.isTrue(! (dx == 0 && dy == 0), "EdgeEnd with identical endpoints found");
}

/*public*/
Edge*
DirectedEdge::parentEdge() const
{
	return m_parentEdge;
}

// [[deprecated]]
Edge*
DirectedEdge::getEdge() const
{
	return m_parentEdge;
}

/*public*/
void
DirectedEdge::set_parentEdge(Edge* newParentEdge)
{
	m_parentEdge = newParentEdge;
}

// [[deprecated]]
void
DirectedEdge::setEdge(Edge* newParentEdge)
{
	m_parentEdge=newParentEdge;
}

/*public*/
int
DirectedEdge::getQuadrant() const
{
	return m_quadrant;
}

/*public*/
const Coordinate&
DirectedEdge::getDirectionPt() const
{
	return m_p1;
}

/*public*/
bool
DirectedEdge::getEdgeDirection() const
{
	return m_edgeDirection;
}

/*public*/
Node*
DirectedEdge::getFromNode() const
{
	return m_from;
}

/*public*/
Node*
DirectedEdge::getToNode() const
{
	return m_to;
}

/*public*/
Coordinate&
DirectedEdge::getCoordinate() const
{
	return m_from->getCoordinate();
}

/*public*/
double
DirectedEdge::getAngle() const
{
	return m_angle;
}

/*public*/
DirectedEdge*
DirectedEdge::getSym() const
{
	return m_sym;
}

/*
 * Sets this DirectedEdge's symmetric DirectedEdge,
 * which runs in the opposite direction.
 */
void
DirectedEdge::setSym(DirectedEdge *newSym)
{
	m_sym = newSym;
}

/*
 * Resests the sym ptr to nullptr
 */
void
DirectedEdge::resetSym()
{
	m_sym = nullptr;
}

/*public*/
int
DirectedEdge::compareTo(const DirectedEdge de) const
{
	return compareDirection(de);
}

/*public*/
int
DirectedEdge::compareDirection(const DirectedEdge e) const
{
// if the rays are in different quadrants, determining the ordering is trivial
	if (m_quadrant > e.m_quadrant) return 1;
	if (m_quadrant < e.m_quadrant) return -1;
	// vectors are in the same quadrant - check relative orientation of direction vectors
	// this is > e if it is CCW of e
	return algorithm::CGAlgorithms::computeOrientation(e.m_p0, e.m_p1, m_p1);
}

/*public*/
string
DirectedEdge::print() const
{
	ostringstream s;
  s << *this;
	return s.str();
}

std::ostream&
operator << (std::ostream& s, const DirectedEdge& de)
{
  s << typeid(de).name() << ": " << de.m_p0 << " - " << de.m_p1;
  s << " " << de.m_quadrant << ":" << de.m_angle;
  return s;
}

} // namespace planargraph
} // namespace geos

