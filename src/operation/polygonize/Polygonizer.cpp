/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2010 Sandro Santilli <strk@kbt.io>
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
 * Last port: operation/polygonize/Polygonizer.java rev. 1.6 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/operation/polygonize/Polygonizer.h>
#include <geos/operation/polygonize/PolygonizeGraph.h>
#include <geos/operation/polygonize/EdgeRing.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Polygon.h>
#include <geos/util/Interrupt.h>
// std
#include <vector>

#ifdef _MSC_VER
#pragma warning(disable:4355)
#endif

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

using namespace geos::geom;

namespace geos {
namespace operation { // geos.operation
namespace polygonize { // geos.operation.polygonize

/* private */
Polygonizer::LineStringAdder::LineStringAdder(Polygonizer *p):
	pol(p)
{
}

/* private */
void
Polygonizer::LineStringAdder::filter_ro(const Geometry *g)
{
	const LineString *ls = dynamic_cast<const LineString *>(g);
	if ( ls ) pol->add(ls);
}


/*
 * Create a polygonizer with the same GeometryFactory
 * as the input Geometry
 */
Polygonizer::Polygonizer():
	lineStringAdder(this),
	graph(),
	dangles(),
	cutEdges(),
	invalidRingLines(),
	m_holeList(),
	shellList(),
	polyList()
{
}

/* C++ interface constructors */
Polygonizer::Polygonizer(const std::vector<geom::Geometry*> geomList) :
	Polygonizer()
{
	for (const auto g : (geomList)) add(g);
}

Polygonizer::Polygonizer(const std::vector<const geom::Geometry*> geomList) :
	Polygonizer()
{
	for (const auto g : (geomList)) add(g);
}

Polygonizer::Polygonizer(const geom::Geometry *g) :
	Polygonizer()
{
	add(g);
}




Polygonizer::~Polygonizer()
{
	clear();
}

void Polygonizer::clear() {
#if 0
	for (auto &p : polyList) delete p;
	polyList.clear();

	for (auto &r : invalidRingLines) delete r;
	invalidRingLines.clear();
#endif
}

void
Polygonizer::add(const std::vector<Geometry*> geomList)
{
	clear();
	for (const auto g : (geomList)) add(g);
}

void
Polygonizer::add(const std::vector<const Geometry*> geomList)
{
	clear();
	for (const auto g : (geomList)) add(g);
}

/* public
 * Add a collection of geometries to be polygonized.
 * May be called multiple times.
 * Any dimension of Geometry may be added;
 * the constituent linework will be extracted and used
 *
 * @param geomList a list of {@link Geometry}s with linework to be polygonized
 */
void
Polygonizer::add(const std::vector<Geometry*> *geomList)
{
	clear();
	for (const auto g : (*geomList)) add(g);
}

/* public
 * Add a collection of geometries to be polygonized.
 * May be called multiple times.
 * Any dimension of Geometry may be added;
 * the constituent linework will be extracted and used
 *
 * @param geomList a list of {@link Geometry}s with linework to be polygonized
 */
void
Polygonizer::add(const std::vector<const Geometry*> *geomList)
{
	clear();
	for (auto &g : (*geomList)) add(g);
}

#ifdef GEOS_USEDEPRECATED
/*
 * Add a geometry to the linework to be polygonized.
 * May be called multiple times.
 * Any dimension of Geometry may be added;
 * the constituent linework will be extracted and used
 *
 * @param g a Geometry with linework to be polygonized
 */
void
Polygonizer::add(Geometry *g)
{
	g->apply_ro(&lineStringAdder);
}
#endif

/* public
 * Add a geometry to the linework to be polygonized.
 * May be called multiple times.
 * Any dimension of Geometry may be added;
 * the constituent linework will be extracted and used
 *
 * @param g a Geometry with linework to be polygonized
 */
void
Polygonizer::add(const Geometry *g)
{
	clear();
	g->apply_ro(&lineStringAdder);
}


/* private
 * Add a linestring to the graph of polygon edges.
 *
 * @param line the LineString to add
 */
void
Polygonizer::add(const LineString *line)
{
	// create a new graph using the factory from the input Geometry
	if (graph.empty())
		graph = PolygonizeGraph(line->getFactory());
	graph.addEdge(line);
}

/*
 * Gets the list of polygons formed by the polygonization.
 * @return a collection of Polygons
 */
std::vector<Polygon*>
Polygonizer::getPolygons()
{
	polygonize();
	return polyList;
}

/* public */
const std::vector<const LineString*>&
Polygonizer::getDangles()
{
	polygonize();
	return dangles;
}

/* public */
const std::vector<const LineString*>&
Polygonizer::getCutEdges()
{
	polygonize();
	return cutEdges;
}

/* public */
const std::vector<LineString*>&
Polygonizer::getInvalidRingLines()
{
	polygonize();
	return invalidRingLines;
}

/* private */
void
Polygonizer::polygonize()
{
	// check if already computed
	if (!polyList.empty()) return;

	polyList.clear();

	// if no geometries were supplied it's possible graph could be null
	if (graph.empty()) return;

	dangles = graph.deleteDangles();
	cutEdges = graph.deleteCutEdges();

	auto edgeRingList = graph.getEdgeRings();

#if GEOS_DEBUG
	cerr<<"Polygonizer::polygonize(): "<<edgeRingList.size()<<" edgeRings in graph"<<endl;
#endif

	auto validEdgeRingList = findValidRings(edgeRingList);

#if GEOS_DEBUG
	cerr<<"                           "<<validEdgeRingList.size()<<" valid"<<endl;
	cerr<<"                           "<<invalidRingLines.size()<<" invalid"<<endl;
#endif

	findShellsAndHoles(validEdgeRingList);

#if GEOS_DEBUG
	cerr<<"                           "<<holeList.size()<<" holes"<<endl;
	cerr<<"                           "<<shellList.size()<<" shells"<<endl;
#endif

	assignHolesToShells(m_holeList, shellList);

	for (const auto &er : shellList)
	{
		polyList.push_back(er->getPolygon());
	}
}

/* private */
std::vector<Polygonizer::EdgeRingPtr>
Polygonizer::findValidRings(std::vector<EdgeRingPtr>& edgeRingList) const
{
  std::vector<EdgeRingPtr> validEdgeRingList;
	for (auto &er : edgeRingList)
	{
		if (er->isValid())
		{
			validEdgeRingList.push_back(std::move(er));
		}
		else
		{
      /* TODO check how the linestirngs are build */
			invalidRingLines.push_back(er.release()->getLineString());
		}
		GEOS_CHECK_FOR_INTERRUPTS();
	}
  edgeRingList.clear();
	return validEdgeRingList;
}

/* private */
void
Polygonizer::findShellsAndHoles(std::vector<EdgeRingPtr>& edgeRingList)
{
	m_holeList.clear();
	shellList.clear();
	for (auto &er : edgeRingList) {
		if (er->isHole())
			m_holeList.push_back(std::move(er));
		else
			shellList.push_back(std::move(er));

		GEOS_CHECK_FOR_INTERRUPTS();
	}
  edgeRingList.clear();
}

/* private */
void
Polygonizer::assignHolesToShells(
    std::vector<EdgeRingPtr>& holeList,
    std::vector<EdgeRingPtr>& shellList)
{
  std::vector<EdgeRingPtr> dangleHoles;
	for (auto &holeER : holeList) {
    auto shell = holeER->findEdgeRingContaining(shellList);

    if (shell != shellList.end()) {
      (*shell)->addHole(holeER.release()->getRingOwnership());
      assert(!holeER);
    } else {
      dangleHoles.push_back(std::move(holeER));
    }
    GEOS_CHECK_FOR_INTERRUPTS();
  }

  holeList.clear();
  holeList = std::move(dangleHoles);
}

#if 0
/* private */
void
Polygonizer::assignHoleToShell(
    EdgeRingPtr holeER,
		const vector<EdgeRingPtr>& shellList)
{
	auto shell = holeER->findEdgeRingContaining(shellList);

	if (shell) {
		shell->addHole(holeER.release()->getRingOwnership());
  }
}
#endif

} // namespace geos.operation.polygonize
} // namespace geos.operation
} // namespace geos

