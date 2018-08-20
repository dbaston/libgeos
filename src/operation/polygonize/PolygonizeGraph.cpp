/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2018 ~    Vicky Vergara
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
 * Modernized to c++11: on 2018
 * Last port: operation/polygonize/PolygonizeGraph.java rev. 6/138 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/operation/polygonize/PolygonizeGraph.h>
#include <geos/operation/polygonize/PolygonizeDirectedEdge.h>
#include <geos/operation/polygonize/PolygonizeEdge.h>
#include <geos/operation/polygonize/EdgeRing.h>
#include <geos/planargraph/Node.h>
#include <geos/planargraph/DirectedEdgeStar.h>
#include <geos/planargraph/DirectedEdge.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/LineString.h>
#include <geos/planargraph/detail.hpp>

#include <cassert>
#include <vector>
#include <set>

using geos::planargraph::Node;
using geos::planargraph::Edge;
using geos::planargraph::DirectedEdge;
using geos::planargraph::DirectedEdgeStar;
using geos::geom::LineString;
using geos::geom::Coordinate;
using geos::geom::CoordinateSequence;
using geos::geom::GeometryFactory;

namespace geos {
namespace operation {  // geos.operation
namespace polygonize {  // geos.operation.polygonize


/*
 * Create a new polygonization graph.
 */
PolygonizeGraph::PolygonizeGraph(const GeometryFactory *newFactory)
	: m_factory(newFactory) {
}

/*
 * Destroy a PolygonizeGraph
 */
PolygonizeGraph::~PolygonizeGraph() {
	//for (auto e : m_newEdgeRings) delete e;
	for (auto c : m_newCoords) delete c;
}

/*
 * Add a LineString forming an edge of the polygon graph.
 * @param line the line to add
 */
void
PolygonizeGraph::addEdge(const LineString *line) {
	if (line->isEmpty()) return;

	// TODO(vicky) fix CoordinateSequence::removeRepeatedPoints function
	auto linePts = CoordinateSequence::removeRepeatedPoints(line->getCoordinatesRO());

	/*
	 * This would catch invalid linestrings
	 * (containing duplicated points only)
	 */
	if ( linePts->getSize() < 2 ) {
		delete linePts;
		return;
	}

	const Coordinate& startPt = linePts->getAt(0);
	const Coordinate& endPt = linePts->getAt(linePts->getSize() - 1);
	auto nStart = getNode(startPt);
	auto nEnd = getNode(endPt);

	auto de0 = std::make_shared<PolygonizeDirectedEdge>(PolygonizeDirectedEdge(
			nStart, nEnd, linePts->getAt(1), true));

	auto de1 = std::make_shared<PolygonizeDirectedEdge>(PolygonizeDirectedEdge(
			nEnd, nStart, linePts->getAt(linePts->getSize() - 2), false));

	Edge *edge = new PolygonizeEdge(line);
	edge->setDirectedEdges(de0, de1);
	add(edge);

	m_newCoords.push_back(linePts);
}


void
PolygonizeGraph::computeNextCWEdges() {
	// set the next pointers for the edges around each node
	for(auto n : m_nodeMap) {
		computeNextCWEdges(n.second);
	}
}

/* private */
void
PolygonizeGraph::convertMaximalToMinimalEdgeRings(
		const DirectedEdges ringEdges) {
	for (const auto de : ringEdges) {
		auto label = de->getLabel();
		auto intNodes = findIntersectionNodes(de, label);
		for (auto &n : intNodes) {
			computeNextCCWEdges(n, label);
		}
	}
}

PolygonizeGraph::NodeVector
PolygonizeGraph::findIntersectionNodes(
		DirectedEdgePtr startDE,
		long label) const {
	NodeVector intNodes;
	auto de = startDE;

	do {
		auto node = de->getFromNode();
		if (node->getDegree<PolygonizeDirectedEdge>(label) > 1) {
			intNodes.push_back(node);
		}
		de = de->getNext();
		assert(de);  // found NULL DE in ring
		assert(de == startDE || !safe_cast<PolygonizeDirectedEdge *>(de)->isInRing());  // found DE already in ring
	} while (de != startDE);

	return intNodes;
}

/* public */
std::vector<PolygonizeGraph::EdgeRingPtr>
PolygonizeGraph::getEdgeRings() {
	std::vector<EdgeRingPtr> edgeRingList;
	// maybe could optimize this, since most of these pointers should
	// be set correctly already
	// by deleteCutEdges()
	// CVVC: the function is public,
	//  Q: is there a guaranty that deleteCutEdges has being called?
	computeNextCWEdges();

	// clear labels of all edges in graph
	label(m_dirEdges, -1);

	convertMaximalToMinimalEdgeRings(findLabeledEdgeRings(m_dirEdges));

	// find all edgerings
	for (auto de : m_dirEdges) {

		if (de->isMarked()) continue;
		if (safe_cast<PolygonizeDirectedEdge*>(de)->isInRing()) continue;

		edgeRingList.emplace_back(findEdgeRing(de));
	}
	return edgeRingList;
}

#if 0
/* public [[deprecated]] */
void
PolygonizeGraph::getEdgeRings(std::vector<EdgeRing*>& edgeRingList) {
	edgeRingList = getEdgeRings();
}
#endif

PolygonizeGraph::DirectedEdges
PolygonizeGraph::findLabeledEdgeRings(
		const DirectedEdges dirEdges) const
{
	DirectedEdges edgeRingStarts;

	// label the edge rings formed
	long currLabel(1);
	for (const auto de : dirEdges) {
		if (de->isMarked()) continue;
		if (de->getLabel() >= 0) continue;

		edgeRingStarts.push_back(de);

		auto edges = findDirEdgesInRing(de);

		label(edges, currLabel);

		++currLabel;
	}
	return edgeRingStarts;
}

/* public */
std::vector<const LineString*>
PolygonizeGraph::deleteCutEdges() {
	std::vector<const LineString*> cutLines;

	computeNextCWEdges();

	// label the current set of edgerings
	/* even that is a find it has side efects on the labels */
	findLabeledEdgeRings(m_dirEdges);  // ignoring the result

	/*
	 * Cut Edges are edges where both dirEdges have the same label.
	 * Delete them, and record them
	 */
	for (auto de : m_dirEdges) {
		if (de->isMarked()) continue;

		auto sym = de->getSym();

		if (de->getLabel() == sym->getLabel()) {
			de->setMarked(true);
			sym->setMarked(true);

			// save the line as a cut edge
			auto ce = dynamic_cast<PolygonizeEdge*>(de->parentEdge());

			cutLines.push_back(ce->getLine());
		}
	}
	return cutLines;
}

/* deprecated */
void
PolygonizeGraph::deleteCutEdges(std::vector<const LineString*> &cutLines) {
	cutLines = deleteCutEdges();
}

void
PolygonizeGraph::label(
		const DirectedEdges &dirEdges,
	 	long label) const {
	for(auto e : dirEdges) e->setLabel(label);
}

void
PolygonizeGraph::computeNextCWEdges(NodePtr node) {
	DirectedEdgePtr startDE = nullptr;
	DirectedEdgePtr prevDE = nullptr;

	// the edges are stored in CCW order around the star
	for(auto de : node->getSortedOutEdges()) {
		auto outDE = de;
		if (outDE->isMarked()) continue;
		if (!startDE) startDE = outDE;
		if (prevDE) {
      prevDE->getSym()->setNext(outDE);
		}
		prevDE = outDE;
	}
	if (prevDE) {
		prevDE->getSym()->setNext(startDE);
	}
}

/**
 * Computes the next edge pointers going CCW around the given node, for the
 * given edgering label.
 * This algorithm has the effect of converting maximal edgerings into
 * minimal edgerings
 */
void
PolygonizeGraph::computeNextCCWEdges(NodePtr node, long label) {

	DirectedEdgePtr firstOutDE = nullptr;
	DirectedEdgePtr prevInDE = nullptr;

	// the edges are stored in CCW order around the star
	auto edges = node->getSortedOutEdges();

	/*
	 * Cycling in reverse order.
	 */
  for(auto it = edges.rbegin(); it != edges.rend(); ++it)
  {
    auto de = *it;
		auto sym = de->getSym();

		auto outDE = (de->getLabel() == label)? de : nullptr;
		auto inDE = (sym->getLabel() == label)? sym : nullptr;

		if (!outDE && !inDE) continue;  // this edge is not in edgering

		if (inDE) prevInDE = inDE;

		if (outDE) {
			if (prevInDE) {
				prevInDE->setNext(outDE);
				prevInDE = nullptr;
			}
			if (!firstOutDE) firstOutDE = outDE;
		}
	}

	if (prevInDE) {
		assert(firstOutDE);
		prevInDE->setNext(firstOutDE);
	}
}

PolygonizeGraph::DirectedEdges
PolygonizeGraph::findDirEdgesInRing(DirectedEdgePtr startDE) const {
	auto de = startDE;
	DirectedEdges edges;

	do {
		edges.push_back(de);
		de = de->getNext();
		assert(de);  // found NULL DE in ring
		assert(de == startDE || !safe_cast<PolygonizeDirectedEdge *>(de)->isInRing());  // found DE already in ring
	} while (de != startDE);

	return edges;
}

PolygonizeGraph::EdgeRingPtr
PolygonizeGraph::findEdgeRing(DirectedEdgePtr startDE) const {
	auto de = startDE;

	EdgeRingPtr er( new EdgeRing(*m_factory) );

	do {
		er->add(de);

		auto e = safe_cast<PolygonizeDirectedEdge *>(de);
		e->setRing(er.get());

		de = de->getNext();
		assert(de);  // found NULL DE in ring
		assert(de == startDE || !safe_cast<PolygonizeDirectedEdge *>(de)->isInRing());  // found DE already in ring
	} while (de != startDE);

	return std::move(er);
}

/* public */
std::vector<const LineString*>
PolygonizeGraph::deleteDangles() {
	std::vector<const LineString*> dangleLines;

	auto nodeStack = findNodesOfDegree(1);

	std::set<const LineString*> uniqueDangles;

	while (!nodeStack.empty()) {
		auto node = nodeStack.back();
		nodeStack.pop_back();
		node->markAll(false);

		/* get sorted edges in ascending order by angle with the positive x-axis */
		for(auto de : node->getSortedOutEdges()) {
			// delete this edge and its sym
			de->setMarked(true);
			auto sym = de->getSym();
			if (sym) sym->setMarked(true);

			// save the line as a dangle
			auto e = dynamic_cast<PolygonizeEdge*>(de->parentEdge());
			auto ls(e->getLine());

			if (uniqueDangles.insert(ls).second) {
				dangleLines.push_back(ls);
			}

			// add the toNode to the list to be processed,
			// if it is now a dangle
			auto toNode = de->getToNode();
      if (toNode->getDegreeNonDeleted() == 1)
      {
        nodeStack.push_back(toNode);
      }
    }
  }
  return dangleLines;
}

/* deprecated */
void
PolygonizeGraph::deleteDangles(std::vector<const LineString*>& dangleLines) {
	dangleLines = deleteDangles();
}

}  // namespace polygonize
}  // namespace operation
}  // namespace geos

