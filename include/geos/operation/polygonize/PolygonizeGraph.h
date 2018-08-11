/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2018 ~    Vicky Vergara
 * Copyright (C) 2006 Refractions Research Inc.
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

#ifndef GEOS_OP_POLYGONIZE_POLYGONIZEGRAPH_H
#define GEOS_OP_POLYGONIZE_POLYGONIZEGRAPH_H

#include <geos/export.h>

#include <geos/planargraph/PlanarGraph.h>  // for inheritance
#include <geos/geom/GeometryFactory.h>

#include <vector>
#include <memory>

#ifdef _MSC_VER
#pragma warning(push)
// warning C4251: needs to have dll-interface to be used by clients of class
#pragma warning(disable: 4251)
#endif

// Forward declarations

namespace geos {
namespace geom {
	class LineString;
	class GeometryFactory;
	class CoordinateSequence;
}  // namespace geom
}  // namespace geos

using geos::geom::GeometryFactory;
using geos::geom::CoordinateSequence;
using geos::geom::LineString;

namespace geos {
namespace operation {
namespace polygonize {


// Forward declarations
class EdgeRing;

/** \brief
 * Represents a planar graph of edges that can be used to compute a
 * polygonization, and implements the algorithms to compute the
 * EdgeRings formed by the graph.
 *
 * The marked flag on DirectedEdge is used to indicate that a directed edge
 * has be logically deleted from the graph.
 *
 */
class GEOS_DLL PolygonizeGraph: public planargraph::PlanarGraph {
public:
	bool empty() const {return m_nodeMap.empty();}
	/**
	 * \brief
	 * Deletes all edges at a node
	 */
	void deleteAllEdges(NodePtr node);

	/**
	 * \brief
	 * Create a new polygonization graph.
	 */
	PolygonizeGraph() {m_factory = new geom::GeometryFactory();};
	explicit PolygonizeGraph(const geom::GeometryFactory *newFactory);

	/**
	 * \brief
	 * Destroy a polygonization graph.
	 */
	~PolygonizeGraph() override;

	/**
	 * \brief
	 * Add a LineString forming an edge of the polygon graph.
	 * @param line the line to add
	 */
	void addEdge(const LineString *line);

	/**
	 * \brief
	 * Computes the EdgeRings formed by the edges in this graph.
	 *
	 * Any old values on edgeRingList will be deleted (not destroyed)
	 *
	 * @note Marked as deprecated, because its public
	 *
	 * @parami[in/out] edgeRingList : the EdgeRings found by the
	 * 	polygonization process will be pushed here.
	 *
	 */
	// [[deprecated]]
	void getEdgeRings(std::vector<EdgeRing*>& edgeRingList);

	/**
	 * \brief
	 * Computes the EdgeRings formed by the edges in this graph.
	 *
	 * @return edgeRingList: the EdgeRings found by the polygonization process.
	 */
	std::vector<EdgeRing*> getEdgeRings();

	/**
	 * \brief
	 * Finds and removes all cut edges from the graph.
	 *
	 * @return  vector with the list of the LineString forming the removed
	 *                   cut edges will be pushed here.
	 */
	std::vector<const LineString*> deleteCutEdges();
	// [[deprecated]]
	void deleteCutEdges(std::vector<const LineString*> &cutLines);

	/** \brief
	 * Marks all edges from the graph which are "dangles".
	 *
	 * Dangles are which are incident on a node with degree 1.
	 * This process is recursive, since removing a dangling edge
	 * may result in another edge becoming a dangle.
	 * In order to handle large recursion depths efficiently,
	 * an explicit recursion stack is used
	 *
	 * @param dangleLines : the LineStrings that formed dangles will
	 *                      be push_back'ed here
	 */
	std::vector<const LineString*> deleteDangles();
	// [[deprecated]]
	void deleteDangles(std::vector<const LineString*> &dangleLines);

 private:
	int getDegreeNonDeleted(NodePtr node) const;

	int getDegree(NodePtr node, long label) const;


	/**
	 * \brief
	 * Convert the maximal edge rings found by the initial graph traversal
	 * into the minimal edge rings required by JTS polygon topology rules.
	 *
	 * @param ringEdges
	 * 	the list of start edges for the edgeRings to convert.
	 *
	 */
	void convertMaximalToMinimalEdgeRings(
			const DirectedEdges ringEdges);

	/**
	 * \brief
	 * Finds all nodes in a maximal edgering
	 * which are self-intersection nodes
	 *
	 * @param startDE
	 * @param label
	 * @param intNodes : intersection nodes found will be pushed here
	 *                   the vector won't be cleared before pushing.
	 */
	NodeVector
	findIntersectionNodes(
		 	DirectedEdgePtr startDE,
			long label) const;

	/**
	 * Finds and labels all edgerings in the graph.
	 *
	 * The edge rings are labelling with unique integers.
	 * The labelling allows detecting cut edges.
	 *
	 * @param dirEdgesIn  a list of the DirectedEdges in the graph
	 * @result vector that contains each ring found
	 */
	DirectedEdges
	findLabeledEdgeRings(
			const DirectedEdges dirEdgesIn) const;

	void label(const DirectedEdges &dirEdges, long label) const;


	/**
	 * \brief
	 * Computes the next edge pointers going CCW around the given node,
	 * for the given edgering label.
	 * This algorithm has the effect of converting maximal edgerings
	 * into minimal edgerings
	 */
	void computeNextCWEdges();
	void computeNextCWEdges(NodePtr node);
	void computeNextCCWEdges(NodePtr node, long label);

	/**
	 * \brief
	 * Traverse a ring of DirectedEdges, accumulating them into a list.
	 * This assumes that all dangling directed edges have been removed
	 * from the graph, so that there is always a next dirEdge.
	 *
	 * @param startDE the DirectedEdge to start traversing at
	 * @result a vector of the DirectedEdge that form a ring
	 */
	DirectedEdges
	findDirEdgesInRing(DirectedEdgePtr startDE) const;

	/* has side effect of saving the Edge Ring found */
	EdgeRing* findEdgeRing(DirectedEdgePtr startDE) const;

	/*
	 * Data members
	 */
	const geom::GeometryFactory *m_factory;


	/*
	 *  These are for memory management
	 */
	/* created as PolygonizeDirectedEdge but saved as DirectedEdge */
	DirectedEdges m_newDirEdges;
	mutable std::vector<EdgeRing *> m_newEdgeRings;
	std::vector<geom::CoordinateSequence *> m_newCoords;
};

}  // namespace polygonize
}  // namespace operation
}  // namespace geos

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif  // GEOS_OP_POLYGONIZE_POLYGONIZEGRAPH_H
