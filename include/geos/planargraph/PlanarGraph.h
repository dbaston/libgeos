/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: planargraph/PlanarGraph.java rev. 107/138 (JTS-1.10)
 *
 **********************************************************************/

#ifndef GEOS_PLANARGRAPH_PLANARGRAPH_H
#define GEOS_PLANARGRAPH_PLANARGRAPH_H

#include <geos/export.h>
#include <geos/planargraph/NodeMap.h> // for composition

#include <vector> // for typedefs
#include <memory>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4251) // warning C4251: needs to have dll-interface to be used by clients of class
#endif

// Forward declarations
namespace geos {
	namespace geom {
		class Coordinate;
	}
	namespace planargraph {
		class DirectedEdge;
		class Edge;
		class Node;
	}
}

namespace geos {
namespace planargraph { // geos.planargraph

/**
 * \class PlanarGraph planargraph.h geos/planargraph.h
 * \brief
 * Represents a directed graph which is embeddable in a planar surface.
 *
 * This class and the other classes in this package serve as a framework for
 * building planar graphs for specific algorithms. This class must be
 * subclassed to expose appropriate methods to construct the graph. This allows
 * controlling the types of graph components (DirectedEdge, Edge and Node)
 * which can be added to the graph. An application which uses the graph
 * framework will almost always provide subclasses for one or more graph
 * components, which hold application-specific data and graph algorithms.
 */
class GEOS_DLL PlanarGraph {
public:
	typedef std::shared_ptr<DirectedEdge> DirectedEdgePtr;
	typedef std::vector<DirectedEdgePtr> DirectedEdges;

	typedef NodeMap::NodePtr NodePtr;
	typedef NodeMap::NodeVector NodeVector;


	typedef std::vector<Edge *> EdgeContainer;
	typedef EdgeContainer::iterator EdgeIterator;


	/**
	 * \brief
	 * Constructs a PlanarGraph without any Edges, DirectedEdges, or Nodes.
	 */
	PlanarGraph() {}

	virtual ~PlanarGraph() {}

	/**
	 * \brief
	 * Returns the Node at the given location,
	 * or null if no Node was there.
	 */
	NodePtr findNode(const geom::Coordinate& pt) {
		return m_nodeMap.find(pt);
	}

	const NodeMap getNodes() const {
		return m_nodeMap;
	}

#if 0
	NodeMap& getNodes() {
		return m_nodeMap;
	}
#endif

	/**
	 * \brief
	 * Returns the Edges that have been added to this PlanarGraph
	 * @see #add(Edge)
	 */
	std::vector<Edge*>& getEdges() {
		return m_edges;
	}

	DirectedEdges& getDirEdges() {
		return m_dirEdges;
	}

	/**
	 * \brief
	 * Removes an Edge and its associated DirectedEdges from their
	 * from-Nodes and from this PlanarGraph.
	 *
	 * Note: This method does not remove the Nodes associated
	 * with the Edge, even if the removal of the Edge reduces the
	 * degree of a Node to zero.
	 */
	void remove(Edge* edge);

	/**
	 * \brief
	 * Removes DirectedEdge from its from-Node and from this PlanarGraph.
	 *
	 * Note:
	 * This method does not remove the Nodes associated with the
	 * DirectedEdge, even if the removal of the DirectedEdge reduces
	 * the degree of a Node to zero.
	 */
	void remove(DirectedEdgePtr de);

	/**
	 * \brief
	 * Removes a node from the graph, along with any associated
	 * DirectedEdges and Edges.
	 */
	void remove(NodePtr node);

	/**
	 * \brief
	 * Returns all Nodes with the given number of Edges around it.
 	 * The return value is a newly allocated vector of existing nodes
	 */
	NodeVector findNodesOfDegree(std::size_t degree) const;

#ifdef GEOS_USEDEPRECATED
	/** @name deprecated */
	///@{
	/**
	 * \brief
	 * Returns an Iterator over the Nodes in this PlanarGraph.
	 */
	// [[deprecated]]
	NodeMap::container::iterator nodeIterator() {
		return m_nodeMap.begin();
	}

	// [[deprecated]]
	NodeMap::container::iterator nodeBegin() {
		return m_nodeMap.begin();
	}

	// [[deprecated]]
	NodeMap::container::const_iterator nodeBegin() const {
		return m_nodeMap.begin();
	}

	// [[deprecated]]
	NodeMap::container::iterator nodeEnd() {
		return m_nodeMap.end();
	}

	// [[deprecated]]
	NodeMap::container::const_iterator nodeEnd() const {
		return m_nodeMap.end();
	}

	// [[deprecated]]
	/**
	 * \brief
	 * Returns the Nodes in this PlanarGraph.
	 *
	 * @param nodes : the nodes are push_back'ed here
	 */
	void getNodes(std::vector<Node*>& nodes) const { nodes = m_nodeMap.getNodes(); }

	/**
	 * \brief
	 * Returns an Iterator over the DirectedEdges in this PlanarGraph,
	 * in the order in which they were added.
	 *
	 * @see add(Edge)
	 * @see add(DirectedEdge)
	 */
	// [[deprecated]]
	std::vector<DirectedEdge*>::iterator dirEdgeIterator() {
		return m_dirEdges.begin();
	}

	/// Alias for edgeBegin()
	std::vector<Edge*>::iterator edgeIterator() {
		return m_edges.begin();
	}

	/// Returns an iterator to first Edge in this graph.
	//
	/// Edges are stored in the order they were added.
	/// @see add(Edge)
	///
	// [[deprecated]]
	std::vector<Edge*>::iterator edgeBegin() {
		return m_edges.begin();
	}

	/// Returns an iterator to one-past last Edge in this graph.
	//
	/// Edges are stored in the order they were added.
	/// @see add(Edge)
	///
	// [[deprecated]]
	std::vector<Edge*>::iterator edgeEnd() {
		return m_edges.end();
	}

	/**
	 * \brief
	 * Get all Nodes with the given number of Edges around it.
	 *
 	 * Found nodes are pushed to the given vector
	 */
	// [[deprecated]]
	void findNodesOfDegree(std::size_t degree, std::vector<Node*>& to) const;
	///@}
#endif

protected:
	/**
	 * \brief
	 * Adds a node to the std::map, replacing any that is already at that
	 * location.
	 *
	 * Only subclasses can add Nodes, to ensure Nodes are
	 * of the right type.
	 * @return the added node
	 */
	void add(NodePtr node) {
		m_nodeMap.add(node);
	}

	/**
	 * \brief
	 * Adds the Edge and its DirectedEdges with this PlanarGraph.
	 *
	 * Assumes that the Edge has already been created with its associated
	 * DirectEdges.
	 * Only subclasses can add Edges, to ensure the edges added are of
	 * the right class.
	 */
	void add(Edge *edge);

	/**
	 * \brief
	 * Adds the Edge to this PlanarGraph.
	 *
	 * Only subclasses can add DirectedEdges,
	 * to ensure the edges added are of the right class.
	 */
	void add(DirectedEdgePtr dirEdge) {
		m_dirEdges.push_back(dirEdge);
	}

	std::vector<Edge*> m_edges;
	DirectedEdges m_dirEdges;
	NodeMap m_nodeMap;

};

} // namespace geos::planargraph
} // namespace geos

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // GEOS_PLANARGRAPH_PLANARGRAPH_H
