/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: planargraph/PlanarGraph.java rev. 107/138 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/planargraph/PlanarGraph.h>
#include <geos/planargraph/Edge.h>
#include <geos/planargraph/NodeMap.h>
#include <geos/planargraph/Node.h>
#include <geos/planargraph/DirectedEdge.h>

#include <algorithm>
#include <vector>
#include <map>

using namespace std;

namespace geos {
namespace planargraph {

template <typename T>
void
find_and_erase(T what, vector<T> &where) {
	if (!what) return;
	where.erase(std::remove(where.begin(), where.end(), what), where.end());
}


/*
 * Adds the Edge and its DirectedEdges with this PlanarGraph.
 * Assumes that the Edge has already been created with its associated
 * DirectEdges.
 * Only subclasses can add Edges, to ensure the edges added are of
 * the right class.
 */
void
PlanarGraph::add(Edge *edge)
{
	m_edges.push_back(edge);
	add(edge->getDirEdge(0));
	add(edge->getDirEdge(1));
}


/*
 * Removes an Edge and its associated DirectedEdges from their from-Nodes and
 * from this PlanarGraph. Note: This method does not remove the Nodes associated
 * with the Edge, even if the removal of the Edge reduces the degree of a
 * Node to zero.
 */
void
PlanarGraph::remove(Edge *edge)
{
	remove(edge->getDirEdge(0));
	remove(edge->getDirEdge(1));
	find_and_erase(edge, m_edges);
}

/*
 * Removes DirectedEdge from its from-Node and from this PlanarGraph. Note:
 * This method does not remove the Nodes associated with the DirectedEdge,
 * even if the removal of the DirectedEdge reduces the degree of a Node to
 * zero.
 */
void
PlanarGraph::remove(DirectedEdge *de)
{
	DirectedEdge *sym = de->getSym();
	if (sym!=nullptr) sym->setSym(nullptr);
	de->getFromNode()->getOutEdges().remove(de);
	find_and_erase(de, m_dirEdges);
}

/*
 * Removes a node from the graph, along with any associated
 * DirectedEdges and Edges.
 */
void
PlanarGraph::remove(Node *node)
{
	// unhook all directed edges
	auto & outEdges=node->getOutEdges().getEdges();
	for(auto de : outEdges) {
		DirectedEdge *sym = de->getSym();
		// remove the diredge that points to this node
		if (sym) remove(sym);
		find_and_erase(de, m_dirEdges);
		find_and_erase(de->getEdge(), m_edges);
	}
	// remove the node from the graph
	m_nodeMap.remove(node->getCoordinate());
	//nodes.remove(node);
}

/*public*/
void
PlanarGraph::findNodesOfDegree(size_t degree, vector<Node*>& nodesFound) const
{
	nodesFound = findNodesOfDegree(degree);
}

/*public*/
vector<Node*>
PlanarGraph::findNodesOfDegree(size_t degree) const
{
	vector<Node*> nodesFound;
	for (const auto &n : m_nodeMap)
	{
		auto node = n.second;
		if (node->getDegree() == degree) nodesFound.push_back(node);
	}
	return nodesFound;
}

} // namespace planargraph
} // namespace geos

