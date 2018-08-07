/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/platform.h>
#include <geos/planargraph/algorithm/ConnectedSubgraphFinder.h>
#include <geos/planargraph/Subgraph.h>
#include <geos/planargraph/Edge.h>
#include <geos/planargraph/Node.h>
#include <geos/planargraph/DirectedEdge.h>
#include <geos/planargraph/DirectedEdgeStar.h>
#include <geos/planargraph/detail.hpp>

#include <cmath>
#include <vector>
#include <stack>

using namespace std;

namespace geos {
namespace planargraph {
namespace algorithm {

void
ConnectedSubgraphFinder::getConnectedSubgraphs(vector<Subgraph *>& subgraphs)
{
	::setVisitedMap(graph.getNodes(), false);
	for (auto &e : graph.getEdges())
	{
		auto node = e->getDirEdge(0)->getFromNode();
		if (!node->isVisited())
	 	{
			subgraphs.push_back(findSubgraph(node));
		}
	}
}

/*private*/
Subgraph*
ConnectedSubgraphFinder::findSubgraph(Node* node)
{
	Subgraph* subgraph = new Subgraph(graph);
	addReachable(node, subgraph);
	return subgraph;
}

/*private*/
void
ConnectedSubgraphFinder::addReachable(Node* startNode,
		Subgraph* subgraph)
{
	stack<Node *> nodeStack;
	nodeStack.push(startNode);
	while ( !nodeStack.empty() )
	{
		Node* node = nodeStack.top();
		nodeStack.pop();
		addEdges(node, nodeStack, subgraph);
	}
}

/*private*/
void
ConnectedSubgraphFinder::addEdges(Node* node,
		stack<Node *>& nodeStack, Subgraph* subgraph)
{
	node->setVisited(true);
	auto des = node->getOutEdges();
	for (auto de : des)
	{
		subgraph->add(de->parentEdge());
		Node *toNode = de->getToNode();
		if ( ! toNode->isVisited() ) nodeStack.push(toNode);
	}
}


} // namespace geos.planargraph.algorithm
} // namespace geos.planargraph
} // namespace geos
