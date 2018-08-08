/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
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
 * Last port: operation/linemerge/LineMerger.java r378 (JTS-1.12)
 *
 **********************************************************************/

#include <geos/operation/linemerge/LineMerger.h>
#include <geos/operation/linemerge/LineMergeDirectedEdge.h>
#include <geos/operation/linemerge/EdgeString.h>
#include <geos/planargraph/DirectedEdge.h>
#include <geos/planargraph/Edge.h>
#include <geos/planargraph/Node.h>
//#include <geos/planargraph/GraphComponent.h>
#include <geos/geom/GeometryComponentFilter.h>
#include <geos/geom/LineString.h>
#include <geos/planargraph/detail.hpp>

#include <cassert>
#include <functional>
#include <vector>

using namespace std;
using namespace geos::planargraph;
using namespace geos::geom;

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

namespace geos {
namespace operation { // geos.operation
namespace linemerge { // geos.operation.linemerge

void
LineMerger::add(vector<Geometry*> *geometries)
{
	for(auto geometry : *geometries) add(geometry);
}

LineMerger::LineMerger():
	mergedLineStrings(nullptr),
	factory(nullptr)
{
}

LineMerger::~LineMerger()
{
	for (auto e : edgeStrings) delete e;
}


struct LMGeometryComponentFilter: public GeometryComponentFilter {
	LineMerger *lm;

	LMGeometryComponentFilter(LineMerger *newLm): lm(newLm) {}

	void filter(const Geometry *geom) {
		const auto ls = dynamic_cast<const LineString *>(geom);
		if ( ls ) lm->add(ls);
	}
};


/**
 * Adds a Geometry to be processed. May be called multiple times.
 * Any dimension of Geometry may be added; the constituent linework will be
 * extracted.
 */
void
LineMerger::add(const Geometry *geometry)
{
	LMGeometryComponentFilter lmgcf(this);
	geometry->applyComponentFilter(lmgcf);
}

void
LineMerger::add(const LineString *lineString)
{
	if (factory==nullptr) factory=lineString->getFactory();
	graph.addEdge(lineString);
}

void
LineMerger::merge()
{
	if (mergedLineStrings!=nullptr) return;

	// reset marks (this allows incremental processing)
	::setMarkedMap(graph.getNodes(), false);
	setMarked(graph.getEdges(), false);

	for (auto e : edgeStrings) delete e;
	edgeStrings.clear();

	buildEdgeStringsForObviousStartNodes();
	buildEdgeStringsForIsolatedLoops();

	mergedLineStrings=new vector<LineString*>(edgeStrings.size());
#if 1
	for (size_t i=0; i<edgeStrings.size(); ++i)
	{
		EdgeString *edgeString=edgeStrings[i];
		(*mergedLineStrings)[i]=edgeString->toLineString();
	}
#else
	for (const auto e : edgeStrings)
	{
		mergedLineStrings->push_back(e->toLineString());
	}

#endif
}

void
LineMerger::buildEdgeStringsForObviousStartNodes()
{
	buildEdgeStringsForNonDegree2Nodes();
}

void
LineMerger::buildEdgeStringsForIsolatedLoops()
{
	buildEdgeStringsForUnprocessedNodes();
}

void
LineMerger::buildEdgeStringsForUnprocessedNodes()
{
#if GEOS_DEBUG
	cerr<<__FUNCTION__<<endl;
#endif

	for (auto &n : graph.getNodes()) {
		Node *node = n.second;
#if GEOS_DEBUG
		cerr<<"Node "<<i<<": "<<*node<<endl;
#endif
		if (!node->isMarked()) {
			assert(node->getDegree()==2);
			buildEdgeStringsStartingAt(node);
			node->setMarked(true);
#if GEOS_DEBUG
			cerr<<" setMarked(true) : "<<*node<<endl;
#endif
		}
	}
}

void
LineMerger::buildEdgeStringsForNonDegree2Nodes()
{
#if GEOS_DEBUG
	cerr<<__FUNCTION__<<endl;
#endif

	for (auto &n : graph.getNodes()) {
		Node* node = n.second;
#if GEOS_DEBUG
		cerr<<"Node "<<i<<": "<<*node<<endl;
#endif
		if (node->getDegree()!=2) {
			buildEdgeStringsStartingAt(node);
			node->setMarked(true);
#if GEOS_DEBUG
			cerr<<" setMarked(true) : "<<*node<<endl;
#endif
		}
	}
}

void
LineMerger::buildEdgeStringsStartingAt(Node *node)
{
	for (auto e : node->getOutEdges())
	{
		auto directedEdge = safe_cast<LineMergeDirectedEdge*> (e);
		if (directedEdge->parentEdge()->isMarked()) {
			continue;
		}
		edgeStrings.push_back(buildEdgeStringStartingWith(directedEdge));
	}
}

EdgeString*
LineMerger::buildEdgeStringStartingWith(LineMergeDirectedEdge *start)
{
	EdgeString *edgeString = new EdgeString(factory);
	LineMergeDirectedEdge *current=start;
	do {
		edgeString->add(current);
		current->parentEdge()->setMarked(true);
		current=current->getNext();
	} while (current!=nullptr && current!=start);
	return edgeString;
}

/**
 * Returns the LineStrings built by the merging process.
 */
vector<LineString*>*
LineMerger::getMergedLineStrings()
{
	merge();

	// Explicitly give ownership to the caller.
	vector<LineString*>* ret = mergedLineStrings;
	mergedLineStrings = nullptr;
	return ret;
}

} // namespace geos.operation.linemerge
} // namespace geos.operation
} // namespace geos
