/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 * Copyright (C) 2005-2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_PLANARGRAPH_NODEMAP_H
#define GEOS_PLANARGRAPH_NODEMAP_H

#include <geos/export.h>
#include <geos/geom/Coordinate.h> // for use in container

#include <map>
#include <vector>
#include <memory>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4251) // warning C4251: needs to have dll-interface to be used by clients of class
#endif

namespace geos {
namespace planargraph {

// Forward declarations
class Node;

/**
 * \brief
 * A map of Node, indexed by the coordinate of the node.
 *
 * Will not:
 *  - create a new node
 *  - destroy the nodes
 *
 * Classes using this NodeMap are responsible to do
 * - pointer cleanup when needed
 * - pass the responability to another class
 */
class GEOS_DLL NodeMap {
public:
#if 0
	typedef Node* NodePtr;
#else
	typedef std::shared_ptr<Node> NodePtr;
#endif
	typedef std::map<geom::Coordinate, NodePtr, geom::CoordinateLessThen> NodeContainer;
	typedef std::vector<NodePtr> NodeVector;

private:
	NodeContainer nodeMap;

public:
	/**
	 * \brief Constructs a NodeMap without any Nodes.
	 */
	NodeMap() = default;
	~NodeMap() = default;

	NodeContainer& getNodeMap();


	/**
	 * \brief
	 * Adds a node to the std::map, replacing any that is already
	 * at that location.
	 * @return the added node
	 */
	NodePtr add(NodePtr n);

	/**
	 * \brief
	 * Removes the Node at the given location, and returns it
	 * (or null if no Node was there).
	 */
	NodePtr remove(geom::Coordinate& pt);

	/**
	 * \brief
	 * Returns the Node at the given location,
	 * or null if no Node was there.
	 */
	NodePtr find(const geom::Coordinate& coord) const;

	/**
	 * \brief
	 * Returns an Iterator over the Nodes in this NodeMap,
	 * sorted in ascending order
	 * by angle with the positive x-axis.
	 */
#ifdef GEOS_USEDEPRECATED
	container::iterator iterator() {
		return nodeMap.begin();
	}
#endif

	NodeContainer::iterator begin() {
		return nodeMap.begin();
	}
	NodeContainer::const_iterator begin() const {
		return nodeMap.begin();
	}

	NodeContainer::iterator end() {
		return nodeMap.end();
	}
	NodeContainer::const_iterator end() const {
		return nodeMap.end();
	}

	bool empty() const {return nodeMap.empty();}

	/**
	 * \brief
	 * Returns the Nodes in this NodeMap, sorted in ascending order
	 * by angle with the positive x-axis.
	 *
	 * @return a vector of Node pointers
	 */
	NodeVector getNodes() const;

	/**
	 * \brief
	 * Returns the Nodes in this NodeMap, sorted in ascending order
	 * by angle with the positive x-axis.
	 *
	 * @param nodes : the nodes are push_back'ed here
	 */
	// [[deprecated]]
	void getNodes(NodeVector& nodes) const;
};


}  // namespace planargraph
}  // namespace geos

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // GEOS_PLANARGRAPH_NODEMAP_H
