/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2019 Daniel Baston
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 *
 **********************************************************************/

#ifndef GEOS_INDEX_KDTREE_KDNODE_H
#define GEOS_INDEX_KDTREE_KDNODE_H

#include <geos/geom/Coordinate.h>

#include <memory>

namespace geos {
namespace index {
namespace kdtree {

template<typename T=void*>
class KdNode {
public:
    /**
     * Create a new kdNode
     *
     * @param p location of node
     * @param _data object to associate with this node
     */
    KdNode(const geom::Coordinate& _p, T& _data) :
            p(_p),
            count(1),
            data(_data)
    {}

    /**
     * Return the X coordinate of the node
     * @return X coordinate of the node
     */
    double getX() const {
        return p.x;
    }

    /**
     * Return the Y coordinate of the node
     * @return Y coordinate of the node
     */
    double getY() const {
        return p.y;
    }

    /**
     * Return the location of the node
     * @return coordinate of the node
     */
    const geom::Coordinate& getCoordinate() const {
        return p;
    }

    /**
     * Return the data associated with the node
     * @return data associated with the node
     */
    const T& getData() const {
        return data;
    }

    KdNode<T>* getLeft() {
        return left.get();
    }

    KdNode<T>* getRight() {
        return right.get();
    }

    void increment() {
        count++;
    }

    /**
     * Get the number of inserted points that are coincident at this location
     *
     * @return number of inserted points that this node represents
     */
    size_t getCount() const {
        return count;
    }

    /**
     * Tests whether more than one point at this location (within the tolerance) has been inserted
     *
     * @return true of more than one point has been inserted at this location
     */
    bool isRepeated() const {
        return count > 1;
    }

    void setLeft(std::unique_ptr<KdNode> && _left) {
        left = std::move(_left);
    }

    void setRight(std::unique_ptr<KdNode> && _right) {
        right = std::move(_right);
    }


private:
    geom::Coordinate p;
    std::unique_ptr<KdNode<T>> left;
    std::unique_ptr<KdNode<T>> right;
    size_t count;
    T data;
};



}
}
}

#endif
