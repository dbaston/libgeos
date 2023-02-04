/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2020 Paul Ramsey <pramsey@cleverelephant.ca>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: index/kdtree/Node.java rev 1.8 (JTS-1.10)
 *
 **********************************************************************/

#pragma once

#include <geos/geom/Coordinate.h> // for composition

namespace geos {
namespace index { // geos::index
namespace kdtree { // geos::index::kdtree

/**
 * A node of a {@link KdTree}, which represents one or more points in the same location.
 */
class GEOS_DLL KdNode {

private:

    geom::CoordinateXY p; // 16
    void* data; // 8
    KdNode* left; // 8
    KdNode* right; // 8
    std::size_t count; // 8

public:

    KdNode(double p_x, double p_y, void* p_data) :
        p(p_x, p_y), data(p_data),
        left(0),
        right(0),
        count(1)
    {}

    KdNode(const geom::CoordinateXY& p_p, void* p_data) :
        p(p_p),
        data(p_data),
        left(0),
        right(0),
        count(1)
    {}

    double getX() const { return p.x; }
    double getY() const { return p.y; }
    const geom::CoordinateXY& getCoordinate() const { return p; }
    void* getData() const { return data; }
    const KdNode* getLeft() const { return left; }
    const KdNode* getRight() const { return right; }
    KdNode* getLeft() { return left; }
    KdNode* getRight() { return right; }
    void increment() { count++; }
    std::size_t getCount() const { return count; }
    bool isRepeated() const { return count > 1; }
    void setLeft(KdNode* p_left) { left = p_left; }
    void setRight(KdNode* p_right) { right = p_right; }

};

} // namespace geos::index::kdtree
} // namespace geos::index
} // namespace geos



