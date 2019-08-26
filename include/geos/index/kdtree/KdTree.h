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

#ifndef GEOS_INDEX_KDTREE_KDTREE_H
#define GEOS_INDEX_KDTREE_KDTREE_H

#include <geos/geom/Envelope.h>
#include <geos/index/kdtree/KdNode.h>
#include <geos/util.h>

namespace geos {
namespace index {
namespace kdtree {

template<typename T=void*>
class KdTree {
public:
    explicit KdTree(double _tolerance = 0.0) : tolerance(_tolerance) {}

    bool isEmpty() const {
        return root == nullptr;
    }

    KdNode<T>* insert(const geom::Coordinate & p) {
        if (std::is_pointer<T>::value) {
            return insert(p, nullptr);
        } else{
            T x;
            return insert(p, x);
        }
    }

    KdNode<T>* insert(const geom::Coordinate & p, T data) {
        if (isEmpty()) {
            root = new KdNode<T>(p, data);
            return root;
        }

        // Check if the point is already in the tree, within the tolerance.
        // If tolerance is zero, this phase of the insertion can be skipped.
        if (tolerance > 0) {
            KdNode<T>* matchNode = findBestMatchNode(p);
            if (matchNode != nullptr) {
                // point already in index - increment counter.
                matchNode->increment();
                return matchNode;
            }
        }

        return insertExact(p, data);
    }

    void query(const geom::Envelope& queryEnv, std::function<void(KdNode<T>*)> visitor) {
        queryNode(root, queryEnv, true, visitor);
    }

    void query(const geom::Envelope & queryEnv, std::vector<const KdNode<T>*> & result) {
        queryNode(root, queryEnv, true, [&result](KdNode<T>* node){ result.push_back(node); });
    }

    std::vector<const KdNode<T>*> query(const geom::Envelope & queryEnv) {
        std::vector<const KdNode<T>*> nodes;
        query(queryEnv, nodes);
        return nodes;
    }

private:
    KdNode<T>* insertExact(const geom::Coordinate & p, T & data) {
        KdNode<T>* currentNode = root;
        KdNode<T>* leafNode = root;
        bool isOddLevel = true;
        bool isLessThan = true;

        // Traverse the tree, first cutting the plane left-right (by X ordinate)
        // then top-bottom (by Y ordinate)
        while (currentNode != nullptr) {
            // check if point is already in tree and if so
            // simply return existing node
            if (p == currentNode->getCoordinate()) {
                currentNode->increment();
                return currentNode;
            }

            if (isOddLevel) {
                isLessThan = p.x < currentNode->getX();
            } else {
                isLessThan = p.y < currentNode->getY();
            }
            leafNode = currentNode;
            if (isLessThan) {
                currentNode = currentNode->getLeft();
            } else {
                currentNode = currentNode->getRight();
            }

            isOddLevel = !isOddLevel;
        }


        // no node found, add new leaf node to tree
        numberOfNodes++;
        auto node = detail::make_unique<KdNode<T>>(p, data);

        if (isLessThan) {
            leafNode->setLeft(std::move(node));
        } else {
            leafNode->setRight(std::move(node));
        }

        return node.get();
    }

    void queryNode(KdNode<T>* currentNode,
            const geom::Envelope& queryEnv,
            bool odd,
            std::function<void(KdNode<T>*)> visitor) {

        if (currentNode == nullptr) {
            return;
        }

        double min, max, discriminant;
        if (odd) {
            min = queryEnv.getMinX();
            max = queryEnv.getMaxX();
            discriminant = currentNode->getX();
        } else {
            min = queryEnv.getMinY();
            max = queryEnv.getMaxY();
            discriminant = currentNode->getY();
        }

        bool searchLeft = min < discriminant;
        bool searchRight = discriminant <= max;

        // search is computed via in-order traversal
        if (searchLeft) {
            queryNode(currentNode->getLeft(), queryEnv, !odd, visitor);
        }
        if (queryEnv.contains(currentNode->getCoordinate())) {
            visitor(currentNode);
        }
        if (searchRight) {
            queryNode(currentNode->getRight(), queryEnv, !odd, visitor);
        }
    }

    KdNode<T>* findBestMatchNode(const geom::Coordinate & p ) {
        BestMatchVisitor visitor(p, tolerance);
        query(visitor.queryEnvelope(), [&visitor](KdNode<T>* node) {
            visitor.visit(node);
        });

        return visitor.getNode();
    }

    class BestMatchVisitor {
    public:
        BestMatchVisitor(const geom::Coordinate & _p, double _tolerance) :
            p(_p),
            tolerance(_tolerance),
            matchNode(nullptr),
            matchDist(std::numeric_limits<double>::infinity()) {}

        KdNode<T>* getNode() {
            return matchNode;
        }

        geom::Envelope queryEnvelope() {
            geom::Envelope queryEnv(p);
            queryEnv.expandBy(tolerance);
            return queryEnv;
        }

        void visit(KdNode<T>* node) {
            double dist = p.distance(node->getCoordinate());

            if (dist > tolerance) {
                return;
            }

            // find closest match to p; break distance ties with Coordinate comparitor
            if (matchNode == nullptr ||
                dist < matchDist ||
                (dist == matchDist && node->getCoordinate().compareTo(matchNode->getCoordinate() < 1))) {
                matchNode = node;
                matchDist = dist;
            }
        }

    private:
        geom::Coordinate p;
        double tolerance;
        KdNode<T>* matchNode;
        double matchDist;
    };


    KdNode<T>* root;
    double tolerance;
    size_t numberOfNodes;
};

}
}
}


#endif
