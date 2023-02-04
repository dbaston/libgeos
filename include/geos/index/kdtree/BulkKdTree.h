#pragma once

#include <geos/geom/Coordinate.h>
#include <geos/geom/Envelope.h>

#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace geos {
namespace index {
namespace kdtree {

class GEOS_DLL BulkKdTree {

public:
    static constexpr bool INITIAL_SPLIT_X = true;

    template<typename ItemType>
    class BulkNode {
    private:

        geom::CoordinateXY p; // 16
        ItemType data; // 8
        std::uint32_t left; // 8
        std::uint32_t right; // 8

    public:

        BulkNode(double p_x, double p_y, ItemType p_data) :
        p(p_x, p_y),
        data(p_data),
        left(0),
        right(0)
        {}

        BulkNode(const geom::CoordinateXY& p_p, ItemType p_data) :
        p(p_p),
        data(p_data),
        left(0),
        right(0)
        {}

        double getX() const { return p.x; }
        double getY() const { return p.y; }
        const geom::CoordinateXY& getCoordinate() const { return p; }
        const ItemType& getData() const { return data; }
        const BulkNode* getLeft(const BulkNode* root) const { return left ? root + static_cast<std::ptrdiff_t>(left) : nullptr; }
        const BulkNode* getRight(const BulkNode* root) const { return right ? root + static_cast<std::ptrdiff_t>(right) : nullptr; }
        BulkNode* getLeft() { return getLeft(); }
        BulkNode* getRight() { return getRight(); }
        void setLeft(const BulkNode* root, BulkNode* p_left) { left = static_cast<std::uint32_t>(p_left - root); }
        void setRight(const BulkNode* root, BulkNode* p_right) { right = static_cast<std::uint32_t>(p_right - root); }
    };

    using NodeType = BulkNode<void*>;


    void insert(const geom::CoordinateXY& coord, void* item) {
        if (built()) {
            throw std::runtime_error("Nope.");
        }
        m_nodes.emplace_back(coord, item);
    }

    void build() {
        // Remove duplicates
        std::sort(m_nodes.begin(), m_nodes.end(), [](const NodeType& a, const NodeType& b) {
            return a.getCoordinate().compareTo(b.getCoordinate()) < 0;
        });
        m_nodes.erase(std::unique(m_nodes.begin(), m_nodes.end(), [](const NodeType& a, const NodeType& b) {
            return a.getCoordinate().equals2D(b.getCoordinate());
        }), m_nodes.end());

        split(m_nodes.begin(), m_nodes.end(), INITIAL_SPLIT_X);
        m_built = true;
    }

    bool built() const {
        return m_built;
    }

    const NodeType*
    query(const geom::CoordinateXY& queryPt) const
    {
        if (!built()) {
            const_cast<BulkKdTree*>(this)->build();
        }

        return queryNodePoint(&m_nodes.front(), queryPt, INITIAL_SPLIT_X);
    }

#if 0
    void
    query(const geom::Envelope& queryEnv, KdNodeVisitor& visitor) const
    {
        query(queryEnv, [&visitor](const NodeType& node) {
            visitor.visit(&node);
        });
    }
#endif

    template<typename F>
    void
    query(const geom::Envelope& queryEnv, F&& visitor) const
    {
        if (!built()) {
            const_cast<BulkKdTree*>(this)->build();
        }

        // Non recursive formulation of in-order traversal from
        // http://web.cs.wpi.edu/~cs2005/common/iterative.inorder
        // Otherwise we may blow up the stack
        // See https://github.com/qgis/QGIS/issues/45226
        typedef std::pair<const NodeType*, bool> Pair;
        std::stack<Pair> activeNodes;
        const NodeType* root = &m_nodes.front();
        const NodeType* currentNode = root;
        bool splitX = INITIAL_SPLIT_X;

        while(true)
        {
            if( currentNode != nullptr )
            {
                double min;
                double discriminant;

                if (splitX) {
                    min = queryEnv.getMinX();
                    discriminant = currentNode->getX();
                } else {
                    min = queryEnv.getMinY();
                    discriminant = currentNode->getY();
                }
                bool searchLeft = min < discriminant;

                activeNodes.emplace(Pair(currentNode, splitX));

                // search is computed via in-order traversal
                const NodeType* leftNode = nullptr;
                if (searchLeft ) {
                    leftNode = currentNode->getLeft(root);
                }
                if( leftNode ) {
                    currentNode = leftNode;
                    splitX = !splitX;
                } else {
                    currentNode = nullptr;
                }
            }
            else if( !activeNodes.empty() )
            {
                currentNode = activeNodes.top().first;
                splitX = activeNodes.top().second;
                activeNodes.pop();

                if (queryEnv.contains(currentNode->getCoordinate())) {
                    visitor(*currentNode);
                }

                double max;
                double discriminant;

                if (splitX) {
                    max = queryEnv.getMaxX();
                    discriminant = currentNode->getX();
                } else {
                    max = queryEnv.getMaxY();
                    discriminant = currentNode->getY();
                }
                bool searchRight = discriminant <= max;

                if (searchRight) {
                    currentNode = currentNode->getRight(root);
                    if( currentNode )
                        splitX = !splitX;
                } else {
                    currentNode = nullptr;
                }
            }
            else
            {
                break;
            }
        }
    }


    void validateConstruction() {
        validateConstruction(m_nodes.front(), INITIAL_SPLIT_X);

        std::size_t hits = 0;
        visitSubtree(m_nodes.front(), [&hits](const NodeType& n) {
            (void) n;
            hits++;
        });

        if (hits != m_nodes.size()) {
            throw std::runtime_error("bad");
        }
    }

    void validateConstruction(const NodeType& node, bool splitX) {
        const NodeType* root = &m_nodes.front();

        if (splitX) {
            double split = node.getX();

            visitSubtree(*node.getLeft(root), [split](const NodeType& n) {
                if (n.getX() > split) {
                    throw std::runtime_error("x > split");
                }
            });
            visitSubtree(*node.getRight(root), [split](const NodeType& n) {
                if (n.getX() <= split) {
                    throw std::runtime_error("x <= split");
                }
            });
        } else {
            double split = node.getX();

            visitSubtree(*node.getLeft(root), [split](const NodeType& n) {
                if (n.getY() > split) {
                    throw std::runtime_error("y > split");
                }
            });
            visitSubtree(*node.getRight(root), [split](const NodeType& n) {
                if (n.getY() <= split) {
                    throw std::runtime_error("y <= split");
                }
            });
        }
    }

    template<typename F>
    void visitSubtree(const NodeType& n, F&& f) {
        const NodeType* root = &m_nodes.front();
        f(n);
        if (n.getLeft(root)) {
            visitSubtree(*n.getLeft(root), f);
        }
        if (n.getRight(root)) {
            visitSubtree(*n.getRight(root), f);
        }
    }

private:
    std::vector<NodeType> m_nodes;
    bool m_built = false;

    using NodeIterator = decltype(m_nodes.begin());

    void split(const NodeIterator& from, const NodeIterator& to, bool splitX) {
        auto distance = std::distance(from, to);

        if (distance <= 1) {
            return;
        }

        // Partially sort the range so that each item beyond midpoint has x/y >= midpoint
        auto midpoint = std::next(from, distance / 2);

#ifdef GEOS_DEBUG
        std::cerr << std::endl;
        std::cerr << "Range is " << std::distance(m_nodes.begin(), from) << " - " << std::distance(m_nodes.begin(), to) - 1 << std::endl;
        std::cerr << "Distance is " << distance << std::endl;
        std::cerr << "Midpoint at " << std::distance(m_nodes.begin(), midpoint) << std::endl;
#endif

        if (splitX) {
            std::nth_element(from, midpoint, to, [](const NodeType& a, const NodeType& b) {
                return a.getX() < b.getX();
            });
        } else {
            std::nth_element(from, midpoint, to, [](const NodeType& a, const NodeType& b) {
                return a.getY() < b.getY();
            });
        }

        // Move the midpoint to the beginning of the range and use it as the root element
        std::swap(*from, *midpoint);
        NodeType& subtree = *from;


        // Set the left and right nodes to the beginning of each partition
        const NodeType* root = &m_nodes.front();
        NodeType* left = &*std::next(from);
        NodeType* right = &*std::next(midpoint);

#ifdef GEOS_DEBUG
        std::cerr << "Split " << (splitX ? "X" : "Y") << " at " << (splitX ? root.getX() : root.getY()) << std::endl;
        std::cerr << "Left from " << std::distance(m_nodes.begin(), std::next(from)) << " - " << std::distance(m_nodes.begin(), std::next(midpoint)) - 1 << std::endl;
#endif

        subtree.setLeft(root, left);
        if (std::distance(std::next(midpoint), to) >= 1) {
#ifdef GEOS_DEBUG
            std::cerr << "Right from " << std::distance(m_nodes.begin(), std::next(midpoint)) << " - " << std::distance(m_nodes.begin(), to) - 1 << std::endl;
#endif
            subtree.setRight(root, right);
        }

        split(std::next(from), std::next(midpoint), !splitX);
        split(std::next(midpoint), to, !splitX);
    }

    const NodeType*
    queryNodePoint(const NodeType* currentNode, const geom::CoordinateXY& queryPt, bool odd) const
    {
        const NodeType* root = &m_nodes.front();

        while (currentNode != nullptr)
        {
            if (currentNode->getCoordinate().equals2D(queryPt))
                return currentNode;

            double ord;
            double discriminant;
            if (odd) {
                ord = queryPt.x;
                discriminant = currentNode->getX();
            }
            else {
                ord = queryPt.y;
                discriminant = currentNode->getY();
            }

            bool searchLeft = (ord < discriminant);
            odd = !odd;
            if (searchLeft) {
                currentNode = currentNode->getLeft(root);
            }
            else {
                currentNode = currentNode->getRight(root);
            }
        }
        return nullptr;
    }
};

}
}
}

