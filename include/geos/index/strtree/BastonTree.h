//
// Created by dan on 11/26/20.
//

#pragma once

#include <geos/geom/Envelope.h>
#include <geos/geom/Geometry.h>
#include <geos/index/SpatialIndex.h> // for inheritance
#include <geos/index/chain/MonotoneChain.h>
#include <geos/index/ItemVisitor.h>
#include <geos/util.h>

namespace geos {
namespace index {
namespace strtree {

    template<typename ItemType>
    class BastonNode {
    private:
        union {
            const ItemType* item;
            const BastonNode<ItemType>* childrenEnd;
        } data;
        const BastonNode<ItemType>* children;

        geom::Envelope bounds;

        double sortVal = std::numeric_limits<double>::quiet_NaN();

    public:
        BastonNode() = delete;

        BastonNode(const ItemType* p_item, const geom::Envelope& env) {
            data.item = p_item;
            children = nullptr;
            bounds = env;
        }

        BastonNode(const BastonNode<ItemType>* begin, const BastonNode<ItemType>* end) {
            children = begin;
            data.childrenEnd = end;
            computeEnvelopeFromChildren();
        }

        void setSortVal(double d) {
            sortVal = d;
        }

        double getSortVal() const {
            return sortVal;
        }

        const BastonNode<ItemType>* beginChildren() const {
            return children;
        }

        const BastonNode<ItemType>* endChildren() const {
            return data.childrenEnd;
        }

        bool isLeaf() const {
            return children == nullptr;
        }

        bool envelopeIntersects(const geom::Envelope& queryEnv) const {
            return getEnvelope().intersects(queryEnv);
        }

        void computeEnvelopeFromChildren() {
            bounds.setToNull();
            for (auto* child = children; child < data.childrenEnd; ++child) {
                bounds.expandToInclude(child->getEnvelope());
            }
        }

        const geom::Envelope& getEnvelope() const {
            return bounds;
        }

        const ItemType* getItem() const {
            return data.item;
        }

        void removeItem() {
            data.item = nullptr;
        }

    };

    template<typename ItemType>
    class BastonTree : public SpatialIndex {
    public:
        explicit BastonTree(size_t p_nodeCapacity = 10) : root(nullptr), nodeCapacity(p_nodeCapacity) {}

        BastonTree(size_t p_nodeCapacity, size_t itemCapacity) : root(nullptr), nodeCapacity(p_nodeCapacity) {
            auto finalSize = treeSize(itemCapacity);
            nodes.reserve(finalSize);
        }

        bool built() const {
            return root != nullptr;
        }

        void insert(const geom::Envelope& itemEnv, const ItemType* item) {
            createLeafNode(item, itemEnv);
        }

        void insert(const geom::Envelope* itemEnv, void* item) override {
            insert(*itemEnv, static_cast<ItemType*>(item));
        }

        template<typename Visitor>
        void query(const geom::Envelope& queryEnv, Visitor&& visitor) {
            if (!built()) {
                build();
            }

            if (root->envelopeIntersects(queryEnv)) {
                if (root->isLeaf()) {
                    visitor(root->getItem());
                } else {
                    query(queryEnv, *root, visitor);
                }
            }
        }

        void query(const geom::Envelope* queryEnv, std::vector<void*> & results) override {
            query(*queryEnv, [&results](const ItemType* x) {
                results.push_back((void*) x);
            });
        }

        void query(const geom::Envelope* queryEnv, ItemVisitor& visitor) override {
            query(*queryEnv, [&visitor](const ItemType* x) {
                visitor.visitItem((void*) x);
            });
        }

        bool remove(const geom::Envelope* itemEnv, void* item) override {
            if (root == nullptr) {
                return false;
            }

            if (root->isLeaf()) {
                if (root->getItem() == item) {
                    root->removeItem();
                    return true;
                }
                return false;
            }

            return remove(*itemEnv, *root, static_cast<ItemType*>(item));
        }

    private:
        using Node = BastonNode<ItemType>;
        using NodeList = std::vector<Node>;
        using NodeListIterator = typename NodeList::iterator;

        NodeList nodes;
        BastonNode<ItemType>* root;
        size_t nodeCapacity;

        void createLeafNode(const ItemType* item, const geom::Envelope& env) {
            nodes.emplace_back(item, env);
        }

        void createBranchNode(const Node* begin, const Node* end) {
            assert(nodes.size() < nodes.capacity());
            nodes.emplace_back(begin, end);
        }

        void build() {
            auto finalSize = treeSize(nodes.size());
            nodes.reserve(finalSize);

            auto begin = nodes.begin();
            auto end = nodes.end();

            while (std::distance(begin, end) > 1) {
                createParentNodes(begin, end);
                begin = end;
                end = nodes.end();
            }

            assert(finalSize == nodes.size());

            root = &nodes.back();
        }

        size_t treeSize(size_t numLeafNodes) {
            size_t nodesInTree = numLeafNodes;

            size_t nodesWithoutParents = numLeafNodes;
            while (nodesWithoutParents > 1) {
                auto numSlices = sliceCount(nodesWithoutParents);
                auto nodesPerSlice = sliceCapacity(nodesWithoutParents, numSlices);

                size_t parentNodesAdded = 0;
                for (size_t j = 0; j < numSlices; j++) {
                    auto nodesInSlice = std::min(nodesWithoutParents, nodesPerSlice);
                    nodesWithoutParents -= nodesInSlice;

                    parentNodesAdded += static_cast<size_t>(std::ceil(static_cast<double>(nodesInSlice) / static_cast<double>(nodeCapacity)));
                }

                nodesInTree += parentNodesAdded;
                nodesWithoutParents = parentNodesAdded;
            }

            return nodesInTree;
        }

        void createParentNodes(const NodeListIterator& begin, const NodeListIterator& end) {
            auto numChildren = std::distance(begin, end);

            auto numSlices = sliceCount(numChildren);
            auto nodesPerSlice = sliceCapacity(numChildren, numSlices);

            setSortValuesX(begin, end);
            //sortNodesX(begin, end);

            auto startOfSlice = begin;
            for (size_t j = 0; j < numSlices; j++) {
                auto nodesRemaining = static_cast<size_t>(std::distance(startOfSlice, end));
                auto nodesInSlice = std::min(nodesRemaining, nodesPerSlice);
                auto endOfSlice = std::next(startOfSlice, nodesInSlice);

                partialSortNodes(startOfSlice, endOfSlice, end);

                addParentNodesFromVerticalSlice(startOfSlice, endOfSlice);

                startOfSlice = endOfSlice;
            }
        }

        void addParentNodesFromVerticalSlice(const NodeListIterator & begin, const NodeListIterator & end) {
            setSortValuesY(begin, end);
            //sortNodesY(begin, end);

            auto firstChild = begin;
            while (firstChild != end) {
                auto childrenRemaining = static_cast<size_t>(std::distance(firstChild, end));
                auto childrenForNode = std::min(nodeCapacity, childrenRemaining);
                auto lastChild = std::next(firstChild, childrenForNode);

                partialSortNodes(firstChild, lastChild, end);

                const Node* ptr_first = &*firstChild;
                const Node* ptr_end = ptr_first + childrenForNode;

                createBranchNode(ptr_first, ptr_end);
                firstChild = lastChild;
            }
        }

        void setSortValuesX(const NodeListIterator& begin, const NodeListIterator& end) {
            std::for_each(begin, end, [](Node& n) {
                const geom::Envelope& e = n.getEnvelope();
                n.setSortVal(e.getMinX() + e.getMaxX());
            });
        }

        void setSortValuesY(const NodeListIterator& begin, const NodeListIterator& end) {
            std::for_each(begin, end, [](Node& n) {
                const geom::Envelope& e = n.getEnvelope();
                n.setSortVal(e.getMinY() + e.getMaxY());
            });
        }

        void sortNodes(NodeListIterator& begin, NodeListIterator& end) {
            std::sort(begin, end, [](const Node &a, const Node &b) {
                return a.getSortVal() < b.getSortVal();
            });
        }

        // Partially sort nodes between `begin` and `end` such that all nodes less than `mid` are placed before `mid`.
        void partialSortNodes(const NodeListIterator& begin, const NodeListIterator& mid, const NodeListIterator& end) {
            std::nth_element(begin, mid, end, [](const Node& a, const Node& b) {
                return a.getSortVal() < b.getSortVal();
            });
        }

        template<typename Visitor>
        void query(const geom::Envelope& queryEnv,
                   const Node& node,
                   Visitor&& visitor) {

            assert(!node.isLeaf());

            for (auto* child = node.beginChildren(); child < node.endChildren(); ++child) {
                if(child->envelopeIntersects(queryEnv)) {
                    if (child->isLeaf() && child->getItem() != nullptr) {
                        visitor(child->getItem());
                    } else {
                        query(queryEnv, *child, visitor);
                    }
                }
            }
        }

        bool remove(const geom::Envelope& queryEnv,
                    const Node& node,
                    const ItemType* item) {

            assert(!node.isLeaf());

            for (auto* child = node.beginChildren(); child < node.endChildren(); ++child) {
                if(child->envelopeIntersects(queryEnv)) {
                    if (child->isLeaf() && child->getItem() == item) {
                        // const cast is ugly, but alternative seems to be to remove all
                        // const qualifiers in Node and open up mutability everywhere?
                        auto mutableChild = const_cast<Node*>(child);
                        mutableChild->removeItem();
                        return true;
                    } else {
                        bool removed = remove(queryEnv, *child, item);
                        if (removed) {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        size_t sliceCount(size_t numNodes) const {
            double minLeafCount = std::ceil(static_cast<double>(numNodes) / static_cast<double>(nodeCapacity));

            return static_cast<size_t>(std::ceil(std::sqrt(minLeafCount)));
        }

        static size_t sliceCapacity(size_t numNodes, size_t numSlices) {
            return static_cast<size_t>(std::ceil(static_cast<double>(numNodes) / static_cast<double>(numSlices)));
        }

    };
}
}
}
