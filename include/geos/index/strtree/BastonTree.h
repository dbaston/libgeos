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

    static const geom::Envelope& getEnvelope(const geom::Geometry* g) {
        return *(g->getEnvelopeInternal());
    }

    static const geom::Envelope& getEnvelope(const index::chain::MonotoneChain* mc) {
        // TODO make MonotoneChain getEnvelope const (using mutable envelope member)
        return const_cast<index::chain::MonotoneChain*>(mc)->getEnvelope();
    }

    template<typename ItemType>
    class BastonNode {
    private:
        const ItemType* item;
        std::vector<const BastonNode*> childNodes; // TODO replace with C array
        mutable std::unique_ptr<geom::Envelope> childBounds;

    public:
        BastonNode() : item(nullptr), childBounds(nullptr) {}

        explicit BastonNode(const ItemType* p_item) : item(p_item) {}

        decltype(childNodes.cbegin()) beginChildren() const {
            return childNodes.cbegin();
        }

        decltype(childNodes.cend()) endChildren() const {
            return childNodes.cend();
        }

        bool isLeaf() const {
            return item != nullptr;
        }

        bool envelopeIntersects(const geom::Envelope& queryEnv) const {
            return getEnvelope().intersects(queryEnv);
        }

        const geom::Envelope& getEnvelope() const {
            if (isLeaf()) {
                return geos::index::strtree::getEnvelope(item);
            } else {
                if (childBounds == nullptr) {
                    childBounds = detail::make_unique<geom::Envelope>();
                    for (const auto& child : childNodes) {
                        childBounds->expandToInclude(child->getEnvelope());
                    }
                }

                return *childBounds;
            }
        }


        const ItemType* getItem() const {
            return item;
        }

        void addChildNode(const BastonNode<ItemType>* child) {
            childNodes.push_back(child);
        }

        size_t numChildren() const {
            return childNodes.size();
        }


    };

    template<typename ItemType>
    class BastonTree : public SpatialIndex {
    public:
        BastonTree() : root(nullptr), nodeCapacity(10) {}

        bool built() const {
            return root != nullptr;
        }

        void insert(const ItemType* x) {
            createLeafNode(x);
        }

        template<typename Visitor>
        void query(const geom::Envelope& queryEnv, Visitor&& visitor) {
            if (!built()) {
                build();
            }

            if (root->envelopeIntersects(queryEnv)) {
                query(queryEnv, *root, visitor);
            }
        }

        void insert(const geom::Envelope* itemEnv, void* item) override {
            insert(static_cast<ItemType*>(item));
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
            throw std::runtime_error("Not implemented.");
        }

    private:
        using Node = BastonNode<ItemType>;
        using NodeList = std::vector<Node>;
        using NodeListIterator = typename NodeList::iterator;

        NodeList nodes;
        BastonNode<ItemType>* root;
        size_t nodeCapacity;


        void createLeafNode(const ItemType* item) {
            nodes.emplace_back(item);
        }

        Node* createBranchNode() {
            assert(nodes.size() < nodes.capacity());

            nodes.emplace_back();
            return &nodes.back();
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

        void createParentNodes(NodeListIterator& begin, NodeListIterator& end) {
            auto numChildren = std::distance(begin, end);

            auto numSlices = sliceCount(numChildren);
            auto nodesPerSlice = sliceCapacity(numChildren, numSlices);

            sortNodesX(begin, end);

            auto startOfSlice = begin;
            for (size_t j = 0; j < numSlices; j++) {
                auto nodesRemaining = static_cast<size_t>(std::distance(startOfSlice, end));
                auto nodesInSlice = std::min(nodesRemaining, nodesPerSlice);

                auto endOfSlice = std::next(startOfSlice, nodesInSlice);
                addParentNodesFromVerticalSlice(startOfSlice, endOfSlice);

                startOfSlice = endOfSlice;
            }
        }

        void addParentNodesFromVerticalSlice(NodeListIterator& begin, NodeListIterator& end) {
            sortNodesY(begin, end);

            Node* parent = nullptr;
            for(auto it = begin; it != end; ++it) {
                if (!parent) {
                    parent = createBranchNode();
                }

                parent->addChildNode(&*it);

                if (parent->numChildren() == nodeCapacity) {
                    parent = nullptr;
                }
            }
        }

        void sortNodesX(NodeListIterator& begin, NodeListIterator& end) {
            std::sort(begin, end, [](const Node& a, const Node& b) {
                const geom::Envelope& ea = a.getEnvelope();
                const geom::Envelope& eb = b.getEnvelope();

                double xa = ea.getMinX() + ea.getMaxX();
                double xb = eb.getMinX() + eb.getMaxX();
                return xa < xb;
            });
        }

        void sortNodesY(NodeListIterator& begin, NodeListIterator& end) {
            std::sort(begin, end, [](const Node& a, const Node& b) {
                const geom::Envelope& ea = a.getEnvelope();
                const geom::Envelope& eb = b.getEnvelope();

                double ya = ea.getMinY() + ea.getMaxY();
                double yb = eb.getMinY() + eb.getMaxY();
                return ya < yb;
            });
        }

        template<typename Visitor>
        void query(const geom::Envelope& queryEnv,
                   const Node& node,
                   Visitor&& visitor) {

            for (auto it = node.beginChildren(); it < node.endChildren(); ++it) {
                const Node* child = *it;

                if(!child->envelopeIntersects(queryEnv)) {
                    continue;
                }

                if (child->isLeaf()) {
                    visitor(child->getItem());
                } else {
                    query(queryEnv, *child, visitor);
                }
            }
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
