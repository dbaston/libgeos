/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2021 Daniel Baston
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_INDEX_STRTREE_RANGESTRTREE_H
#define GEOS_INDEX_STRTREE_RANGESTRTREE_H

#include <geos/geom/Envelope.h>
#include <geos/index/strtree/STRtreeUtil.h>
#include <geos/index/SpatialIndex.h>
#include <geos/index/ItemVisitor.h>
#include <geos/util.h>

// For adapters...move elsewhere.
#include <geos/index/chain/MonotoneChain.h>

#include <array>
#include <vector>

#include <immintrin.h>

namespace geos {
namespace index {
namespace strtree {

template<typename BoundsType, size_t NodeSize>
struct Branch {

    constexpr size_t capacity() const {
        return NodeSize;
    }

    template<typename T>
    Branch(const T& begin, const T& end) {
        assert(static_cast<size_t>(std::distance(begin, end)) <= NodeSize);

        std::size_t k = 0;
        for(auto it = begin; it != end; ++it) {
            xmin[k] = detail::nearest_down<BoundsType>(it->getXMin());
            xmax[k] = detail::nearest_up<BoundsType>(it->getXMax());
            ymin[k] = detail::nearest_down<BoundsType>(it->getYMin());
            ymax[k] = detail::nearest_up<BoundsType>(it->getYMax());
            children[k] = std::addressof(*it);

            xminp = std::min(xminp, xmin[k]);
            xmaxp = std::max(xmaxp, xmax[k]);
            yminp = std::min(yminp, ymin[k]);
            ymaxp = std::max(ymaxp, ymax[k]);

            k++;
        }
        for (; k < NodeSize; k++) {
            xmin[k] = xmax[k] = ymin[k] = ymax[k] = std::numeric_limits<BoundsType>::quiet_NaN();
        }

        midx = getXMin() + getXMax();
        midy = getYMin() + getYMax();
    }

    BoundsType getXMin() const {
        return xminp;
    }

    BoundsType getYMin() const {
        return yminp;
    }

    BoundsType getXMax() const {
        return xmaxp;
    }

    BoundsType getYMax() const {
        return ymaxp;
    }

    std::array<BoundsType, NodeSize> xmin;
    std::array<BoundsType, NodeSize> xmax;
    std::array<BoundsType, NodeSize> ymin;
    std::array<BoundsType, NodeSize> ymax;

    std::array<void*, NodeSize> children;

    BoundsType xminp = std::numeric_limits<BoundsType>::infinity();
    BoundsType xmaxp = -std::numeric_limits<BoundsType>::infinity();
    BoundsType yminp = std::numeric_limits<BoundsType>::infinity();
    BoundsType ymaxp = -std::numeric_limits<BoundsType>::infinity();

    BoundsType midx;
    BoundsType midy;
};

template<typename T>
class EnvelopeAdapter {};

template<>
struct EnvelopeAdapter<const geom::Envelope*> {
    using InternalBoundsType = double;

    static InternalBoundsType getXMin(const geom::Envelope* item) {
        return item->getMinX();
    }

    static InternalBoundsType getXMax(const geom::Envelope* item) {
        return item->getMaxX();
    }

    static InternalBoundsType getYMin(const geom::Envelope* item) {
        return item->getMinY();
    }

    static InternalBoundsType getYMax(const geom::Envelope* item) {
        return item->getMaxY();
    }
};


template<>
struct EnvelopeAdapter<const index::chain::MonotoneChain*> {
    using InternalBoundsType = double;

    static InternalBoundsType getXMin(const index::chain::MonotoneChain* item) {
        return item->getEnvelope().getMinX();
    }

    static InternalBoundsType getXMax(const index::chain::MonotoneChain* item) {
        return item->getEnvelope().getMaxX();
    }

    static InternalBoundsType getYMin(const index::chain::MonotoneChain* item) {
        return item->getEnvelope().getMinY();
    }

    static InternalBoundsType getYMax(const index::chain::MonotoneChain* item) {
        return item->getEnvelope().getMaxY();
    }
};


template<typename ItemType, typename BoundsType>
class RangeSTRtree : public SpatialIndex {
public:
    void insert(const ItemType& i) {
        leaves.emplace_back(i);
    }

    void insert(const geom::Envelope& e, const ItemType& i) {
        // TODO check that e is redundant with i
        (void) e;
        insert(i);
    }

    void insert(const ItemType& i, const ItemType& j) {
        (void) j;
        leaves.emplace_back(i);
    }

    bool built() const {
        return !branches.empty();
    }

    void build() {
        // pre-allocate branch vector so we get stable addresses
        branches.reserve(STRtreeUtil::numBranches(leaves.size(), BranchSize));

        // begin and end define a range of nodes needing parents
        createParentNodes(leaves.begin(), leaves.end());

        auto begin = branches.begin();
        auto end = branches.end();
        // TODO avoid iterator invalidation
        while (std::distance(begin, end) > 1) {
            createParentNodes(begin, end);
            begin = end; // parents just added become children in the next round
            end = branches.end();
        }
    }

    template<typename Visitor>
    void query(const geom::Envelope& e, Visitor&& visitor) {
        if (!built()) {
            build();
        }

        if (e.isNull()) {
            return;
        }

        query(branches.back(),
              detail::nearest_down<BoundsType>(e.getMinX()),
              detail::nearest_up<BoundsType>(e.getMaxX()),
              detail::nearest_down<BoundsType>(e.getMinY()),
              detail::nearest_up<BoundsType>(e.getMaxY()),
              visitor);
    }

    void query(const geom::Envelope* queryEnv, ItemVisitor& visitor) override {
        query(*queryEnv, [&visitor](const ItemType& x) {
            visitor.visitItem(const_cast<void*>(static_cast<const void*>(x)));
        });
    }

    void query(const geom::Envelope* e, std::vector<void*>& hits) override {
        query(*e, [&hits](const ItemType& item) {
            hits.push_back((void*) item);
        });
    }

    void insert(const geom::Envelope* itemEnv, void* item) override {
        static_assert(std::is_pointer<ItemType>::value, "Only available for pointer types.");

        insert(itemEnv, std::move(static_cast<ItemType*>(item)));
    }

private:
    static const size_t BranchSize = 32u / sizeof(BoundsType);

    struct Leaf {
        ItemType item;

        explicit Leaf (const ItemType& i) : item(i) {}

        using InternalBoundsType = typename EnvelopeAdapter<ItemType>::InternalBoundsType;

        InternalBoundsType getXMin() const {
            return EnvelopeAdapter<ItemType>::getXMin(item);
        }

        InternalBoundsType getXMax() const {
            return EnvelopeAdapter<ItemType>::getXMax(item);
        }

        InternalBoundsType getYMin() const {
            return EnvelopeAdapter<ItemType>::getYMin(item);
        }

        InternalBoundsType getYMax() const {
            return EnvelopeAdapter<ItemType>::getYMax(item);
        }

        InternalBoundsType getX() const {
            return getXMin() + getXMax();
        }

        InternalBoundsType getY() const {
            return getYMin() + getYMax();
        }
    };

    using BranchT = Branch<BoundsType, BranchSize>;

    std::vector<BranchT> branches;
    std::vector<Leaf> leaves;

    using LeafIterator = typename decltype(leaves)::iterator;
    using BranchIterator = typename decltype(branches)::iterator;

    template<typename T>
    void createBranchNode(const T& begin, const T& end) {
        branches.emplace_back(begin, end);
    }

    void sortNodesX(const LeafIterator& begin, const LeafIterator& end) {
        std::sort(begin, end, [](const Leaf& a, const Leaf& b) {
            return a.getX() < b.getX();
        });
    }

    void sortNodesY(const LeafIterator& begin, const LeafIterator& end) {
        std::sort(begin, end, [](const Leaf& a, const Leaf& b) {
            return a.getY() < b.getY();
        });
    }

    void sortNodesX(const BranchIterator& begin, const BranchIterator& end) {
        std::sort(begin, end, [](const BranchT &a, const BranchT &b) {
            return a.midx < b.midx;
        });
    }

    void sortNodesY(const BranchIterator& begin, const BranchIterator& end) {
        std::sort(begin, end, [](const BranchT &a, const BranchT &b) {
            return a.midy < b.midy;
        });
    }

    template<typename NodeIterator>
    void createParentNodes(const NodeIterator& begin, const NodeIterator& end) {
        // Arrange child nodes in two dimensions.
        // First, divide them into vertical slices of a given size (left-to-right)
        // Then create nodes within those slices (bottom-to-top)
        auto numChildren = static_cast<size_t>(std::distance(begin, end));
        auto numSlices = STRtreeUtil::sliceCount(numChildren, BranchSize);
        auto nodesPerSlice = STRtreeUtil::sliceCapacity(numChildren, numSlices);

        sortNodesX(begin, end);

        auto startOfSlice = begin;
        for (size_t j = 0; j < numSlices; j++) {
            auto nodesRemaining = static_cast<size_t>(std::distance(startOfSlice, end));
            auto nodesInSlice = std::min(nodesRemaining, nodesPerSlice);
            auto endOfSlice = std::next(startOfSlice, static_cast<long>(nodesInSlice));

            addParentNodesFromVerticalSlice(startOfSlice, endOfSlice);

            startOfSlice = endOfSlice;
        }
    }

    template<typename NodeIterator>
    void addParentNodesFromVerticalSlice(const NodeIterator& begin, const NodeIterator& end) {
        static constexpr size_t lBranchSize = 32u / sizeof(BoundsType);

        sortNodesY(begin, end);

        // Arrange the nodes vertically and fill up parent nodes sequentially until they're full.
        // A possible improvement would be to rework this such so that if we have 81 nodes we
        // put 9 into each parent instead of 10 or 1.
        auto firstChild = begin;
        while (firstChild != end) {
            auto childrenRemaining = static_cast<size_t>(std::distance(firstChild, end));
            auto childrenForNode = std::min(lBranchSize, childrenRemaining);
            auto lastChild = std::next(firstChild, static_cast<long>(childrenForNode));

            createBranchNode(firstChild, lastChild);
            firstChild = lastChild;
        }

    }

    template<typename Visitor>
    void query(const Branch<float, 8>& branch, float p_xmin, float p_xmax, float p_ymin, float p_ymax, Visitor&& visitor) {
        auto child_xmin = _mm256_loadu_ps(branch.xmin.data());
        auto child_xmax = _mm256_loadu_ps(branch.xmax.data());
        auto child_ymin = _mm256_loadu_ps(branch.ymin.data());
        auto child_ymax = _mm256_loadu_ps(branch.ymax.data());

        auto xmin = _mm256_set1_ps(p_xmin);
        auto xmax = _mm256_set1_ps(p_xmax);
        auto ymin = _mm256_set1_ps(p_ymin);
        auto ymax = _mm256_set1_ps(p_ymax);

        auto cmp = _mm256_cmp_ps(child_xmin, xmax, _CMP_LE_OQ);
        cmp = _mm256_and_ps(cmp, _mm256_cmp_ps(child_xmax, xmin, _CMP_GE_OQ));
        cmp = _mm256_and_ps(cmp, _mm256_cmp_ps(child_ymin, ymax, _CMP_LE_OQ));
        cmp = _mm256_and_ps(cmp, _mm256_cmp_ps(child_ymax, ymin, _CMP_GE_OQ));

        int cmpint = _mm256_movemask_ps(cmp);

        if (isBranch(branch.children[0])) {
            if (cmpint & 1)
                query(*((const BranchT *) branch.children[0]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 2)
                query(*((const BranchT *) branch.children[1]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 4)
                query(*((const BranchT *) branch.children[2]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 8)
                query(*((const BranchT *) branch.children[3]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 16)
                query(*((const BranchT *) branch.children[4]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 32)
                query(*((const BranchT *) branch.children[5]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 64)
                query(*((const BranchT *) branch.children[6]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 128)
                query(*((const BranchT *) branch.children[7]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
        } else {
            if (cmpint & 1)
                visitor(*((const ItemType *) branch.children[0]));
            if (cmpint & 2)
                visitor(*((const ItemType *) branch.children[1]));
            if (cmpint & 4)
                visitor(*((const ItemType *) branch.children[2]));
            if (cmpint & 8)
                visitor(*((const ItemType *) branch.children[3]));
            if (cmpint & 16)
                visitor(*((const ItemType *) branch.children[4]));
            if (cmpint & 32)
                visitor(*((const ItemType *) branch.children[5]));
            if (cmpint & 64)
                visitor(*((const ItemType *) branch.children[6]));
            if (cmpint & 128)
                visitor(*((const ItemType *) branch.children[7]));
        }
    }

    template<typename Visitor>
    void query(const Branch<double, 4>& branch, double p_xmin, double p_xmax, double p_ymin, double p_ymax, Visitor&& visitor) {
        auto child_xmin = _mm256_loadu_pd(branch.xmin.data());
        auto child_xmax = _mm256_loadu_pd(branch.xmax.data());
        auto child_ymin = _mm256_loadu_pd(branch.ymin.data());
        auto child_ymax = _mm256_loadu_pd(branch.ymax.data());

        auto xmin = _mm256_set1_pd(p_xmin);
        auto xmax = _mm256_set1_pd(p_xmax);
        auto ymin = _mm256_set1_pd(p_ymin);
        auto ymax = _mm256_set1_pd(p_ymax);

        auto cmp = _mm256_cmp_pd(child_xmin, xmax, _CMP_LE_OQ);
        cmp = _mm256_and_pd(cmp, _mm256_cmp_pd(child_xmax, xmin, _CMP_GE_OQ));
        cmp = _mm256_and_pd(cmp, _mm256_cmp_pd(child_ymin, ymax, _CMP_LE_OQ));
        cmp = _mm256_and_pd(cmp, _mm256_cmp_pd(child_ymax, ymin, _CMP_GE_OQ));

        int cmpint = _mm256_movemask_pd(cmp);

        if (isBranch(branch.children[0])) {
            if (cmpint & 1)
                query(*((const BranchT *) branch.children[0]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 2)
                query(*((const BranchT *) branch.children[1]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 4)
                query(*((const BranchT *) branch.children[2]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
            if (cmpint & 8)
                query(*((const BranchT *) branch.children[3]), p_xmin, p_xmax, p_ymin, p_ymax, visitor);
        } else {
            if (cmpint & 1)
                visitor(*((const ItemType *) branch.children[0]));
            if (cmpint & 2)
                visitor(*((const ItemType *) branch.children[1]));
            if (cmpint & 4)
                visitor(*((const ItemType *) branch.children[2]));
            if (cmpint & 8)
                visitor(*((const ItemType *) branch.children[3]));
        }
    }

    bool isBranch(void* p) {
        return (p >= branches.data()) && (p < (branches.data() + branches.size()));
    }

};

}
}
}

#endif //GEOS_RANGESTRTREE_H
