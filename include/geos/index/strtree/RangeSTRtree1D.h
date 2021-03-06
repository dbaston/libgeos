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

#ifndef GEOS_INDEX_STRTREE_RANGESTRTREE1D_H
#define GEOS_INDEX_STRTREE_RANGESTRTREE1D_H

#include <geos/geom/Envelope.h>
#include <geos/index/strtree/STRtreeUtil.h>
#include <geos/index/strtree/Interval.h>
#include <geos/util.h>

#include <array>
#include <vector>

#include <immintrin.h>

namespace geos {
namespace index {
namespace strtree {

template<typename BoundsType, size_t NodeSize>
struct Branch1D {
    constexpr size_t capacity() const {
        return NodeSize;
    }

    template<typename T>
    Branch1D(const T& begin, const T& end) {
        assert(static_cast<size_t>(std::distance(begin, end)) <= NodeSize);

        std::size_t k = 0;
        for(auto it = begin; it != end; ++it) {
            min[k] = detail::nearest_down<BoundsType>(it->getMin());
            max[k] = detail::nearest_up<BoundsType>(it->getMax());
            children[k] = std::addressof(*it);

            minp = std::min(minp, min[k]);
            maxp = std::max(maxp, max[k]);

            k++;
        }
        for (; k < NodeSize; k++) {
            min[k] = max[k] = std::numeric_limits<BoundsType>::quiet_NaN();
        }
    }

    BoundsType getMin() const {
        return minp;
    }

    BoundsType getMax() const {
        return maxp;
    }

    std::array<BoundsType, NodeSize> min;
    std::array<BoundsType, NodeSize> max;

    std::array<void*, NodeSize> children;

    BoundsType minp = std::numeric_limits<BoundsType>::infinity();
    BoundsType maxp = -std::numeric_limits<BoundsType>::infinity();

    //BoundsType mid;
};

template<typename ItemType, typename BoundsType>
class RangeSTRtree1D {
public:
    static const size_t BranchSize = 32u / sizeof(BoundsType);

    void insert(const ItemType& i) {
        leaves.push_back(i);
    }

    void insert(const ItemType& i, const ItemType& j) {
        (void) j;
        leaves.push_back(i);
    }

    bool built() const {
        return !branches.empty();
    }

    void build() {
        // pre-allocate branch vector so we get stable addresses
        branches.reserve(STRtreeUtil::numBranches(leaves.size(), BranchSize));

        sortNodes(leaves.begin(), leaves.end());
        createParentNodes(leaves.begin(), leaves.end());

        // begin and end define a range of nodes needing parents
        auto begin = branches.begin();
        auto end = branches.end();
        while (std::distance(begin, end) > 1) {
            createParentNodes(begin, end);
            begin = end; // parents just added become children in the next round
            end = branches.end();
        }
    }

    template<typename Visitor>
    void query(const index::strtree::Interval& e, Visitor&& visitor) {
        if (!built()) {
            build();
        }

        //if (e.isNull()) {
        //    return;
        //}

        query(branches.back(),
              detail::nearest_down<BoundsType>(e.getMin()),
              detail::nearest_up<BoundsType>(e.getMax()),
              visitor);
    }

    void query(const geom::Envelope* e, std::vector<void*> hits) {
        query(*e, [&hits](const ItemType& item) {
            hits.push_back((void*) item);
        });
    }

private:
    struct Leaf {
        ItemType item;

        Leaf (const ItemType& i) : item(i) {}

        using InternalBoundsType = decltype(item->getMin());

        InternalBoundsType getMin() const {
            return item->getMin();
        }

        InternalBoundsType getMax() const {
            return item->getMax();
        }

        InternalBoundsType getMid() const {
            return getMin() + getMax();
        }
    };

    using BranchT = Branch1D<BoundsType, BranchSize>;

    std::vector<BranchT> branches;
    std::vector<Leaf> leaves;

    using LeafIterator = typename decltype(leaves)::iterator;
    using BranchIterator = typename decltype(branches)::iterator;

    template<typename T>
    void createBranchNode(const T& begin, const T& end) {
        branches.emplace_back(begin, end);
    }

    void sortNodes(const LeafIterator& begin, const LeafIterator& end) {
        std::sort(begin, end, [](const Leaf& a, const Leaf& b) {
            return a.getMid() < b.getMid();
        });
    }

#if 0
    void sortNodes(const BranchIterator& begin, const BranchIterator& end) {
        std::sort(begin, end, [](const BranchT &a, const BranchT &b) {
            return a.mid < b.mid;
        });
    }
#endif

    template<typename NodeIterator>
    void createParentNodes(const NodeIterator& begin, const NodeIterator& end) {
        static constexpr size_t lBranchSize = 32u / sizeof(BoundsType);

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
    void query(const Branch1D<float, 8>& branch, float p_min, float p_max, Visitor&& visitor) {
        auto child_min = _mm256_loadu_ps(branch.min.data());
        auto child_max = _mm256_loadu_ps(branch.max.data());

        auto min = _mm256_set1_ps(p_min);
        auto max = _mm256_set1_ps(p_max);

        auto cmp = _mm256_cmp_ps(child_min, max, _CMP_LE_OQ);
        cmp = _mm256_and_ps(cmp, _mm256_cmp_ps(child_max, min, _CMP_GE_OQ));

        int cmpint = _mm256_movemask_ps(cmp);

        if (isBranch(branch.children[0])) {
            if (cmpint & 1)
                query(*((const BranchT *) branch.children[0]), p_min, p_max, visitor);
            if (cmpint & 2)
                query(*((const BranchT *) branch.children[1]), p_min, p_max, visitor);
            if (cmpint & 4)
                query(*((const BranchT *) branch.children[2]), p_min, p_max, visitor);
            if (cmpint & 8)
                query(*((const BranchT *) branch.children[3]), p_min, p_max, visitor);
            if (cmpint & 16)
                query(*((const BranchT *) branch.children[4]), p_min, p_max, visitor);
            if (cmpint & 32)
                query(*((const BranchT *) branch.children[5]), p_min, p_max, visitor);
            if (cmpint & 64)
                query(*((const BranchT *) branch.children[6]), p_min, p_max, visitor);
            if (cmpint & 128)
                query(*((const BranchT *) branch.children[7]), p_min, p_max, visitor);
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
    void query(const Branch1D<double, 4>& branch, double p_min, double p_max, Visitor&& visitor) {
        auto child_min = _mm256_loadu_pd(branch.min.data());
        auto child_max = _mm256_loadu_pd(branch.max.data());

        auto min = _mm256_set1_pd(p_min);
        auto max = _mm256_set1_pd(p_max);

        auto cmp = _mm256_cmp_pd(child_min, max, _CMP_LE_OQ);
        cmp = _mm256_and_pd(cmp, _mm256_cmp_pd(child_max, min, _CMP_GE_OQ));

        int cmpint = _mm256_movemask_pd(cmp);

        if (isBranch(branch.children[0])) {
            if (cmpint & 1)
                query(*((const BranchT *) branch.children[0]), p_min, p_max, visitor);
            if (cmpint & 2)
                query(*((const BranchT *) branch.children[1]), p_min, p_max, visitor);
            if (cmpint & 4)
                query(*((const BranchT *) branch.children[2]), p_min, p_max, visitor);
            if (cmpint & 8)
                query(*((const BranchT *) branch.children[3]), p_min, p_max, visitor);
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
