/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2020-2021 Daniel Baston
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/index/strtree/STRtreeUtil.h>

#include <algorithm>
#include <cmath>

namespace geos {
namespace index {
namespace strtree {

std::size_t
STRtreeUtil::numBranches(std::size_t numLeafNodes, std::size_t branchSize) {
    size_t nodesInTree = 0;

    size_t nodesWithoutParents = numLeafNodes;
    while (nodesWithoutParents > 1) {
        auto numSlices = sliceCount(nodesWithoutParents, branchSize);
        auto nodesPerSlice = sliceCapacity(nodesWithoutParents, numSlices);

        size_t parentNodesAdded = 0;
        for (size_t j = 0; j < numSlices; j++) {
            auto nodesInSlice = std::min(nodesWithoutParents, nodesPerSlice);
            nodesWithoutParents -= nodesInSlice;

            parentNodesAdded += static_cast<size_t>(std::ceil(
                    static_cast<double>(nodesInSlice) / static_cast<double>(branchSize)));
        }

        nodesInTree += parentNodesAdded;
        nodesWithoutParents = parentNodesAdded;
    }

    return nodesInTree;
}

std::size_t
STRtreeUtil::sliceCount(size_t numNodes, std::size_t branchSize) {
    double minLeafCount = std::ceil(static_cast<double>(numNodes) / static_cast<double>(branchSize));

    return static_cast <size_t> (std::ceil(std::sqrt(minLeafCount)));
}

std::size_t
STRtreeUtil::sliceCapacity(size_t numNodes, size_t numSlices) {
    return static_cast<size_t>(std::ceil(static_cast<double>(numNodes) / static_cast<double>(numSlices)));
}

}
}
}
