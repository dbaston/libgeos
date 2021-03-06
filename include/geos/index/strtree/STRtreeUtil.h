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

#ifndef GEOS_INDEX_STRTREE_STRTREE_UTIL_H
#define GEOS_INDEX_STRTREE_STRTREE_UTIL_H

#include <cstddef>

namespace geos {
namespace index {
namespace strtree {

class STRtreeUtil {

public:
    // calculate what the tree size will be when it is built. This is simply
    // a version of createParentNodes that doesn't actually create anything.
    static std::size_t numBranches(std::size_t numLeafNodes, std::size_t branchSize);

    static std::size_t sliceCount(size_t numNodes, std::size_t branchSize);

    static std::size_t sliceCapacity(size_t numNodes, size_t numSlices);

};

}
}
}

#endif
