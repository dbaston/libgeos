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
 **********************************************************************/

#ifndef GEOS_COORDINATEBUFFERSEQUENCE_H
#define GEOS_COORDINATEBUFFERSEQUENCE_H

#include <geos/geom/CoordinateSequence.h>

namespace geos {
namespace geom{

class CoordinateBufferSequence : public CoordinateSequence {

public:
    CoordinateBufferSequence(void* buff, size_t dim) : m_buff{buff}, m_dim{dim} {}

private:
    void* m_buff;
    size_t m_dim;
};

}
}


#endif //GEOS_COORDINATEBUFFERSEQUENCE_H
