/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2022 ISciences, LLC
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/


#pragma once

#include <cstddef>
#include <iterator>

namespace geos {
namespace geom {

template<typename CoordinateType>
class CoordinateSequenceIterator {
    CoordinateType* m_ptr;

public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = CoordinateType;
    using reference = CoordinateType&;
    using pointer = CoordinateType;
    using difference_type = std::ptrdiff_t;

    CoordinateSequenceIterator(CoordinateType* c) : m_ptr(c) {}

    reference operator*() const {
        return *m_ptr;
    }

    pointer operator->() const {
        return m_ptr;
    }

    CoordinateSequenceIterator& operator++() {
        m_ptr++;
        return *this;
    }

    CoordinateSequenceIterator operator++(int) {
        CoordinateSequenceIterator ret = *this;
        m_ptr++;
        return ret;
    }

    CoordinateSequenceIterator& operator--() {
        m_ptr--;
        return *this;
    }

    CoordinateSequenceIterator operator--(int) {
        CoordinateSequenceIterator ret = *this;
        m_ptr--;
        return ret;
    }

    difference_type operator-(const CoordinateSequenceIterator& other) const {
        return this->m_ptr - other.m_ptr;
    }

    CoordinateSequenceIterator operator+(difference_type n) const {
        return CoordinateSequenceIterator(this->m_ptr + n);
    }

    CoordinateSequenceIterator operator+=(difference_type n) {
        this->m_ptr += n;
        return *this;
    }

    CoordinateSequenceIterator operator-(difference_type n) const {
        return CoordinateSequenceIterator(this->m_ptr - n);
    }

    CoordinateSequenceIterator operator-=(difference_type n) {
        this->m_ptr -= n;
        return *this;
    }

    CoordinateType& operator[](difference_type n) const {
        return *(*this + n);
    }

    bool operator==(const CoordinateSequenceIterator& other) const {
        return this->m_ptr == other.m_ptr;
    }

    bool operator!=(const CoordinateSequenceIterator& other) const {
        return !(*this == other);
    }

    bool operator<(const CoordinateSequenceIterator& other) const {
        return this->m_ptr < other.m_ptr;
    }

    bool operator<=(const CoordinateSequenceIterator& other) const {
        return this->m_ptr <= other.m_ptr;
    }

    bool operator>(const CoordinateSequenceIterator& other) const {
        return this->m_ptr > other.m_ptr;
    }

    bool operator>=(const CoordinateSequenceIterator& other) const {
        return this->m_ptr >= other.m_ptr;
    }


};

}
}
