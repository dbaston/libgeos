/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/profiler.h>
#include <geos/geom/CoordinateFilter.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Envelope.h>
#include <geos/util/IllegalArgumentException.h>
#include <geos/util.h>

#include <cstdio>
#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>
#include <sstream>

namespace geos {
namespace geom { // geos::geom

#if PROFILE
static Profiler* profiler = Profiler::instance();
#endif

CoordinateSequence::CoordinateSequence(std::size_t size, std::size_t dim) :
    m_vect(size*3),
    m_type(DataType::VECTOR),
    dimension(dim),
    m_stride(3)
{
    assert(dimension == 0 || dimension == 2 || dimension == 3);

    Coordinate defaultCoord;

    for (std::size_t i = 0; i < size; i++) {
        setAt(defaultCoord, i);
    }
}

CoordinateSequence::CoordinateSequence(const std::initializer_list<Coordinate>& list) :
    m_vect(),
    m_type(DataType::VECTOR),
    dimension(0),
    m_stride(3)
{
    reserve(list.size());
    add(list.begin(), list.end());
}

CoordinateSequence::CoordinateSequence(const Coordinate& c) :
    m_coord(c),
    m_type(DataType::SINGLE),
    dimension(std::isnan(c.z) ? 2 : 3),
    m_stride(3)
{
}

CoordinateSequence::CoordinateSequence(double* buf, std::size_t size, std::uint8_t stride, std::size_t dim) :
    m_buf{buf, size},
    m_type(DataType::BUFFER),
    dimension(dim),
    m_stride(stride)
{
}

CoordinateSequence::CoordinateSequence(const CoordinateSequence& other) :
    m_vect(),
    m_type(DataType::VECTOR),
    dimension(other.dimension),
    m_stride(3)
{
    // TODO set this type to single-point if possible

    add(other.cbegin(), other.cend());
}

CoordinateSequence&
CoordinateSequence::operator=(const CoordinateSequence& other) {
    // TODO set this type to single-point if possible

    dimension = other.dimension;
    m_stride = other.m_stride;

    clear();
    add(other.cbegin(), other.cend());
    return *this;
}

CoordinateSequence::CoordinateSequence(CoordinateSequence&& other) :
    m_vect(),
    m_type(DataType::VECTOR),
    dimension(other.dimension),
    m_stride(other.m_stride)
{
    // TODO set this type to single-point if possible

    if (other.m_type == DataType::VECTOR) {
        m_vect = std::move(other.m_vect);
    } else {
        add(other.cbegin(), other.cend());
    }
}

CoordinateSequence&
CoordinateSequence::operator=(CoordinateSequence&& other) {
    // TODO set this type to single-point if possible

    clear();

    if (other.m_type == DataType::VECTOR) {
        m_vect = std::move(other.m_vect);
    } else {
        add(other.cbegin(), other.cend());
    }

    return *this;
}

CoordinateSequence::~CoordinateSequence()
{
    if (m_type == DataType::VECTOR) {
        m_vect.~vector();
    }
}

#if 0
CoordinateSequence::CoordinateSequence(std::vector<Coordinate> && coords, std::size_t dim) :
    vect(std::move(coords)),
    dimension(dim)
{
    assert(dimension == 2 || dimension == 3);
}

CoordinateSequence::CoordinateSequence(std::unique_ptr<std::vector<Coordinate>> && coords, std::size_t dim) :
    vect(std::move(*coords)),
    dimension(dim)
{
}
#endif

void
CoordinateSequence::add(const Coordinate& c, bool allowRepeated)
{
    convertToVector();

    if(!allowRepeated && !isEmpty()) {
        const Coordinate& last = back();
        if(last.equals2D(c)) {
            return;
        }
    }
    add(c);
}

void
CoordinateSequence::add(const CoordinateSequence* cl, bool allowRepeated, bool direction)
{
    convertToVector();
    // FIXME:  don't rely on negative values for 'j' (the reverse case)

    const auto npts = cl->size();
    if(direction) {
        for(std::size_t i = 0; i < npts; ++i) {
            add(cl->getAt(i), allowRepeated);
        }
    }
    else {
        for(auto j = npts; j > 0; --j) {
            add(cl->getAt(j - 1), allowRepeated);
        }
    }
}

/*public*/
void
CoordinateSequence::add(std::size_t i, const Coordinate& coord, bool allowRepeated)
{
    convertToVector();

    // don't add duplicate coordinates
    if(! allowRepeated) {
        std::size_t sz = size();
        if(sz > 0) {
            if(i > 0) {
                const Coordinate& prev = getAt(i - 1);
                if(prev.equals2D(coord)) {
                    return;
                }
            }
            if(i < sz) {
                const Coordinate& next = getAt(i);
                if(next.equals2D(coord)) {
                    return;
                }
            }
        }
    }

    m_vect.insert(
                std::next(m_vect.begin(), static_cast<std::ptrdiff_t>(i * m_stride)),
                &coord.x,
                &coord.x + m_stride);
}

void
CoordinateSequence::apply_rw(const CoordinateFilter* filter)
{
    for(auto it = begin(); it != end(); ++it) {
        filter->filter_rw(&*it);
    }
    dimension = 0; // re-check (see http://trac.osgeo.org/geos/ticket/435)
}

void
CoordinateSequence::apply_ro(CoordinateFilter* filter) const
{
    for(auto it = cbegin(); it != cend(); ++it) {
        filter->filter_ro(&*it);
    }
}

std::unique_ptr<CoordinateSequence>
CoordinateSequence::clone() const
{
    return detail::make_unique<CoordinateSequence>(*this);
}

void
CoordinateSequence::closeRing()
{
    if(!isEmpty() && front() != back()) {
        add(front());
    }
}

std::size_t
CoordinateSequence::getDimension() const
{
    if(dimension != 0) {
        return dimension;
    }

    if(isEmpty()) {
        return 3;
    }

    if(std::isnan(getAt(0).z)) {
        dimension = 2;
    }
    else {
        dimension = 3;
    }

    return dimension;
}


double
CoordinateSequence::getOrdinate(std::size_t index, std::size_t ordinateIndex) const
{
    switch(ordinateIndex) {
        case CoordinateSequence::X:
            return getAt(index).x;
        case CoordinateSequence::Y:
            return getAt(index).y;
        case CoordinateSequence::Z:
            return getAt(index).z;
        default:
            return DoubleNotANumber;
    }
}

bool
CoordinateSequence::hasRepeatedPoints() const
{
    const std::size_t p_size = getSize();
    for(std::size_t i = 1; i < p_size; i++) {
        if(getAt(i - 1) == getAt(i)) {
            return true;
        }
    }
    return false;
}

/*
 * Returns either the given coordinate array if its length is greater than the
 * given amount, or an empty coordinate array.
 */
CoordinateSequence*
CoordinateSequence::atLeastNCoordinatesOrNothing(std::size_t n,
        CoordinateSequence* c)
{
    if(c->getSize() >= n) {
        return c;
    }
    else {
        // FIXME: return NULL rather then empty coordinate array
        return new CoordinateSequence(0, c->getDimension());
    }
}


bool
CoordinateSequence::hasRepeatedPoints(const CoordinateSequence* cl)
{
    const std::size_t size = cl->getSize();
    for(std::size_t i = 1; i < size; i++) {
        if(cl->getAt(i - 1) == cl->getAt(i)) {
            return true;
        }
    }
    return false;
}


const Coordinate*
CoordinateSequence::minCoordinate() const
{
    const Coordinate* minCoord = nullptr;
    const std::size_t p_size = getSize();
    for(std::size_t i = 0; i < p_size; i++) {
        if(minCoord == nullptr || minCoord->compareTo(getAt(i)) > 0) {
            minCoord = &getAt(i);
        }
    }
    return minCoord;
}

size_t
CoordinateSequence::indexOf(const Coordinate* coordinate,
                            const CoordinateSequence* cl)
{
    std::size_t p_size = cl->size();
    for(std::size_t i = 0; i < p_size; ++i) {
        if((*coordinate) == cl->getAt(i)) {
            return i;
        }
    }
    return std::numeric_limits<std::size_t>::max();
}

void
CoordinateSequence::scroll(CoordinateSequence* cl,
                           const Coordinate* firstCoordinate)
{
    // FIXME: use a standard algorithm instead
    std::size_t i, j = 0;
    std::size_t ind = indexOf(firstCoordinate, cl);
    if(ind < 1) {
        return;    // not found or already first
    }

    const std::size_t length = cl->getSize();
    std::vector<Coordinate> v(length);
    for(i = ind; i < length; i++) {
        v[j++] = cl->getAt(i);
    }
    for(i = 0; i < ind; i++) {
        v[j++] = cl->getAt(i);
    }
    cl->setPoints(v);
}

int
CoordinateSequence::increasingDirection(const CoordinateSequence& pts)
{
    std::size_t ptsize = pts.size();
    for(std::size_t i = 0, n = ptsize / 2; i < n; ++i) {
        std::size_t j = ptsize - 1 - i;
        // skip equal points on both ends
        int comp = pts[i].compareTo(pts[j]);
        if(comp != 0) {
            return comp;
        }
    }
    // array must be a palindrome - defined to be in positive direction
    return 1;
}

/* public */
bool
CoordinateSequence::isRing() const
{
    if (size() < 4)
        return false;

    if (getAt(0) != getAt(size()-1))
        return false;

    return true;
}

void
CoordinateSequence::reverse(CoordinateSequence* cl)
{
    // FIXME: use a standard algorithm
    auto last = cl->size() - 1;
    auto mid = last / 2;
    for(std::size_t i = 0; i <= mid; i++) {
        const Coordinate tmp = cl->getAt(i);
        cl->setAt(cl->getAt(last - i), i);
        cl->setAt(tmp, last - i);
    }
}

bool
CoordinateSequence::equals(const CoordinateSequence* cl1,
                           const CoordinateSequence* cl2)
{
    // FIXME: use std::equals()

    if(cl1 == cl2) {
        return true;
    }
    if(cl1 == nullptr || cl2 == nullptr) {
        return false;
    }
    std::size_t npts1 = cl1->getSize();
    if(npts1 != cl2->getSize()) {
        return false;
    }
    for(std::size_t i = 0; i < npts1; i++) {
        if(!(cl1->getAt(i) == cl2->getAt(i))) {
            return false;
        }
    }
    return true;
}

void
CoordinateSequence::expandEnvelope(Envelope& env) const
{
    const std::size_t p_size = getSize();
    for(std::size_t i = 0; i < p_size; i++) {
        env.expandToInclude(getAt(i));
    }
}

Envelope
CoordinateSequence::getEnvelope() const {
    Envelope e;
    expandEnvelope(e);
    return e;
}

CoordinateSequence::iterator
CoordinateSequence::begin() {
    return {this};
}

CoordinateSequence::const_iterator
CoordinateSequence::begin() const {
    return {this};
}

CoordinateSequence::iterator
CoordinateSequence::end() {
    return {this, getSize()};
}

CoordinateSequence::const_iterator
CoordinateSequence::end() const {
    return {this, getSize()};
}

CoordinateSequence::const_iterator
CoordinateSequence::cbegin() const {
    return {this};
}

CoordinateSequence::const_iterator
CoordinateSequence::cend() const {
    return {this, getSize()};
}

void
CoordinateSequence::setOrdinate(std::size_t index, std::size_t ordinateIndex, double value)
{
    switch(ordinateIndex) {
        case CoordinateSequence::X:
        getAt(index).x = value;
        break;
        case CoordinateSequence::Y:
        getAt(index).y = value;
        break;
        case CoordinateSequence::Z:
        getAt(index).z = value;
        break;
        default: {
            std::stringstream ss;
            ss << "Unknown ordinate index " << ordinateIndex;
            throw util::IllegalArgumentException(ss.str());
            break;
        }
    }
}

void
CoordinateSequence::setPoints(const std::vector<Coordinate>& v)
{
    assert(m_stride == 3);
    const double* cbuf = reinterpret_cast<const double*>(v.data());
    m_vect.assign(cbuf, cbuf + v.size()*sizeof(Coordinate)/sizeof(double));
}

void
CoordinateSequence::toVector(std::vector<Coordinate>& out) const
{
    if (m_stride == 3) {
        out.insert(out.end(),
                   reinterpret_cast<const Coordinate*>(data()),
                   reinterpret_cast<const Coordinate*>(data() + m_stride * size()));
    } else {
        out.insert(out.end(), cbegin(), cend());
    }
}

void
CoordinateSequence::pop_back()
{
    convertToVector();

    switch (m_stride) {
    case 4: m_vect.pop_back(); // fall through
    case 3: m_vect.pop_back(); // fall through
    case 2: m_vect.pop_back();
            m_vect.pop_back();
        break;
    default:
        assert(0);
    }
}


std::ostream&
operator<< (std::ostream& os, const CoordinateSequence& cs)
{
    os << "(";
    for(std::size_t i = 0, n = cs.size(); i < n; ++i) {
        const Coordinate& c = cs[i];
        if(i) {
            os << ", ";
        }
        os << c;
    }
    os << ")";

    return os;
}

std::string
CoordinateSequence::toString() const
{
    std::ostringstream ss;
    ss << *this;
    return ss.str();
}

bool
operator== (const CoordinateSequence& s1, const CoordinateSequence& s2)
{
    return CoordinateSequence::equals(&s1, &s2);
}

bool
operator!= (const CoordinateSequence& s1, const CoordinateSequence& s2)
{
    return ! CoordinateSequence::equals(&s1, &s2);
}

void
CoordinateSequence::convertToVector() {
    switch (m_type) {
        case DataType::SINGLE: {
            const double* from = &m_coord.x;
            std::vector<double> vect(m_stride);
            vect.assign(from, from + m_stride);
            new(&m_vect) std::vector<double>(std::move(vect));
            m_type = DataType::VECTOR;
            return;
        }
        case DataType::BUFFER: {
            std::vector<double> vect(m_buf.m_buf_size);
            vect.assign(m_buf.m_buf, m_buf.m_buf + m_buf.m_buf_size);
            new(&m_vect) std::vector<double>(std::move(vect));
            m_type = DataType::VECTOR;
            return;
        }
        case DataType::VECTOR: return;
    }
}

} // namespace geos::geom
} // namespace geos
