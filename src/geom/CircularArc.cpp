/**********************************************************************
*
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2024-2025 ISciences, LLC
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/geom/CircularArc.h>
#include <sstream>

namespace geos::geom {

CircularArc::CircularArc(const CoordinateXY& q0, const CoordinateXY& q1, const CoordinateXY& q2)
: m_seq(new CoordinateSequence(3, false, false))
, m_pos(0)
, m_own_coordinates(true)
{
    m_seq->setAt(q0, 0);
    m_seq->setAt(q1, 1);
    m_seq->setAt(q2, 2);
}

#if 0
CircularArc::CircularArc(const Coordinate& q0, const Coordinate& q1, const Coordinate& q2)
    : m_seq(new CoordinateSequence(3, true, false))
    , m_pos(0)
    , m_own_coordinates(true)
{
    m_seq->setAt(q0, 0);
    m_seq->setAt(q1, 1);
    m_seq->setAt(q2, 2);
}

CircularArc::CircularArc(const CoordinateXYM& q0, const CoordinateXYM& q1, const CoordinateXYM& q2)
    : m_seq(new CoordinateSequence(3, false, true))
    , m_pos(0)
    , m_own_coordinates(true) {
    m_seq->setAt(q0, 0);
    m_seq->setAt(q1, 1);
    m_seq->setAt(q2, 2);
}

CircularArc::CircularArc(const CoordinateXYZM& q0, const CoordinateXYZM& q1, const CoordinateXYZM& q2)
        : m_seq(new CoordinateSequence(3, true, true))
        , m_pos(0)
        , m_own_coordinates(true) {
    m_seq->setAt(q0, 0);
    m_seq->setAt(q1, 1);
    m_seq->setAt(q2, 2);
}
#endif

CircularArc::CircularArc(double theta0, double theta2, const CoordinateXY& center, double radius, int orientation)
        : m_seq(new CoordinateSequence(3, false, false)),
          m_pos(0),
          m_center(center),
          m_radius(radius),
          m_orientation(orientation),
          m_center_known(true),
          m_radius_known(true),
          m_orientation_known(true),
          m_own_coordinates(true)
{
    m_seq->setAt(algorithm::CircularArcs::createPoint(center, radius, theta0), 0);
    m_seq->setAt(algorithm::CircularArcs::createPoint(center, radius, algorithm::CircularArcs::getMidpointAngle(theta0, theta2, orientation==algorithm::Orientation::COUNTERCLOCKWISE)), 1);
    m_seq->setAt(algorithm::CircularArcs::createPoint(center, radius, theta2), 2);
}

CircularArc::CircularArc(const CoordinateXY& q0, const CoordinateXY& q2, const CoordinateXY& center, double radius, int orientation)
    : m_seq(new CoordinateSequence(3, false, false)),
      m_pos(0),
      m_center(center),
      m_radius(radius),
      m_orientation(orientation),
      m_center_known(true),
      m_radius_known(true),
      m_orientation_known(true),
      m_own_coordinates(true)
{
    m_seq->setAt(q0, 0);
    m_seq->setAt(algorithm::CircularArcs::getMidpoint(q0, q2, center, radius, orientation==algorithm::Orientation::COUNTERCLOCKWISE), 1);
    m_seq->setAt(q2, 2);
}

CircularArc::CircularArc(CoordinateSequence& seq, std::size_t pos) :
    m_seq(&seq),
    m_pos(pos),
    m_own_coordinates(false) {}

CircularArc::CircularArc(std::unique_ptr<CoordinateSequence> seq, std::size_t pos) :
    CircularArc(*seq, pos)
{
    m_own_coordinates = true;
    seq.release();
}

CircularArc::CircularArc(CoordinateSequence& seq, std::size_t pos, const CoordinateXY& center, double radius, int orientation) :
    m_seq(&seq),
    m_pos(pos),
    m_center(center),
    m_radius(radius),
    m_orientation(orientation),
    m_center_known(true),
    m_radius_known(true),
    m_orientation_known(true),
    m_own_coordinates(false)
{}

CircularArc::CircularArc(std::unique_ptr<CoordinateSequence> seq, std::size_t pos, const CoordinateXY& center, double radius, int orientation) :
    CircularArc(*seq, pos, center, radius, orientation)
{
    m_own_coordinates = true;
    seq.release();
}

CircularArc::CircularArc(const CircularArc& other) :
    m_seq(new CoordinateSequence(0, other.getCoordinateSequence()->hasZ(), other.getCoordinateSequence()->hasM())),
    m_pos(0),
    m_center(other.m_center),
    m_radius(other.m_radius),
    m_orientation(other.m_orientation),
    m_center_known(other.m_center_known),
    m_radius_known(other.m_radius_known),
    m_orientation_known(other.m_orientation_known),
    m_own_coordinates(true)
{
    m_seq->reserve(3);
    m_seq->add(*other.getCoordinateSequence(), other.getCoordinatePosition(), other.getCoordinatePosition() + 2);
}

CircularArc::CircularArc(CircularArc&& other) noexcept {
    m_seq = other.m_seq;
    m_pos = other.m_pos;
    m_center = other.m_center;
    m_radius = other.m_radius;
    m_orientation = other.m_orientation;
    m_center_known = other.m_center_known;
    m_radius_known = other.m_radius_known;
    m_orientation_known = other.m_orientation_known;
    m_own_coordinates = other.m_own_coordinates;

    if (other.m_own_coordinates) {
        other.m_own_coordinates = false;
    }
}

CircularArc&
CircularArc::operator=(const CircularArc& other)
{
    if (m_own_coordinates) {
        delete m_seq;
    }

    m_seq = new CoordinateSequence(0, other.getCoordinateSequence()->hasZ(), other.getCoordinateSequence()->hasM());
    m_pos = other.m_pos;
    m_own_coordinates = true;
    m_orientation = other.m_orientation;
    m_orientation_known = other.m_orientation_known;
    m_center = other.m_center;
    m_center_known = other.m_center_known;
    m_radius = other.m_radius;
    m_radius_known = other.m_radius_known;

    m_seq->reserve(3);
    m_seq->add(*other.getCoordinateSequence(), other.getCoordinatePosition(), other.getCoordinatePosition() + 2);

    return *this;
}

CircularArc&
CircularArc::operator=(CircularArc&& other) noexcept
{
    if (m_own_coordinates) {
        delete m_seq;
    }

    m_seq = other.m_seq;
    m_pos = other.m_pos;
    m_own_coordinates = other.m_own_coordinates;
    m_orientation = other.m_orientation;
    m_orientation_known = other.m_orientation_known;
    m_center = other.m_center;
    m_center_known = other.m_center_known;
    m_radius = other.m_radius;
    m_radius_known = other.m_radius_known;

    if (m_own_coordinates) {
        other.m_own_coordinates = false;
    }

    return *this;
}

CircularArc::~CircularArc()
{
    if (m_own_coordinates) {
        delete m_seq;
    }
}

bool
CircularArc::containsAngle(double theta) const {
    auto t0 = theta0();
    auto t2 = theta2();

    if (theta == t0 || theta == t2) {
        return true;
    }

    if (getOrientation() == algorithm::Orientation::COUNTERCLOCKWISE) {
        std::swap(t0, t2);
    }

    t2 -= t0;
    theta -= t0;

    if (t2 < 0) {
        t2 += 2*MATH_PI;
    }
    if (theta < 0) {
        theta += 2*MATH_PI;
    }

    return theta >= t2;
}

bool
CircularArc::containsPoint(const CoordinateXY& q) const
{
    if (q == p0() || q == p1() || q == p2()) {
        return true;
    }

    //auto dist = std::abs(q.distance(getCenter()) - getRadius());

    //if (dist > 1e-8) {
    //    return false;
    //}

    if (triangulate::quadedge::TrianglePredicate::isInCircleRobust(p0(), p1(), p2(), q) != geom::Location::BOUNDARY) {
        return false;
    }

    return containsPointOnCircle(q);
}

double
CircularArc::getAngle() const
{
    if (isCircle()) {
        return 2*MATH_PI;
    }

    /// Even Rouault:
    /// potential optimization?: using crossproduct(p0 - center, p2 - center) = radius * radius * sin(angle)
    /// could yield the result by just doing a single asin(), instead of 2 atan2()
    /// actually one should also likely compute dotproduct(p0 - center, p2 - center) = radius * radius * cos(angle),
    /// and thus angle = atan2(crossproduct(p0 - center, p2 - center) , dotproduct(p0 - center, p2 - center) )
    auto t0 = theta0();
    auto t2 = theta2();

    if (getOrientation() == algorithm::Orientation::COUNTERCLOCKWISE) {
        std::swap(t0, t2);
    }

    if (t0 < t2) {
        t0 += 2*MATH_PI;
    }

    auto diff = t0-t2;

    return diff;
}

double
CircularArc::getArea() const {
    if (isLinear()) {
        return 0;
    }

    auto R = getRadius();
    auto theta = getAngle();
    return R*R/2*(theta - std::sin(theta));
}

double
CircularArc::getLength() const {
    if (isLinear()) {
        return p0().distance(p2());
    }

    return getAngle()*getRadius();
}

bool
CircularArc::isUpwardAtPoint(const CoordinateXY& q) const {
    auto quad = geom::Quadrant::quadrant(getCenter(), q);
    bool isUpward;

    if (getOrientation() == algorithm::Orientation::CLOCKWISE) {
        isUpward = (quad == geom::Quadrant::SW || quad == geom::Quadrant::NW);
    } else {
        isUpward = (quad == geom::Quadrant::SE || quad == geom::Quadrant::NE);
    }

    return isUpward;
}

void
CircularArc::reverse() {
    m_seq->swap(m_pos, m_pos+2);
    if (m_orientation_known) {
        if (m_orientation == algorithm::Orientation::COUNTERCLOCKWISE) {
            m_orientation = algorithm::Orientation::CLOCKWISE;
        } else if (m_orientation == algorithm::Orientation::CLOCKWISE) {
            m_orientation = algorithm::Orientation::COUNTERCLOCKWISE;
        }
    }
}

std::pair<CircularArc, CircularArc>
CircularArc::splitAtPoint(const CoordinateXY& q) const {
    return std::make_pair(
        CircularArc(p0(), q, getCenter(), getRadius(), getOrientation()),
        CircularArc(q, p2(), getCenter(), getRadius(), getOrientation())
    );
}

std::string
CircularArc::toString() const {
    std::stringstream ss;
    ss << "CIRCULARSTRING ";
    if (m_seq->hasZ()) {
        ss << "Z";
    }
    if (m_seq->hasM()) {
        ss << "M";
    }
    if (m_seq->hasZ() || m_seq->hasM()) {
        ss << " ";
    }
    ss << "(";
    m_seq->applyAt(m_pos, [&ss](const auto& pt) {
        ss << pt << ", " << *(&pt + 1) << ", " << *(&pt + 2);
    });
    ss << ")";
    return ss.str();
}

}
