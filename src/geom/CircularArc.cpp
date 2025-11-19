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
    ss << "CIRCULARSTRING (" << p0() << ", " << p1() << ", " << p2() << ")";
    return ss.str();
}

}
