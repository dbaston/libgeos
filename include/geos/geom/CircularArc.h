/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2024 ISciences, LLC
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#pragma once

#include <geos/export.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/LineSegment.h>
#include <geos/geom/Quadrant.h>
#include <geos/algorithm/CircularArcs.h>
#include <geos/algorithm/Orientation.h>
#include <geos/triangulate/quadedge/TrianglePredicate.h>

namespace geos {
namespace geom {

/// A CircularArc is a reference to three points that define a circular arc.
/// It provides for the lazy calculation of various arc properties such as the center, radius, and orientation
class GEOS_DLL CircularArc {
public:

    using CoordinateXY = geom::CoordinateXY;

    CircularArc() : CircularArc({0, 0}, {0, 0}, {0, 0}) {}

    CircularArc(const CoordinateXY& q0, const CoordinateXY& q1, const CoordinateXY& q2);

    CircularArc(double theta0, double theta2, const CoordinateXY& center, double radius, int orientation);

    CircularArc(const CoordinateXY& q0, const CoordinateXY& q2, const CoordinateXY& center, double radius, int orientation);

    CircularArc(CoordinateSequence&, std::size_t pos);

    CircularArc(CoordinateSequence&, std::size_t pos, const CoordinateXY& center, double radius, int orientation);

    ~CircularArc();

    CircularArc(const CircularArc& other) :
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



    CircularArc(CircularArc&& other) noexcept :
        m_seq(other.m_seq),
        m_pos(other.m_pos),
        m_center(other.m_center),
        m_radius(other.m_radius),
        m_orientation(other.m_orientation),
        m_center_known(other.m_center_known),
        m_radius_known(other.m_radius_known),
        m_orientation_known(other.m_orientation_known),
        m_own_coordinates(other.m_own_coordinates)
    {
        if (m_own_coordinates) {
            other.m_own_coordinates = false;
        }
    }

    CircularArc& operator=(CircularArc&& other) noexcept
    {
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


    /// Return the orientation of the arc as one of:
    /// - algorithm::Orientation::CLOCKWISE,
    /// - algorithm::Orientation::COUNTERCLOCKWISE
    /// - algorithm::Orientation::COLLINEAR
    int getOrientation() const {
        if (!m_orientation_known) {
            m_orientation = algorithm::Orientation::index(p0(), p1(), p2());
            m_orientation_known = true;
        }
        return m_orientation;
    }

    bool isCCW() const {
        return getOrientation() == algorithm::Orientation::COUNTERCLOCKWISE;
    }

    /// Return the center point of the circle associated with this arc
    const CoordinateXY& getCenter() const {
        if (!m_center_known) {
            m_center = algorithm::CircularArcs::getCenter(p0(), p1(), p2());
            m_center_known = true;
        }

        return m_center;
    }

    /// Return the radius of the circle associated with this arc
    double getRadius() const {
        if (!m_radius_known) {
            m_radius = getCenter().distance(p0());
            m_radius_known = true;
        }

        return m_radius;
    }

    /// Return whether this arc forms a complete circle
    bool isCircle() const {
        return p0().equals(p2());
    }

    /// Returns whether this arc forms a straight line (p0, p1, and p2 are collinear)
    bool isLinear() const {
        return !std::isfinite(getRadius());
    }

    /// Return the inner angle of the sector associated with this arc
    double getAngle() const;

    /// Return the length of the arc
    double getLength() const;

    /// Return the area enclosed by the arc p0-p1-p2 and the line segment p2-p0
    double getArea() const;

    /// Return the distance from the centerpoint of the arc to the line segment formed by the end points of the arc.
    double getSagitta() const {
        CoordinateXY midpoint = algorithm::CircularArcs::getMidpoint(p0(), p2(), getCenter(), getRadius(), isCCW());
        return algorithm::Distance::pointToSegment(midpoint, p0(), p2());
    }

    /// Return the angle of p0
    double theta0() const {
        return algorithm::CircularArcs::getAngle(p0(), getCenter());
    }

    /// Return the angle of p2
    double theta2() const {
        return algorithm::CircularArcs::getAngle(p2(), getCenter());
    }

    /// Check to see if a coordinate lies on the arc
    /// Only the angle is checked, so it is assumed that the point lies on
    /// the circle of which this arc is a part.
    bool containsPointOnCircle(const CoordinateXY& q) const {
        double theta = std::atan2(q.y - getCenter().y, q.x - getCenter().x);
        return containsAngle(theta);
    }

    /// Check to see if a coordinate lies on the arc, after testing whether
    /// it lies on the circle.
    bool containsPoint(const CoordinateXY& q) const;

    /// Check to see if a given angle lies on this arc
    bool containsAngle(double theta) const;

    /// Return true if the arc is pointing positive in the y direction
    /// at the location of a specified point. The point is assumed to
    /// be on the arc.
    bool isUpwardAtPoint(const CoordinateXY& q) const;

    void reverse();

    // Split an arc at a specified point.
    // The point is assumed to be on the arc.
    std::pair<CircularArc, CircularArc> splitAtPoint(const CoordinateXY& q) const;

    class Iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = geom::CoordinateXY;
        using pointer = const geom::CoordinateXY*;
        using reference = const geom::CoordinateXY&;

        Iterator(const CircularArc& arc, int i) : m_arc(arc), m_i(i) {}

        reference operator*() const {
            return m_i == 0 ? m_arc.p0() : (m_i == 1 ? m_arc.p1() : m_arc.p2());
        }

        Iterator& operator++() {
            m_i++;
            return *this;
        }

        Iterator operator++(int) {
            Iterator ret = *this;
            m_i++;
            return ret;
        }

        bool operator==(const Iterator& other) const {
            return m_i == other.m_i;
        }

        bool operator!=(const Iterator& other) const {
            return !(*this == other);
        }

    private:
        const CircularArc& m_arc;
        int m_i;

    };

    Iterator begin() const {
        return Iterator(*this, 0);
    }

    Iterator end() const {
        return Iterator(*this, 3);
    }

    const CoordinateXY& p0() const {
        return m_seq->getAt<CoordinateXY>(m_pos);
    }

    const CoordinateXY& p1() const {
        return m_seq->getAt<CoordinateXY>(m_pos + 1);
    }

    const CoordinateXY& p2() const {
        return m_seq->getAt<CoordinateXY>(m_pos + 2);
    }

    std::string toString() const;

    const CoordinateSequence* getCoordinateSequence() const {
        return m_seq;
    }

    std::size_t getCoordinatePosition() const {
        return m_pos;
    }

private:
    CoordinateSequence* m_seq;
    std::size_t m_pos;

    //CoordinateXY m_p0;
    //CoordinateXY m_p1;
    //CoordinateXY m_p2;

    mutable CoordinateXY m_center;
    mutable double m_radius;
    mutable int m_orientation;
    mutable bool m_center_known = false;
    mutable bool m_radius_known = false;
    mutable bool m_orientation_known = false;
    bool m_own_coordinates;
};

}
}
