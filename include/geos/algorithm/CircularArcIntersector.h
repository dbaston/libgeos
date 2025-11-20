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

#include <array>
#include <cstdint>

#include <geos/export.h>
#include <geos/algorithm/LineIntersector.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CircularArc.h>
#include <geos/geom/LineSegment.h>

namespace geos::algorithm {

class GEOS_DLL CircularArcIntersector {
public:
    using CoordinateXY = geom::CoordinateXY;
    using CoordinateXYZM = geom::CoordinateXYZM;
    using CircularArc = geom::CircularArc;
    using Envelope = geom::Envelope;

    enum intersection_type : uint8_t {
        NO_INTERSECTION = 0,
        ONE_POINT_INTERSECTION = 1,
        TWO_POINT_INTERSECTION = 2,
        COCIRCULAR_INTERSECTION = 3,
    };

    intersection_type getResult() const
    {
        return result;
    }

    const CoordinateXYZM& getPoint(std::uint8_t i) const
    {
        return intPt[i];
    }

    const CircularArc& getArc(std::uint8_t i) const
    {
        return intArc[i];
    }

    std::uint8_t getNumPoints() const
    {
        return nPt;
    }

    std::uint8_t getNumArcs() const
    {
        return nArc;
    }

#if 0
    void computeIntersection(const geom::CoordinateSequence& p, std::size_t p0, const CoordinateXY& centerP, double radiusP,
                             const geom::CoordinateSequence& q, std::size_t q0, const CoordinateXY& centerQ, double radiusQ);
#endif

    /// Determines whether and where a circular arc intersects a line segment.
    ///
    /// Sets the appropriate value of intersection_type and stores the intersection
    /// points, if any.
    //void intersects(const CircularArc& arc, const CoordinateXY& p0, const CoordinateXY& p1);
    void intersects(const CircularArc& arc, const geom::CoordinateSequence& seq, std::size_t pos0, std::size_t pos1);

#if 0
    template<typename C1, typename C2>
    void computeIntersectionArcArc(const CircularArc& arc0, const CircularArc& arc1);
    //void computeIntersectionArcArc(const C1& p0, const C1& p2, const CoordinateXY& centerP, double radiusP, int orientationP,
    //                               const C2& q0, const C2& q2, const CoordinateXY& centerQ, double radiusQ, int orientationQ);

    template<typename C1, typename C2>
    void computeIntersectArcSegment(const C1& p0, const C1& p1, const C1& p2, const CoordinateXY& center, double radius,
                                    const C2& q0, const C2& q1);
#endif

    //void intersects(const CircularArc& arc, const geom::LineSegment& seg)
    //{
    //    intersects(arc, seg.p0, seg.p1);
    //}

    /// Determines whether and where two circular arcs intersect.
    ///
    ///	Sets the appropriate value of intersection_type and stores the intersection
    /// points and/or arcs, if any.
    void intersects(const CircularArc& arc1, const CircularArc& arc2);

    static int
    circleIntersects(const CoordinateXY& center, double r, const CoordinateXY& p0, const CoordinateXY& p1, CoordinateXY& isect0, CoordinateXY& isect1);


    template<typename C1, typename C2>
    void intersects(const C1& p0, const C1& p1, const C2& q0, const C2& q1)
    {
        reset();

        LineIntersector li;
        li.computeIntersection(p0, p1, q0, q1);
        if (li.getIntersectionNum() == 2) {
            // FIXME this means a collinear intersection, so we should report as cocircular?
            intPt[0] = li.getIntersection(0);
            intPt[1] = li.getIntersection(1);
            result = TWO_POINT_INTERSECTION;
        } else if (li.getIntersectionNum() == 1) {
            intPt[0] = li.getIntersection(0);
            nPt = 1;
            result = ONE_POINT_INTERSECTION;
        } else {
            result = NO_INTERSECTION;
        }
    }

private:
    void reset() {
        nPt = 0;
        nArc = 0;
    }

    std::array<CoordinateXYZM, 2> intPt;
    std::array<CircularArc, 2> intArc;
    intersection_type result = NO_INTERSECTION;
    std::uint8_t nPt = 0;
    std::uint8_t nArc = 0;

};

//template CircularArcIntersector::computeIntersectionArcArc<geom::CoordinateXY, geom::CoordinateXY>(const geom::CoordinateXY& p,);

}
