/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2016 Daniel Baston
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: operation/distance/FacetSequence.java (aff11591)
 *
 **********************************************************************/

#include <geos/geom/CircularArc.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/LineSegment.h>
#include <geos/algorithm/CircularArcs.h>
#include <geos/algorithm/Distance.h>
#include <geos/operation/distance/FacetSequence.h>
#include <geos/util/IllegalArgumentException.h>

#include <memory>

using namespace geos::geom;
using namespace geos::operation::distance;
using namespace geos::algorithm;


FacetSequence::FacetSequence(const CoordinateSequence* p_pts, bool p_isCurved, std::size_t p_start, std::size_t p_end)
    : pts(p_pts)
    , start(p_start)
    , end(p_end)
    , curved(p_isCurved)
{
    if (isCurved() && (size() % 2 == 0)) {
        throw util::IllegalArgumentException("Curved FacetSequence must have an odd number of vertices."); }
    computeEnvelope();
}

size_t
FacetSequence::size() const
{
    return end - start;
}

bool
FacetSequence::isCurved() const
{
    return curved;
}

bool
FacetSequence::isPoint() const
{
    return end - start == 1;
}

double
FacetSequence::distance(const FacetSequence& facetSeq) const
{
    return nearestLocations(facetSeq, nullptr);
}

double
FacetSequence::nearestLocations(const FacetSequence& facetSeq, std::vector<CoordinateXY>* locs) const
{
    const bool isCurvedThis = isCurved();
    const bool isCurvedOther = facetSeq.isCurved();
    const bool isPointThis = isPoint();
    const bool isPointOther = facetSeq.isPoint();

    if (isPointThis && isPointOther) {
        const CoordinateXY& pt = pts->getAt<CoordinateXY>(start);
        const CoordinateXY& seqPt = facetSeq.pts->getAt<CoordinateXY>(facetSeq.start);
        return pt.distance(seqPt);
    }
    else if (isPointThis) {
        const CoordinateXY& pt = pts->getAt<CoordinateXY>(start);
        if (isCurvedOther) {
            return computeDistancePointArc(pt, facetSeq, locs);
        } else {
            return computeDistancePointLine(pt, facetSeq, locs);
        }
    }
    else if (isPointOther) {
        const CoordinateXY& seqPt = facetSeq.pts->getAt<CoordinateXY>(facetSeq.start);
        if (isCurvedThis) {
            double d = computeDistancePointArc(seqPt, *this, locs);
            if (locs) {
                std::reverse(locs->begin(), locs->end());
            }
            return d;
        } else {
            double d = computeDistancePointLine(seqPt, *this, locs);
            if (locs) {
                std::reverse(locs->begin(), locs->end());
            }
            return d;
        }
    }
    else if (isCurvedThis) {
        if (isCurvedOther) {
            return computeDistanceArcArc(facetSeq, locs);
        } else {
            return computeDistanceArcLine(facetSeq, locs);
        }
    } else if (isCurvedOther) {
        double d = facetSeq.computeDistanceArcLine(*this, locs);
        if (locs) {
            std::reverse(locs->begin(), locs->end());
        }
        return d;
    } else {
        return computeDistanceLineLine(facetSeq, locs);
    }
}

/*
* Rather than get bent out of shape about returning a pointer
* just return the whole mess, since it only ends up holding two
* locations.
*/
std::vector<CoordinateXY>
FacetSequence::nearestLocations(const FacetSequence& facetSeq) const
{
    std::vector<CoordinateXY> locs;
    nearestLocations(facetSeq, &locs);
    return locs;
}

double
FacetSequence::computeDistancePointLine(const CoordinateXY& pt,
                                        const FacetSequence& facetSeq,
                                        std::vector<CoordinateXY> *locs)
{
    double minDistance = DoubleInfinity;

    for (std::size_t i = facetSeq.start; i < facetSeq.end - 1; i++) {
        const CoordinateXY& q0 = facetSeq.pts->getAt<CoordinateXY>(i);
        const CoordinateXY& q1 = facetSeq.pts->getAt<CoordinateXY>(i + 1);
        double dist = Distance::pointToSegment(pt, q0, q1);
        if (dist < minDistance || (locs != nullptr && locs->empty())) {
            minDistance = dist;
            if (locs != nullptr) {
                updateNearestLocationsPointLine(pt, q0, q1, locs);
            }
            if (minDistance <= 0.0) {
                return minDistance;
            }
        }
    }

    return minDistance;
}

double
FacetSequence::computeDistancePointArc(const CoordinateXY& pt,
                                       const FacetSequence& facetSeq,
                                       std::vector<CoordinateXY>* locs)
{
    double minDistance = DoubleInfinity;

    for (std::size_t i = facetSeq.start; i < facetSeq.end - 2; i+= 2) {
        CircularArc arc(*facetSeq.pts, i);
        CoordinateXY closestPt = arc.closestPoint(pt);
        const double dist = pt.distance(closestPt);

        if (dist < minDistance) {
            minDistance = dist;
            updateNearestLocations(pt, closestPt, locs);

            if (minDistance <= 0.0) {
                return minDistance;
            }
        }
    }

    return minDistance;
}

double
FacetSequence::computeDistanceArcLine(const FacetSequence& facetSeq, std::vector<CoordinateXY>* locs) const
{
    double minDistance = DoubleInfinity;
    const Envelope* lineEnv1 = facetSeq.getEnvelope();

    // brute force approach!
    for (std::size_t i = start; i < end - 2; i+= 2) {
        CircularArc arc(*pts, i);
        const Envelope arcEnv = arc.getEnvelope();

        if (arcEnv.distanceSquared(*lineEnv1) > minDistance*minDistance) {
            continue;
        }

        for(std::size_t j = facetSeq.start; j < facetSeq.end - 1; j++) {
            const CoordinateXY& q0 = facetSeq.pts->getAt<CoordinateXY>(j);
            const CoordinateXY& q1 = facetSeq.pts->getAt<CoordinateXY>(j+1);

            Envelope segEnv(q0, q1);

            if (arcEnv.distanceSquared(segEnv) > minDistance*minDistance) {
                continue;
            }

            const auto closestPts = arc.closestPoints(q0, q1);
            const double dist = closestPts[0].distance(closestPts[1]);
            if(dist < minDistance) {
                minDistance = dist;
                updateNearestLocations(closestPts[0], closestPts[1], locs);
            }

            if(minDistance <= 0) {
                return minDistance;
            }
        }
    }

    return minDistance;
}

double
FacetSequence::computeDistanceArcArc(const FacetSequence& facetSeq, std::vector<CoordinateXY>* locs) const
{
    double minDistance = DoubleInfinity;

    for(std::size_t i = start; i < end - 2; i += 2) {
        CircularArc arc0(*pts, i);

        Envelope pEnv = arc0.getEnvelope();
        if (pEnv.distanceSquared(*facetSeq.getEnvelope()) > minDistance*minDistance) {
            continue;
        }

        for(std::size_t j = facetSeq.start; j < facetSeq.end - 2; j += 2) {
            CircularArc arc1(*facetSeq.pts, j);

            Envelope qEnv = arc1.getEnvelope();
            if (pEnv.distanceSquared(qEnv) > minDistance*minDistance) {
                continue;
            }

            const auto closestPoints = arc0.closestPoints(arc1);
            const double dist = closestPoints[0].distance(closestPoints[1]);

            if(dist <= minDistance) {
                minDistance = dist;
                updateNearestLocations(closestPoints[0], closestPoints[1], locs);

                if(minDistance <= 0.0) {
                    return minDistance;
                }
            }
        }
    }

    return minDistance;

}

void
FacetSequence::updateNearestLocations(const geom::CoordinateXY& p0, const geom::CoordinateXY& p1, std::vector<geom::CoordinateXY> *locs)
{
    if (!locs) {
        return;
    }

    locs->clear();
    locs->push_back(p0);
    locs->push_back(p1);
}

void
FacetSequence::updateNearestLocationsPointLine(
    const CoordinateXY& pt, const CoordinateXY& q0, const CoordinateXY& q1,
    std::vector<CoordinateXY> *locs)
{
    CoordinateXY segClosestPoint = LineSegment::closestPoint(q0, q1, pt);
    updateNearestLocations(pt, segClosestPoint, locs);
}

double
FacetSequence::computeDistanceLineLine(const FacetSequence& facetSeq, std::vector<CoordinateXY>* locs) const
{
    double minDistance = DoubleInfinity;

    for(std::size_t i = start; i < end - 1; i++) {
        const CoordinateXY& p0 = pts->getAt<CoordinateXY>(i);
        const CoordinateXY& p1 = pts->getAt<CoordinateXY>(i + 1);

        // Avoid calculating distance from zero-length segment
        if (p0 == p1)
            continue;

        Envelope pEnv(p0, p1);
        if (pEnv.distanceSquared(*facetSeq.getEnvelope()) > minDistance*minDistance) {
            continue;
        }

        for(std::size_t j = facetSeq.start; j < facetSeq.end - 1; j++) {
            const CoordinateXY& q0 = facetSeq.pts->getAt<CoordinateXY>(j);
            const CoordinateXY& q1 = facetSeq.pts->getAt<CoordinateXY>(j + 1);

            // Avoid calculating distance to zero-length segment
            if (q0 == q1)
                continue;

            Envelope qEnv(q0, q1);
            if (pEnv.distanceSquared(qEnv) > minDistance*minDistance) {
                continue;
            }

            double dist = Distance::segmentToSegment(p0, p1, q0, q1);
            if(dist <= minDistance) {
                minDistance = dist;
                if(locs != nullptr) {
                    updateNearestLocationsLineLine(p0, p1, q0, q1, locs);
                }
                if(minDistance <= 0.0) return minDistance;
            }
        }
    }

    return minDistance;
}

void
FacetSequence::updateNearestLocationsLineLine(
    const CoordinateXY& p0, const CoordinateXY& p1,
    const CoordinateXY& q0, const CoordinateXY& q1,
    std::vector<CoordinateXY> *locs)
{
    const auto closestPts = LineSegment::closestPoints(p0, p1, q0, q1);
    updateNearestLocations(closestPts[0], closestPts[1], locs);
}

void
FacetSequence::computeEnvelope()
{
    env = Envelope();
    if (isCurved()) {
        for (std::size_t i = start; i+2 < end; i+= 2) {
            algorithm::CircularArcs::expandEnvelope(
                env,
                pts->getAt<CoordinateXY>(i),
                pts->getAt<CoordinateXY>(i+1),
                pts->getAt<CoordinateXY>(i+2));
        }
    } else {
        for(std::size_t i = start; i < end; i++) {
            env.expandToInclude(pts->getAt<CoordinateXY>(i));
        }
    }
}

const Envelope*
FacetSequence::getEnvelope() const
{
    return &env;
}

const CoordinateXY*
FacetSequence::getCoordinate(std::size_t index) const
{
    return &(pts->getAt<CoordinateXY>(start + index));
}

std::size_t
FacetSequence::normalize(const CoordinateSequence& pts, std::size_t index)
{
    if (index >= pts.size() - 1 && pts.isRing()) {
        return 0;
    }
    return index;
}

CoordinateSequenceLocation
FacetSequence::nearestLocation(const CoordinateXY& p) const
{
    if (isPoint()) {
        // JTS uses index 0 / pts[0]; use start so a one-point facet
        // that is a subsequence of a longer CoordinateSequence is correct.
        return CoordinateSequenceLocation(pts, start, pts->getAt(start));
    }
    return nearestLocationOnLine(p);
}

CoordinateSequenceLocation
FacetSequence::nearestLocationOnLine(const CoordinateXY& queryPt) const
{
    double minDistance = DoubleInfinity;
    std::size_t index = start;

    CoordinateXY nearestPt;
    for (std::size_t i = start; i < end - 1; i++) {
        const CoordinateXY& q0 = pts->getAt<CoordinateXY>(i);
        const CoordinateXY& q1 = pts->getAt<CoordinateXY>(i + 1);
        double dist = Distance::pointToSegment(queryPt, q0, q1);
        if (dist < minDistance) {
            minDistance = dist;
            nearestPt = LineSegment::closestPoint(q0, q1, queryPt);
            index = i;
            //-- segments are half-open, so 2nd endpoint belongs to next segment
            //-- except for last segment on non-closed sequence
            if (dist == 0.0 && queryPt.equals2D(q1)) {
                if (index < pts->size() - 1) {
                    index++;
                }
                //-- normalize index for a ring
                index = normalize(*pts, index);
            }
            if (minDistance <= 0.0) {
                break;
            }
        }
    }
    return CoordinateSequenceLocation(pts, index, nearestPt);
}

