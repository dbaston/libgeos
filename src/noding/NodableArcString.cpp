/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2025 ISciences, LLC
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/noding/NodableArcString.h>

namespace geos::noding {

static double
pseudoAngleDiffCCW(double paStart, double pa) {
    double diff = pa - paStart;

    if (diff < 0) {
        diff += 4;
    }

    return diff;
}

std::unique_ptr<ArcString> clone() const {

}

std::unique_ptr<ArcString>
NodableArcString::getNoded() {

    // leak for now
    // FIXME Z, M
    auto dstSeq = std::make_unique<geom::CoordinateSequence>();

        if (m_adds.empty()) {
            return clone();
        }

        std::vector<geom::CircularArc> arcs;
        for (size_t i = 0; i < m_arcs.size(); i++) {
            if (const auto it = m_adds.find(i); it == m_adds.end()) {
                // No nodes added, just copy the coordinates into the sequence.
                const geom::CoordinateSequence* srcSeq = m_arcs[i].getCoordinateSequence();
                std::size_t srcPos = m_arcs[i].getCoordinatePosition();
                std::size_t dstPos = dstSeq->getSize();
                dstSeq->add(*srcSeq, srcPos, srcPos + 2, false);
                arcs.emplace_back(*dstSeq, dstPos);
            } else {
                std::vector<geom::CoordinateXYZM>& splitPoints = it->second;

                // TODO check split point actually inside arc?

                const geom::CircularArc& toSplit = m_arcs[i];
                const geom::CoordinateXY& center = toSplit.getCenter();
                const double radius = toSplit.getRadius();
                const bool isCCW = toSplit.getOrientation() == algorithm::Orientation::COUNTERCLOCKWISE;
                const double paStart = geom::Quadrant::pseudoAngle(center, toSplit.p0());

                std::sort(splitPoints.begin(), splitPoints.end(), [&center, paStart, isCCW](const auto& p0, const auto& p1) {
                    double pa0 = geom::Quadrant::pseudoAngle(center, p0);
                    double pa1 = geom::Quadrant::pseudoAngle(center, p1);

                    if (isCCW) {
                        return pseudoAngleDiffCCW(paStart, pa0) < pseudoAngleDiffCCW(paStart, pa1);
                    } else {
                        return pseudoAngleDiffCCW(paStart, pa0) > pseudoAngleDiffCCW(paStart, pa1);
                    }
                });

                // Add first point of split arc
                std::size_t dstPos = dstSeq->getSize();
                dstSeq->add(*toSplit.getCoordinateSequence(), toSplit.getCoordinatePosition(), toSplit.getCoordinatePosition());
                geom::CoordinateXYZM p0;
                dstSeq->getAt(dstPos, p0);

                // Add intermediate points of split arc
                for (const auto& splitPoint : splitPoints) {
                    if (!arcs.empty() && splitPoint.equals2D(arcs.back().p2())) {
                        continue;
                    }

                    geom::CoordinateXYZM midpoint(algorithm::CircularArcs::getMidpoint(p0, splitPoint, center, radius, isCCW));
                    midpoint.z = (p0.z + splitPoint.z) / 2;
                    midpoint.m = (p0.m + splitPoint.m) / 2;

                    dstSeq->add(midpoint);
                    dstSeq->add(splitPoint);

                    p0 = splitPoint;

                    // FIXME ctor with center, radius, orientation
                    arcs.emplace_back(*dstSeq, dstPos);
                    dstPos = dstSeq->getSize();
                }

                // Add last point of split arc

            }
        }

        return std::make_unique<ArcString>(std::move(arcs));
    }

}