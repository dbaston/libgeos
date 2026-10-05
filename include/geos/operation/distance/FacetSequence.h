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

#pragma once

#include <geos/geom/Envelope.h>
#include <geos/geom/Coordinate.h>
#include <geos/operation/distance/CoordinateSequenceLocation.h>

#include <vector>

namespace geos {
namespace geom {
class CoordinateSequence;
}
}

namespace geos {
namespace operation {
namespace distance {

class FacetSequence {

private:

    const geom::CoordinateSequence* pts;
    const std::size_t start;
    const std::size_t end;
    bool curved;

    /*
    * Unlike JTS, we store the envelope in the FacetSequence so
    * that it has a clear owner.  This is helpful when making a
    * tree of FacetSequence objects (FacetSequenceTreeBuilder)
    */
    geom::Envelope env;

    double computeDistanceLineLine(const FacetSequence& facetSeq,
                                   std::vector<geom::CoordinateXY> *locs) const;

    double computeDistanceArcArc(const FacetSequence& facetSeq,
                                 std::vector<geom::CoordinateXY> *locs) const;

    double computeDistanceArcLine(const FacetSequence& facetSeq,
                                  std::vector<geom::CoordinateXY> *locs) const;

    static double computeDistancePointLine(const geom::CoordinateXY& pt,
                                           const FacetSequence& facetSeq,
                                           std::vector<geom::CoordinateXY> *locs);

    static double computeDistancePointArc(const geom::CoordinateXY& pt,
                                          const FacetSequence& facetSeq,
                                          std::vector<geom::CoordinateXY> *locs);

    static void updateNearestLocations(const geom::CoordinateXY& p0, const geom::CoordinateXY& p1, std::vector<geom::CoordinateXY> *locs);

    static void updateNearestLocationsPointLine(const geom::CoordinateXY& pt,
                                                const geom::CoordinateXY& q0, const geom::CoordinateXY& q1,
                                                std::vector<geom::CoordinateXY> *locs);

    static void updateNearestLocationsLineLine(const geom::CoordinateXY& p0, const geom::CoordinateXY& p1,
                                               const geom::CoordinateXY& q0, const geom::CoordinateXY& q1,
                                               std::vector<geom::CoordinateXY> *locs);

    CoordinateSequenceLocation nearestLocationOnLine(const geom::CoordinateXY& pt) const;

    static std::size_t normalize(const geom::CoordinateSequence& pts, std::size_t index);

    void computeEnvelope();

    double nearestLocations(const FacetSequence& facetSeq, std::vector<geom::CoordinateXY>* locs) const;

public:

    const geom::Envelope* getEnvelope() const;

    const geom::CoordinateXY* getCoordinate(std::size_t index) const;

    std::size_t size() const;

    bool isCurved() const;

    bool isPoint() const;

    double distance(const FacetSequence& facetSeq) const;

    FacetSequence(const geom::CoordinateSequence* pts, bool isCurved, std::size_t start, std::size_t end);

    /**
     * Computes the locations of the nearest points between this sequence
     * and another sequence.
     * The locations are presented in the same order as the input sequences.
     *
     * @return a pair of Coordinates for the nearest points
     */
    std::vector<geom::CoordinateXY> nearestLocations(const FacetSequence& facetSeq) const;

    /**
     * Computes the location of the nearest point of this sequence
     * to a query point.
     *
     * @param p the query point
     * @return the nearest location on this sequence
     */
    CoordinateSequenceLocation nearestLocation(const geom::CoordinateXY& p) const;

};

}
}
}
