/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2012  Sandro Santilli <strk@kbt.io>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * NOTE: this is not in JTS. JTS has a snapround/GeometryNoder though
 *
 **********************************************************************/

#include <geos/algorithm/CircularArcIntersector.h>
#include <geos/noding/GeometryNoder.h>
#include <geos/noding/SegmentString.h>
#include <geos/noding/NodedSegmentString.h>
#include <geos/noding/OrientedCoordinateArray.h>
#include <geos/noding/Noder.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CompoundCurve.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/CircularString.h>
#include <geos/geom/MultiCurve.h>
#include <geos/geom/LineString.h>

#include <geos/noding/ArcIntersectionAdder.h>
#include <geos/noding/IteratedNoder.h>
#include <geos/noding/NodableArcString.h>
#include <geos/noding/SimpleNoder.h>

#include <geos/algorithm/LineIntersector.h>
#include <geos/noding/IntersectionAdder.h>
#include <geos/noding/MCIndexNoder.h>

#include <geos/noding/snapround/MCIndexSnapRounder.h>

#include <memory> // for unique_ptr
#include <iostream>


namespace geos {
namespace noding { // geos.noding

namespace {

/**
 * Add every linear element in a geometry into PathString vector
 */
class PathStringExtractor: public geom::GeometryComponentFilter {
public:
    PathStringExtractor(std::vector<std::unique_ptr<PathString>> & to,
                           bool constructZ,
                           bool constructM,
                           const void* context)
        : _to(to)
        , _constructZ(constructZ)
        , _constructM(constructM)
        , _context(context)
    {}

    void
    filter_ro(const geom::Geometry* g) override
    {
        if (g->isEmpty()) {
            return;
        }

        if(const auto* ls = dynamic_cast<const geom::LineString*>(g)) {
            auto coord = ls->getSharedCoordinates();
            auto ss = std::make_unique<NodedSegmentString>(coord, _constructZ, _constructM, _context);
            _to.push_back(std::move(ss));
        } else if (const auto* cs = dynamic_cast<const geom::CircularString*>(g)) {
            const auto& coords = cs->getSharedCoordinates();
            auto arcs = cs->getArcs();

            auto as = std::make_unique<NodableArcString>(std::move(arcs), coords, _constructZ, _constructM, _context);
            _to.push_back(std::move(as));
        } else if (const auto* cc = dynamic_cast<const geom::CompoundCurve*>(g)) {
            for (std::size_t i = 0; i < cc->getNumCurves(); i++) {
                filter_ro(cc->getCurveN(i));
            }
        }
    }
private:
    std::vector<std::unique_ptr<PathString>>& _to;
    bool _constructZ;
    bool _constructM;
    const void* _context;

    PathStringExtractor(PathStringExtractor const&); /*= delete*/
    PathStringExtractor& operator=(PathStringExtractor const&); /*= delete*/
};

}


/* public static */
std::unique_ptr<geom::Geometry>
GeometryNoder::node(const geom::Geometry& geom)
{
    GeometryNoder noder(geom);
    return noder.getNoded();
}

std::unique_ptr<geom::Geometry>
GeometryNoder::node(const geom::Geometry& geom1, const geom::Geometry& geom2)
{
    GeometryNoder noder(geom1, geom2);
    return noder.getNoded();
}

/* public */
GeometryNoder::GeometryNoder(const geom::Geometry& g)
    :
    argGeom1(&g),
    argGeom2(nullptr),
    argGeomHasCurves(g.hasCurvedComponents()),
    onlyFirstGeomEdges(false)
{}

GeometryNoder::GeometryNoder(const geom::Geometry& g1, const geom::Geometry& g2)
    :
    argGeom1(&g1),
    argGeom2(&g2),
    argGeomHasCurves(g1.hasCurvedComponents() || g2.hasCurvedComponents()),
    onlyFirstGeomEdges(false)
{}

GeometryNoder::~GeometryNoder() = default;

/* private */
std::unique_ptr<geom::Geometry>
GeometryNoder::toGeometry(std::vector<std::unique_ptr<PathString>>& nodedEdges) const
{
    const geom::GeometryFactory* geomFact = argGeom1->getFactory();

    std::set< OrientedCoordinateArray > ocas;

    // Create a geometry out of the noded substrings.
    std::vector<std::unique_ptr<geom::Geometry>> lines;
    lines.reserve(nodedEdges.size());

    bool resultArcs = false;
    for(auto& path :  nodedEdges) {
        if (onlyFirstGeomEdges && path->getData() != argGeom1) {
            continue;
        }

        const auto& coords = path->getCoordinates();

        bool isLinear = dynamic_cast<SegmentString*>(path.get());

        // TODO: Make OrientedCoordinateArray not require strict equality of arc control points

        OrientedCoordinateArray oca1(*coords);
        // Check if an equivalent edge is known
        if(ocas.insert(oca1).second) {
            if (isLinear) {
                lines.push_back(geomFact->createLineString(coords));
            } else {
                resultArcs = true;
                lines.push_back(geomFact->createCircularString(coords));
            }
        }
    }

    if (resultArcs) {
        return geomFact->createMultiCurve(std::move(lines));
    } else {
        return geomFact->createMultiLineString(std::move(lines));
    }
}

/* public */
std::unique_ptr<geom::Geometry>
GeometryNoder::getNoded()
{
    if (argGeom1->isEmpty() && (argGeom2 == nullptr || argGeom2->isEmpty()))
        return argGeom1->clone();

    std::vector<std::unique_ptr<PathString>> lineList;

    extractPathStrings(*argGeom1, lineList);
    if (argGeom2 != nullptr) {
        extractPathStrings(*argGeom2, lineList);
    }

    Noder& p_noder = getNoder();
    p_noder.computePathNodes(PathString::toRawPointerVector(lineList));
    auto nodedEdges = p_noder.getNodedPaths();

    std::unique_ptr<geom::Geometry> noded = toGeometry(nodedEdges);

    return noded;
}

/* private static */
void
GeometryNoder::extractPathStrings(const geom::Geometry& g,
                                  std::vector<std::unique_ptr<PathString>>& to)
{
    PathStringExtractor ex(to, g.hasZ(), g.hasM(), &g);
    g.apply_ro(&ex);
}

/* private */
Noder&
GeometryNoder::getNoder()
{
    if(!noder) {
        const geom::PrecisionModel* pm = argGeom1->getFactory()->getPrecisionModel();
        if (argGeomHasCurves) {
            noder = std::make_unique<SimpleNoder>();

            m_cai = std::make_unique<algorithm::CircularArcIntersector>(argGeom1->getPrecisionModel());
            m_aia = std::make_unique<ArcIntersectionAdder>(*m_cai);
            detail::down_cast<SimpleNoder*>(noder.get())->setArcIntersector(*m_aia);
        } else {
            noder = std::make_unique<IteratedNoder>(pm);
        }
    }
    return *noder;
}

void
GeometryNoder::setOnlyFirstGeomEdges(bool p_onlyFirstGeomEdges)
{
    onlyFirstGeomEdges = p_onlyFirstGeomEdges;
}

} // namespace geos.noding
} // namespace geos
