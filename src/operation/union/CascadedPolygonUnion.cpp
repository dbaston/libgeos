/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2011 Sandro Santilli <strk@kbt.io>
 * Copyright (C) 2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: operation/union/CascadedPolygonUnion.java r487 (JTS-1.12+)
 * Includes custom code to deal with https://trac.osgeo.org/geos/ticket/837
 *
 **********************************************************************/


#include <geos/operation/union/CascadedPolygonUnion.h>
#include <geos/operation/union/OverlapUnion.h>
#include <geos/operation/overlay/OverlayOp.h>
#include <geos/geom/HeuristicOverlay.h>
#include <geos/geom/Dimension.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/util/PolygonExtracter.h>
#include <geos/index/strtree/TemplateSTRtree.h>

// std
#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>
#include <sstream>

#include <geos/operation/valid/IsValidOp.h>
#include <geos/operation/IsSimpleOp.h>
#include <geos/algorithm/BoundaryNodeRule.h>
#include <geos/util/TopologyException.h>
#include <string>
#include <iomanip>

//#define GEOS_DEBUG_CASCADED_UNION 1
//#define GEOS_DEBUG_CASCADED_UNION_PRINT_INVALID 1

namespace {

#if GEOS_DEBUG
inline bool
check_valid(const geos::geom::Geometry& g, const std::string& label, bool doThrow = false, bool validOnly = false)
{
    using namespace geos;

    if(g.isLineal()) {
        if(! validOnly) {
            operation::IsSimpleOp sop(g, algorithm::BoundaryNodeRule::getBoundaryEndPoint());
            if(! sop.isSimple()) {
                if(doThrow) {
                    throw geos::util::TopologyException(
                        label + " is not simple");
                }
                return false;
            }
        }
    }
    else {
        operation::valid::IsValidOp ivo(&g);
        if(! ivo.isValid()) {
            using operation::valid::TopologyValidationError;
            TopologyValidationError* err = ivo.getValidationError();
#ifdef GEOS_DEBUG_CASCADED_UNION
            std::cerr << label << " is INVALID: "
                      << err->toString()
                      << " (" << std::setprecision(20)
                      << err->getCoordinate() << ")"
                      << std::endl
#ifdef GEOS_DEBUG_CASCADED_UNION_PRINT_INVALID
                      << "<a>" << std::endl
                      << g.toString() << std::endl
                      << "</a>" << std::endl
#endif
                      ;
#endif // GEOS_DEBUG_CASCADED_UNION
            if(doThrow) {
                throw geos::util::TopologyException(
                    label + " is invalid: " + err->toString(),
                    err->getCoordinate());
            }
            return false;
        }
    }
    return true;
}
#endif

} // anonymous namespace


namespace geos {
namespace operation { // geos.operation
namespace geounion {  // geos.operation.geounion

// ////////////////////////////////////////////////////////////////////////////
void
GeometryListHolder::deleteItem(geom::Geometry* item)
{
    delete item;
}

// ////////////////////////////////////////////////////////////////////////////
std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::Union(std::vector<geom::Polygon*>* polys)
{
    CascadedPolygonUnion op(polys);
    return op.Union();
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::Union(std::vector<geom::Polygon*>* polys, UnionStrategy* unionFun)
{
    CascadedPolygonUnion op(polys, unionFun);
    return op.Union();
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::Union(const geom::MultiPolygon* multipoly)
{
    std::vector<geom::Polygon*> polys;

    for(const auto& g : *multipoly) {
        polys.push_back(dynamic_cast<geom::Polygon*>(g.get()));
    }

    CascadedPolygonUnion op(&polys);
    return op.Union();
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::Union()
{
    if(inputPolys->empty()) {
        return nullptr;
    }

    geomFactory = inputPolys->front()->getFactory();

    /*
     * A spatial index to organize the collection
     * into groups of close geometries.
     * This makes unioning more efficient, since vertices are more likely
     * to be eliminated on each round.
     */

    index::strtree::TemplateSTRtree<const geom::Geometry*> index(10, inputPolys->size());
    for (const auto& p : *inputPolys) {
        index.insert(p);
    }

    // TODO avoid creating this vector and run binaryUnion off the iterators directly
    std::vector<const geom::Geometry*> geoms(index.items().begin(), index.items().end());

    return binaryUnion(geoms, 0, geoms.size());
}


std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::binaryUnion(const std::vector<const geom::Geometry*> & geoms,
                                  std::size_t start, std::size_t end)
{
    if(end - start <= 1) {
        return unionSafe(geoms[start], nullptr);
    }
    else if(end - start == 2) {
        return unionSafe(geoms[start], geoms[start + 1]);
    }
    else {
        // recurse on both halves of the list
        std::size_t mid = (end + start) / 2;
        std::unique_ptr<geom::Geometry> g0(binaryUnion(geoms, start, mid));
        std::unique_ptr<geom::Geometry> g1(binaryUnion(geoms, mid, end));
        return unionSafe(std::move(g0), std::move(g1));
    }
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::unionSafe(const geom::Geometry* g0, const geom::Geometry* g1) const
{
    if(g0 == nullptr && g1 == nullptr) {
        return nullptr;
    }

    if(g0 == nullptr) {
        return g1->clone();
    }
    if(g1 == nullptr) {
        return g0->clone();
    }

    return unionActual(g0, g1);
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::unionSafe(std::unique_ptr<geom::Geometry> && g0, std::unique_ptr<geom::Geometry> && g1)
{
    if(g0 == nullptr && g1 == nullptr) {
        return nullptr;
    }

    if(g0 == nullptr) {
        return std::move(g1);
    }
    if(g1 == nullptr) {
        return std::move(g0);
    }

    return unionActual(std::move(g0), std::move(g1));
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::unionActual(const geom::Geometry* g0, const geom::Geometry* g1) const
{
    std::unique_ptr<geom::Geometry> ug;
    ug = unionFunction->Union(g0, g1);
    return restrictToPolygons(std::move(ug));
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::unionActual(std::unique_ptr<geom::Geometry> && g0, std::unique_ptr<geom::Geometry> && g1) const
{
    std::unique_ptr<geom::Geometry> ug;
    ug = unionFunction->Union(std::move(g0), std::move(g1));
    return restrictToPolygons(std::move(ug));
}

std::unique_ptr<geom::Geometry>
CascadedPolygonUnion::restrictToPolygons(std::unique_ptr<geom::Geometry> g)
{
    using namespace geom;

    if(g->isPolygonal()) {
        return g;
    }

    auto gfact = g->getFactory();
    auto coordDim = g->getCoordinateDimension();

    auto coll = dynamic_cast<GeometryCollection*>(g.get());
    if (coll) {
        // Release polygons from the collection and re-form into MultiPolygon
        auto components = coll->releaseGeometries();
        components.erase(std::remove_if(components.begin(), components.end(), [](const std::unique_ptr<Geometry> & cmp) {
            return !cmp->isPolygonal();
        }), components.end());

        return gfact->createMultiPolygon(std::move(components));
    } else {
        // Not polygonal and not a collection? No polygons here.
        return gfact->createPolygon(coordDim);
    }
}

/************************************************************************/

using operation::overlay::OverlayOp;

std::unique_ptr<geom::Geometry>
ClassicUnionStrategy::Union(const geom::Geometry* g0, const geom::Geometry* g1)
{
    // TODO make an rvalue overload for this that can consume its inputs.
    // At a minimum, a copy in the buffer fallback can be eliminated.
    try {
        // return SnapIfNeededOverlayOp.union(g0, g1);
        return geom::HeuristicOverlay(g0, g1, overlay::OverlayOp::opUNION);
    }
    catch (const util::TopologyException &ex) {
        // union-by-buffer only works for polygons
        if (g0->getDimension() != 2 || g1->getDimension() != 2)
          throw ex;
        return unionPolygonsByBuffer(g0, g1);
    }
}

bool
ClassicUnionStrategy::isFloatingPrecision() const
{
  return true;
}

/*private*/
std::unique_ptr<geom::Geometry>
ClassicUnionStrategy::unionPolygonsByBuffer(const geom::Geometry* g0, const geom::Geometry* g1)
{
    std::vector<std::unique_ptr<geom::Geometry>> geoms;
    geoms.push_back(g0->clone());
    geoms.push_back(g1->clone());
    std::unique_ptr<geom::GeometryCollection> coll = g0->getFactory()->createGeometryCollection(std::move(geoms));
    return coll->buffer(0);
}





} // namespace geos.operation.union
} // namespace geos.operation
} // namespace geos
