#pragma once

#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/GeometryCollection.h>

namespace geos {
namespace geom {
namespace util {

/**
 * The GeometryDimensionFilter removes components of a geometry that do not
 * have the specified topological dimension. It can be used as a postprocessing
 * operation, for example to remove lower-dimensionality components from a
 * mixed-dimension overlay result.
 */
class GeometryDimensionFilter {
public:
    /**
     * Remove components of the input that do not have the specified dimensionality.
     * A dimension value of -1 can be used to indicate that only the highest-dimension
     * components should be kept.
     *
     * @param g geometry to filter
     * @param dimension Desired output dimension. A value of -1 indicates that the
     *        highest-dimension components should be kept.
     * @return filtered geometry
     */
    static std::unique_ptr<Geometry>
    filterDimension(std::unique_ptr<Geometry> g, int dimension)
    {
        const GeometryFactory& gfact = *g->getFactory();
        std::unique_ptr<Geometry> ret;

        if (dimension == -1) {
            dimension = g->getDimension();
        }

        if (g->getGeometryTypeId() == geos::geom::GeometryTypeId::GEOS_GEOMETRYCOLLECTION) {
            auto components = static_cast<GeometryCollection*>(g.get())->releaseGeometries();
            std::remove_if(components.begin(), components.end(), [dimension](const Geometry* component) {
                return component->getDimension() == dimension;
            });
            ret = gfact.buildGeometry(std::move(components));

            return filterDimension(std::move(ret), dimension);
        } else if (g->getDimension() == dimension) {
            return g;
        } else {
            return gfact.createEmpty(dimension);
        }
    }

    static std::unique_ptr<Geometry>
    filterDimension(const Geometry& g, int dimension)
    {
        return filterDimension(g.clone(), dimension);
    }

};

}
}
}
