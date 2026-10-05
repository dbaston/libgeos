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
 * Last port: operation/distance/FacetSequenceTreeBuilder.java (f6187ee2 JTS-1.14)
 *
 **********************************************************************/

#include <geos/operation/distance/FacetSequenceTreeBuilder.h>
#include <geos/geom/CompoundCurve.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Point.h>

using namespace geos::geom;
using namespace geos::index::strtree;

namespace geos {
namespace operation {
namespace distance {

std::unique_ptr<TemplateSTRtree<const FacetSequence*>>
FacetSequenceTreeBuilder::build(const Geometry* g)
{
    auto tree = std::make_unique<FacetSequenceTree>(computeFacetSequences(g));

    tree->build();
    return tree;
}

std::vector<FacetSequence>
FacetSequenceTreeBuilder::computeFacetSequences(const Geometry* g)
{
    std::vector<FacetSequence> sections;

    class FacetSequenceAdder : public geom::GeometryComponentFilter {
        std::vector<FacetSequence>&  m_sections;

        public :
            explicit FacetSequenceAdder(std::vector<FacetSequence> & p_sections) :
                m_sections(p_sections) {}
            void
            filter_ro(const Geometry* geom) override
            {
                if(const SimpleCurve* sc = dynamic_cast<const SimpleCurve*>(geom)) {
                    const CoordinateSequence* seq = sc->getCoordinatesRO();
                    addFacetSequences(seq, sc->isCurved(), m_sections);
                }
                else if(const Point* pt = dynamic_cast<const Point*>(geom)) {
                    const CoordinateSequence* seq = pt->getCoordinatesRO();
                    addFacetSequences(seq, false, m_sections);
                }
                else if (const CompoundCurve* cc = dynamic_cast<const CompoundCurve*>(geom)) {
                    for (std::size_t i = 0; i < cc->getNumCurves(); i++) {
                        const SimpleCurve* part = cc->getCurveN(i);
                        const CoordinateSequence* seq = part->getCoordinatesRO();
                        addFacetSequences(seq, part->isCurved(), m_sections);
                    }
                }
            }
    };

    FacetSequenceAdder facetSequenceAdder(sections);
    g->apply_ro(&facetSequenceAdder);

    return sections;
}

void
FacetSequenceTreeBuilder::addFacetSequences(
    const CoordinateSequence* pts,
    bool isCurved,
    std::vector<FacetSequence>& sections)
{
    std::size_t i = 0;
    std::size_t size = pts->size();
    if ( size == 0 ) return;

    // FIXME dont break arcs
    while(i <= size - 1) {
        std::size_t end = i + FACET_SEQUENCE_SIZE + 1;
        if (end % 2 == 0) {
            end++;
        }
        // if only one point remains after this section, include it in this
        // section
        if(end >= size - 1) {
            end = size;
        }
        sections.emplace_back(pts, isCurved, i, end);
        i += FACET_SEQUENCE_SIZE;
    }
}

}
}
}

