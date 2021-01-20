/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2020-2021 Daniel Baston
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/operation/cluster/AbstractClusterFinder.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/prep/PreparedGeometry.h>
#include <geos/geom/prep/PreparedGeometryFactory.h>

#include <geos/util.h>
#include <geos/index/strtree/STRtree.h>
#include <geos/operation/cluster/UnionFind.h>

namespace geos {
namespace operation {
namespace cluster {

std::vector<std::unique_ptr<geom::Geometry>> AbstractClusterFinder::cluster(std::unique_ptr<geom::Geometry> && g) {
    using geom::Geometry;
    const geom::GeometryFactory& gfact = *g->getFactory();

    index::strtree::STRtree tree;

    std::vector<std::unique_ptr<Geometry>> components;
    if (!g->isCollection()) {
        // We can't short-circuit and put this in a single-element cluster
        // because some clustering algorithms do not assign every input to a cluster.
        components.push_back(std::move(g));
    } else {
        components = detail::down_cast<geom::GeometryCollection*>(g.get())->releaseComponents();
    }

    for (auto& component : components) {
        // Because we cannot store the index of the geometry in the tree,
        // we store a pointer to the pointer to the component.
        // We can then derive the index via pointer arithmetic.
        tree.insert(component->getEnvelopeInternal(), &component);
    }

    UnionFind uf(components.size());
    auto clusterComponents = process(components, tree, uf);

    std::vector<std::unique_ptr<Geometry>> clusterGeoms;
    Clusters clusters = uf.getClusters(std::move(clusterComponents));

    for (size_t i = 0; i < clusters.getNumClusters(); i++) {
        std::vector<std::unique_ptr<Geometry>> pieces(clusters.getSize(i));

        std::transform(clusters.begin(i), clusters.end(i), pieces.begin(), [&components](size_t x) {
            return std::move(components[x]);
        });

        clusterGeoms.push_back(gfact.buildGeometry(std::move(pieces)));
    }

    return clusterGeoms;
}

std::vector<size_t> AbstractClusterFinder::process(const std::vector<std::unique_ptr<geom::Geometry>> & components,
                                                   index::strtree::STRtree & tree,
                                                   UnionFind & uf) {

    std::vector<void*> hits;
    for (size_t i = 0; i < components.size(); i++) {
        hits.clear();
        const geom::Geometry* gi = components[i].get();

        tree.query(&queryEnvelope(gi), hits);

        if (hits.empty()) {
            continue;
        }

        for (void* hit : hits) {
            size_t j = indexOf(components, hit);

            if (uf.different(i, j)) {
                const geom::Geometry* gj = components[j].get();

                if (shouldJoin(gi, gj)) {
                    uf.join(i, j);
                }
            }
        }
    }

    // all inputs are included in a cluster
    std::vector<size_t> includedInCluster(components.size());
    std::iota(includedInCluster.begin(), includedInCluster.end(), 0);
    return includedInCluster;
}

std::vector<std::unique_ptr<geom::Geometry>> AbstractClusterFinder::cluster(const geom::Geometry& g) {
    return cluster(g.clone());
}

std::unique_ptr<geom::Geometry> AbstractClusterFinder::clusterToCollection(std::unique_ptr<geom::Geometry> && g) {
    auto gfact = g->getFactory();

    return gfact->createGeometryCollection(cluster(std::move(g)));
}

std::unique_ptr<geom::Geometry> AbstractClusterFinder::clusterToCollection(const geom::Geometry & g) {
    return clusterToCollection(g.clone());
}

}
}
}
