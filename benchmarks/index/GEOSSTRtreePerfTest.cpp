#include <benchmark/benchmark.h>

#include <geos_c.h>

#include <memory>
#include <random>
#include <vector>
#include <geos/index/strtree/TemplateSTRtree.h>

//using GeomPtr = std::unique_ptr<GEOSGeometry, decltype(GEOSGeom_destroy)>;
using GeomPtr = GEOSGeometry*;

std::vector<GeomPtr>
random_points(std::default_random_engine & e, std::size_t n) {
    std::uniform_real_distribution<> x(0, 1);
    std::uniform_real_distribution<> y(0, 1);

    std::vector<GeomPtr> geoms;
    for (std::size_t i = 0; i < n; i++) {
        geoms.emplace_back(GEOSGeom_createPointFromXY(x(e), y(e)));
    }

    return geoms;
}

template<typename T>
GEOSSTRtree* buildTree(const T& coll, bool preSorted)
{
    GEOSSTRtree* tree = GEOSSTRtree_create(10);
    for (const auto& g : coll) {
        GEOSSTRtree_insert(tree, g, g);
    }
    if (preSorted) {
        (reinterpret_cast<geos::index::strtree::TemplateSTRtree<void*>*>(tree))->needsSorting = false;
    }

    GEOSSTRtree_query(tree, &*coll[0], [](void* item, void* userdata) {}, nullptr);

    return tree;
}

static void BM_BuildTree(benchmark::State& state) {
    std::default_random_engine eng(12345);
    initGEOS(nullptr, nullptr);
    auto geoms = random_points(eng, 10000);

    for (auto _ : state) {
        GEOSSTRtree* tree = buildTree(geoms, false);
        GEOSSTRtree_destroy(tree);
    }
}

static void BM_BuildTreePresorted(benchmark::State& state) {
    std::default_random_engine eng(12345);
    initGEOS(nullptr, nullptr);
    auto geoms = random_points(eng, 10000);

    GEOSSTRtree* tree = buildTree(geoms, false);

    // Retrieve items in tree order
    std::vector<GEOSGeometry*> sortedGeoms;
    GEOSSTRtree_iterate(tree, [](void* item, void* userdata) {
        auto* vec = reinterpret_cast<decltype(&sortedGeoms)>(userdata);
        vec->push_back(reinterpret_cast<GEOSGeometry*>(item));
    }, &sortedGeoms);
    GEOSSTRtree_destroy(tree);

    for (auto _ : state) {
        GEOSSTRtree* tree = buildTree(sortedGeoms, true);
        GEOSSTRtree_destroy(tree);
    }
}

BENCHMARK(BM_BuildTree);
BENCHMARK(BM_BuildTreePresorted);

BENCHMARK_MAIN();
