#include <tut/tut.hpp>
// geos
#include <geos/index/kdtree/BulkKdTree.h>
#include <geos/geom/Envelope.h>

#include <random>

using namespace geos::index::kdtree;
using namespace geos::geom;

namespace tut {


// dummy data, not used
struct test_bulkkdtree_data {

};

using group = test_group<test_bulkkdtree_data>;
using object = group::object;

group test_bulkkdtree_group("geos::index::kdtree::BulkKdTree");

template<>
template<>
void object::test<1>()
{
    BulkKdTree index;
    index.insert(CoordinateXY{1, 2}, nullptr);

    index.build();
}

template<>
template<>
void object::test<2>()
{
    std::default_random_engine e;
    std::uniform_real_distribution<double> u;

    BulkKdTree index;

    std::vector<CoordinateXY> coords;
    for (std::size_t i = 0; i < 100; i++) {
        coords.emplace_back(u(e), u(e));
    }

    for (const auto& c : coords) {
        index.insert(c, nullptr);
    }

    index.build();

    index.validateConstruction();

    Envelope env(0.2, 0.4, 0.7, 0.9);
    std::set<CoordinateXY> hits;
    index.query(env, [&hits](const auto& n) {
        hits.emplace(n.getCoordinate());
    });

    for (const auto& c: coords) {
        if (env.contains(c)) {
            ensure(hits.find(c) != hits.end());
        } else {
            ensure(hits.find(c) == hits.end());
        }
    }
}




} // namespace tut

