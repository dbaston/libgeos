#include <tut/tut.hpp>
// geos
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/index/kdtree/KdTree.h>

using namespace geos::index::kdtree;

namespace tut {
// dummy data, not used
struct test_kdtree_data {
    KdTree<> build(const geos::geom::Geometry & g, double tolerance) {
        KdTree<> index(tolerance);
        auto coords = g.getCoordinates();
        for (size_t i = 0; i < coords->size(); i++) {
            index.insert(coords[i]);
        }

        return index;
    }
};

using group = test_group<test_kdtree_data>;
using object = group::object;

group test_kdtree_group("geos::index::kdtree::KdTree");

//
// Test Cases
//

template<>
template<>
void object::test<1>
()
{
    KdTree<> index(0.001);

    index.insert({1, 1});
    index.insert({1, 1});

    geos::geom::Envelope queryEnv(0, 10, 0, 10);

    auto result = index.query(queryEnv);

    ensure_equals(result.size(), 1ul);

    ensure(result[0]->isRepeated());
    ensure_equals(result[0]->getCount(), 2ul);
}


} // namespace tut

