#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>
#include <utility.h>

#include <geos/geom/CircularString.h>
#include <geos/geom/LineString.h>
#include <geos/operation/split/SplitGeometryAtVertex.h>
#include <geos/io/WKTReader.h>

using geos::geom::CoordinateXY;
using geos::geom::CircularString;
using geos::geom::LineString;
using geos::operation::split::SplitGeometryAtVertex;

namespace tut {

struct test_splitgeometryatvertex_data {
    const geos::io::WKTReader reader_;
};

typedef test_group<test_splitgeometryatvertex_data, 255> group;
typedef group::object object;

group test_splitgeometryatvertextest_group("geos::operation::split::SplitGeometryAtVertex");

template<>
template<>
void object::test<1>()
{
    set_test_name("Split LineString ZM at vertex");

    auto input = reader_.read<LineString>("LINESTRING ZM (0 3 2 3, 5 8 3 4, 2 2 4 5, 6 1 5 6)");

    {
        auto splitAtStart = SplitGeometryAtVertex::splitSimpleCurveAtVertex(*input, 0);

        ensure(splitAtStart.first->isEmpty());
        ensure(splitAtStart.first->hasZ());
        ensure(splitAtStart.first->hasM());
        ensure(splitAtStart.second->equalsIdentical(input.get()));
    }

    {
        auto splitAtEnd = SplitGeometryAtVertex::splitSimpleCurveAtVertex(*input, input->getNumPoints() - 1);

        ensure(splitAtEnd.first->equalsIdentical(input.get()));
        ensure(splitAtEnd.second->isEmpty());
        ensure(splitAtEnd.second->hasZ());
        ensure(splitAtEnd.second->hasM());
    }

    {
        auto splitInMiddle = SplitGeometryAtVertex::splitSimpleCurveAtVertex(*input, 2);

        auto expectedFirst = reader_.read("LINESTRING ZM (0 3 2 3, 5 8 3 4, 2 2 4 5)");
        auto expectedSecond = reader_.read("LINESTRING ZM (2 2 4 5, 6 1 5 6)");

        ensure(splitInMiddle.first->equalsIdentical(expectedFirst.get()));
        ensure(splitInMiddle.second->equalsIdentical(expectedSecond.get()));
    }
}

template<>
template<>
void object::test<2>()
{
    set_test_name("Split CircularString ZM at vertex");

    auto input = reader_.read<CircularString>("CIRCULARSTRING ZM (-5 0 1 2, 0 5 2 3, 5 0 3 4, 10 -5 4 5, 15 0 5 6)");

    {
        auto splitAtStart = SplitGeometryAtVertex::splitSimpleCurveAtVertex(*input, 0);

        ensure(splitAtStart.first->isEmpty());
        ensure(splitAtStart.first->hasZ());
        ensure(splitAtStart.first->hasM());
        ensure(splitAtStart.second->equalsIdentical(input.get()));
    }

    {
        auto splitAtEnd = SplitGeometryAtVertex::splitSimpleCurveAtVertex(*input, input->getNumPoints() - 1);

        ensure(splitAtEnd.first->equalsIdentical(input.get()));
        ensure(splitAtEnd.second->isEmpty());
        ensure(splitAtEnd.second->hasZ());
        ensure(splitAtEnd.second->hasM());
    }

    {
        auto splitInMiddle = SplitGeometryAtVertex::splitSimpleCurveAtVertex(*input, 2);

        auto expectedFirst = reader_.read("CIRCULARSTRING ZM (-5 0 1 2, 0 5 2 3, 5 0 3 4)");
        auto expectedSecond = reader_.read("CIRCULARSTRING ZM (5 0 3 4, 10 -5 4 5, 15 0 5 6)");

        ensure(splitInMiddle.first->equalsIdentical(expectedFirst.get()));
        ensure(splitInMiddle.second->equalsIdentical(expectedSecond.get()));
    }

    ensure_THROW(SplitGeometryAtVertex::splitSimpleCurveAtVertex(*input, 1), geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("Split LineString ZM at new point");

    auto input = reader_.read<LineString>("LINESTRING ZM (0 3 2 3, 5 8 3 4, 2 2 4 5, 6 1 5 6)");

    CoordinateXY pt{2, 3};

    ensure_THROW(SplitGeometryAtVertex::splitLineStringAtPoint(*input, 3, pt), geos::util::IllegalArgumentException);

    // Split first segment
    {
        auto [first, second] = SplitGeometryAtVertex::splitLineStringAtPoint(*input, 0, pt);

        auto expectedFirst = reader_.read("LINESTRING ZM (0 3 2 3, 2 3 2.282842712474619 3.282842712474619)");
        auto expectedSecond = reader_.read("LINESTRING ZM (2 3 2.282842712474619 3.282842712474619, 5 8 3 4, 2 2 4 5, 6 1 5 6)");

        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }

    // Split second segment
    {
        auto [first, second] = SplitGeometryAtVertex::splitLineStringAtPoint(*input, 1, pt);

        auto expectedFirst = reader_.read("LINESTRING ZM (0 3 2 3, 5 8 3 4, 2 3 3.8692269873603533 4.869226987360353)");
        auto expectedSecond = reader_.read("LINESTRING ZM (2 3 3.8692269873603533 4.869226987360353, 2 2 4 5, 6 1 5 6)");

        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }


}

}
