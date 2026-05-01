#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>
#include <utility.h>

#include <geos/geom/CircularString.h>
#include <geos/geom/LineString.h>
#include <geos/operation/split/SplitLinealAtPoint.h>
#include <geos/io/WKTReader.h>

using geos::geom::CoordinateXY;
using geos::geom::CircularString;
using geos::geom::LineString;
using geos::operation::split::SplitLinealAtPoint;

namespace tut {

struct test_splitlinealatpoint_data {
    const geos::io::WKTReader reader_;
};

typedef test_group<test_splitlinealatpoint_data, 255> group;
typedef group::object object;

group test_splitlinealatpointtest_group("geos::operation::split::SplitLinealAtPoint");

template<>
template<>
void object::test<1>()
{
    set_test_name("Split LineString ZM at vertex");

    auto input = reader_.read<LineString>("LINESTRING ZM (0 3 2 3, 5 8 3 4, 2 2 4 5, 6 1 5 6)");

    {
        auto splitAtStart = SplitLinealAtPoint::splitSimpleCurveAtVertex(*input, 0);

        ensure(splitAtStart.first->isEmpty());
        ensure(splitAtStart.first->hasZ());
        ensure(splitAtStart.first->hasM());
        ensure(splitAtStart.second->equalsIdentical(input.get()));
    }

    {
        auto splitAtEnd = SplitLinealAtPoint::splitSimpleCurveAtVertex(*input, input->getNumPoints() - 1);

        ensure(splitAtEnd.first->equalsIdentical(input.get()));
        ensure(splitAtEnd.second->isEmpty());
        ensure(splitAtEnd.second->hasZ());
        ensure(splitAtEnd.second->hasM());
    }

    {
        auto splitInMiddle = SplitLinealAtPoint::splitSimpleCurveAtVertex(*input, 2);

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
        auto splitAtStart = SplitLinealAtPoint::splitSimpleCurveAtVertex(*input, 0);

        ensure(splitAtStart.first->isEmpty());
        ensure(splitAtStart.first->hasZ());
        ensure(splitAtStart.first->hasM());
        ensure(splitAtStart.second->equalsIdentical(input.get()));
    }

    {
        auto splitAtEnd = SplitLinealAtPoint::splitSimpleCurveAtVertex(*input, input->getNumPoints() - 1);

        ensure(splitAtEnd.first->equalsIdentical(input.get()));
        ensure(splitAtEnd.second->isEmpty());
        ensure(splitAtEnd.second->hasZ());
        ensure(splitAtEnd.second->hasM());
    }

    {
        auto splitInMiddle = SplitLinealAtPoint::splitSimpleCurveAtVertex(*input, 2);

        auto expectedFirst = reader_.read("CIRCULARSTRING ZM (-5 0 1 2, 0 5 2 3, 5 0 3 4)");
        auto expectedSecond = reader_.read("CIRCULARSTRING ZM (5 0 3 4, 10 -5 4 5, 15 0 5 6)");

        ensure(splitInMiddle.first->equalsIdentical(expectedFirst.get()));
        ensure(splitInMiddle.second->equalsIdentical(expectedSecond.get()));
    }

    ensure_THROW(SplitLinealAtPoint::splitSimpleCurveAtVertex(*input, 1), geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("Split LineString ZM at new point");

    auto input = reader_.read<LineString>("LINESTRING ZM (0 3 2 3, 5 8 3 4, 2 2 4 5, 6 1 5 6)");

    CoordinateXY pt{2, 3};

    ensure_THROW(SplitLinealAtPoint::splitLineStringAtPoint(*input, 3, pt), geos::util::IllegalArgumentException);

    // Split first segment
    {
        auto [first, second] = SplitLinealAtPoint::splitLineStringAtPoint(*input, 0, pt);

        auto expectedFirst = reader_.read("LINESTRING ZM (0 3 2 3, 2 3 2.282842712474619 3.282842712474619)");
        auto expectedSecond = reader_.read("LINESTRING ZM (2 3 2.282842712474619 3.282842712474619, 5 8 3 4, 2 2 4 5, 6 1 5 6)");

        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }

    // Split second segment
    {
        auto [first, second] = SplitLinealAtPoint::splitLineStringAtPoint(*input, 1, pt);

        auto expectedFirst = reader_.read("LINESTRING ZM (0 3 2 3, 5 8 3 4, 2 3 3.8692269873603533 4.869226987360353)");
        auto expectedSecond = reader_.read("LINESTRING ZM (2 3 3.8692269873603533 4.869226987360353, 2 2 4 5, 6 1 5 6)");

        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }

}

template<>
template<>
void object::test<4>()
{
    set_test_name("Split CompoundCurve at existing vertices");

    auto input = reader_.read<geos::geom::CompoundCurve>("COMPOUNDCURVE (CIRCULARSTRING(2 8, 4 7, 1 9, 3 15, 8 16), (8 16, 8 10, 6 14, 4 12), CIRCULARSTRING (4 12, 4 10, 2 8))");

    // Split at first point
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 0, 0, CoordinateXY{2, 8});

        ensure(first->isEmpty());
        ensure(second->equalsIdentical(input.get()));
    }

    // Split at intermediate point of first curve
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 0, 2, CoordinateXY{1, 9});

        auto expectedFirst = reader_.read("CIRCULARSTRING (2 8, 4 7, 1 9)");
        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);

        auto expectedSecond = reader_.read("COMPOUNDCURVE (CIRCULARSTRING(1 9, 3 15, 8 16), (8 16, 8 10, 6 14, 4 12), CIRCULARSTRING (4 12, 4 10, 2 8))");
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }

    // Split at last point of first curve
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 0, 2, CoordinateXY{8, 16});

        ensure(first->equalsIdentical(input->getCurveN(0)));

        auto expectedSecond = reader_.read("COMPOUNDCURVE ((8 16, 8 10, 6 14, 4 12), CIRCULARSTRING (4 12, 4 10, 2 8))");
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }

    // Split at first point of second curve
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 1, 0, CoordinateXY{8, 16});

        ensure(first->equalsIdentical(input->getCurveN(0)));

        auto expectedSecond = reader_.read("COMPOUNDCURVE ((8 16, 8 10, 6 14, 4 12), CIRCULARSTRING (4 12, 4 10, 2 8))");
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }

    // Split at intermate point of second curve
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 1, 1, CoordinateXY{6, 14});

        auto expectedFirst = reader_.read("COMPOUNDCURVE (CIRCULARSTRING(2 8, 4 7, 1 9, 3 15, 8 16), (8 16, 8 10, 6 14))");
        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);

        auto expectedSecond = reader_.read("COMPOUNDCURVE ((6 14, 4 12), CIRCULARSTRING (4 12, 4 10, 2 8))");
        ensure_equals_exact_geometry_xyzm(second.get(), expectedSecond.get(), 0.0);
    }

    // Split at last point of second curve
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 1, 3, CoordinateXY{4, 12});

        auto expectedFirst = reader_.read("COMPOUNDCURVE (CIRCULARSTRING(2 8, 4 7, 1 9, 3 15, 8 16), (8 16, 8 10, 6 14, 4 12))");
        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);

        ensure(second->equalsIdentical(input->getCurveN(2)));
    }

    // Split at first point of third curve
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 2, 0, CoordinateXY{4, 12});

        auto expectedFirst = reader_.read("COMPOUNDCURVE (CIRCULARSTRING(2 8, 4 7, 1 9, 3 15, 8 16), (8 16, 8 10, 6 14, 4 12))");
        ensure_equals_exact_geometry_xyzm(first.get(), expectedFirst.get(), 0.0);

        ensure(second->equalsIdentical(input->getCurveN(2)));
    }

    // Split at last point of third curve
    {
        auto [first, second] = SplitLinealAtPoint::splitCompoundCurveAtPoint(*input, 2, 2, CoordinateXY{2, 8});

        ensure(first->equalsIdentical(input.get()));
        ensure(second->isEmpty());
    }

}

}
