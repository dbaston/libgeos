#include <tut/tut.hpp>
#include "utility.h"

#include <geos/noding/NodableArcString.h>

using geos::algorithm::Orientation;
using geos::geom::CircularArc;
using geos::geom::Ordinate;
using geos::geom::CoordinateXY;
using geos::geom::CoordinateXYZM;
using geos::noding::NodableArcString;

namespace tut {

struct test_nodablearcstring_data {

    using XY = CoordinateXY;
    using XYZ = Coordinate;
    using XYM = geos::geom::CoordinateXYM;
    using XYZM = CoordinateXYZM;

    static void test_add_points(const CircularArc& arc, const std::vector<CoordinateXY>& coords,
                                const std::vector<CircularArc>& expected, bool reversed=false) {
        std::vector<CircularArc> arcs;
        arcs.push_back(arc);
        NodableArcString nas(arcs, nullptr, false, false, nullptr);

        for (const auto& coord : coords) {
            nas.addInt(coord, 0);
        }

        auto noded = nas.getNoded();

        ensure_equals(noded->getSize(), expected.size());

        for (std::size_t i = 0; i < expected.size(); i++) {
            ensure_arc_equals(noded->getArc(i), expected[i]);
        }

        if (!reversed) {
            const auto revArc = arc.reverse();

            std::vector<CircularArc> revExpected;
            for (const auto& x : expected) {
                revExpected.push_back(x.reverse());
            }
            std::reverse(revExpected.begin(), revExpected.end());

            test_add_points(revArc, coords, revExpected, true);
        }
    }

    static void ensure_arc_equals(const CircularArc& actual, const CircularArc& expected) {
        // FIXME check Z and M
        ensure_equals_xy(actual.p0(), expected.p0());
        ensure_equals_xy(actual.p2(), expected.p2());
        ensure_equals_xy(actual.getCenter(), expected.getCenter());
        ensure_equals(actual.getRadius(), expected.getRadius());
        ensure_equals(actual.getOrientation(), expected.getOrientation());
    }
};

typedef test_group<test_nodablearcstring_data> group;
typedef group::object object;

group test_nodablearcstring_group("geos::noding::NodableArcString");

template<>
template<>
void object::test<1>()
{
    set_test_name("CW half-circle, upper half-plane");

    CircularArc in = CircularArc::create(CoordinateXY{-5, 0}, CoordinateXY{0, 5}, CoordinateXY{5, 0});

    std::vector<CoordinateXY> coords;
    coords.emplace_back(4, 3);
    coords.emplace_back(3, 4);
    coords.emplace_back(-3, 4);
    coords.emplace_back(-4, 3);

    std::vector<CircularArc> expected;
    expected.push_back(CircularArc::create(XY{-5, 0}, {-4, 3}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{-4, 3}, {-3, 4}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{-3, 4}, {3, 4}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{3, 4}, {4, 3}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{4, 3}, {5, 0}, {0, 0}, 5, Orientation::CLOCKWISE));

    test_add_points(in, coords, expected);
}

template<>
template<>
void object::test<2>()
{
    set_test_name("CW half-circle, right half-plane");

    CircularArc in = CircularArc::create(CoordinateXY{0, 5}, CoordinateXY{5, 0}, CoordinateXY{0, -5});

    std::vector<CoordinateXY> coords;
    coords.emplace_back(4, -3);
    coords.emplace_back(4, 3);
    coords.emplace_back(3, -4);
    coords.emplace_back(3, 4);
    coords.emplace_back(5, 0);

    std::vector<CircularArc> expected;
    expected.push_back(CircularArc::create(XY{0, 5}, {3, 4}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{3, 4}, {4, 3}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{4, 3}, {5, 0}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{5, 0}, {4, -3}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{4, -3}, {3, -4}, {0, 0}, 5, Orientation::CLOCKWISE));
    expected.push_back(CircularArc::create(XY{3, -4}, {0, -5}, {0, 0}, 5, Orientation::CLOCKWISE));

    test_add_points(in, coords, expected);
}

template<>
template<>
void object::test<3>()
{
    set_test_name("no points added");
    CircularArc in = CircularArc::create(CoordinateXY{-1, 0}, CoordinateXY{0, 1}, CoordinateXY{1, 0});

    std::vector<CoordinateXY> coords;
    std::vector<CircularArc> expected;
    expected.push_back(in);
    test_add_points(in, coords, expected);
}

template<>
template<>
void object::test<4>()
{
    set_test_name("Z/M interpolated");

    CoordinateSequence seq = CoordinateSequence::XYZM(3);
    CoordinateXYZM p0{0, 5, 6, 2};
    CoordinateXYZM p1{5, 0, 7, 3};
    CoordinateXYZM p2{4, -3, 9, 1};

    seq.setAt(p0, 0);
    seq.setAt(p1, 1);
    seq.setAt(p2, 2);

    CircularArc arc (seq, 0);

    CoordinateXYZM intPt{4, 3, 13, 5};

    std::vector<CircularArc> in { arc };
    NodableArcString nas(std::move(in), nullptr, true, true, nullptr);

    nas.addIntersection( intPt, 0);

    auto noded = nas.getNoded();

    auto nodedCoords = noded->releaseCoordinates();

    ensure_equals(nodedCoords->getSize(), 5u);

    ensure_equals_xyzm(nodedCoords->getAt<CoordinateXYZM>(0), p0);

    CoordinateXYZM midpoint0 = nodedCoords->getAt<CoordinateXYZM>(1);
    ensure_equals(midpoint0.z, (p0.z + intPt.z) / 2);
    ensure_equals(midpoint0.m, (p0.m + intPt.m) / 2);

    ensure_equals_xyzm(nodedCoords->getAt<CoordinateXYZM>(2), intPt);

    CoordinateXYZM midpoint1 = nodedCoords->getAt<CoordinateXYZM>(3);
    ensure_equals(midpoint1.z, (intPt.z + p2.z) / 2);
    ensure_equals(midpoint1.m, (intPt.m + p2.m) / 2);

    ensure_equals_xyzm(nodedCoords->getAt<CoordinateXYZM>(4), p2);
}

}