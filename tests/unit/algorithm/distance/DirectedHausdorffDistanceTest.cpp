//
// Test Suite for geos::algorithm::distance::DirectedHausdorffDistance
// Ported from JTS DirectedHausdorffDistanceTest (locationtech/jts#1182)
// plus the discrete-under-estimate witness from DiscreteHausdorffDistance.

#include <tut/tut.hpp>
#include <tut/tut_macros.hpp>

#include <geos/algorithm/distance/DirectedHausdorffDistance.h>
#include <geos/algorithm/distance/DiscreteHausdorffDistance.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/io/WKTReader.h>
#include <geos/util/IllegalArgumentException.h>

#include <cmath>
#include <memory>
#include <string>

#include "utility.h"

using geos::algorithm::distance::DirectedHausdorffDistance;
using geos::algorithm::distance::DiscreteHausdorffDistance;
using geos::geom::CoordinateXY;
using geos::geom::Geometry;
using geos::geom::GeometryFactory;
using geos::geom::PrecisionModel;

namespace tut {

struct test_directedhausdorffdistance_data {
    test_directedhausdorffdistance_data()
        : pm()
        , gf(GeometryFactory::create(&pm))
        , reader(gf.get())
    {}

    static constexpr double TOLERANCE = 0.001;

    std::unique_ptr<Geometry>
    read(const std::string& wkt) const
    {
        return reader.read(wkt);
    }

    void
    checkDistance(const std::string& wkt1, const std::string& wkt2,
                  double expectedDistance) const
    {
        auto g1 = read(wkt1);
        auto g2 = read(wkt2);

        double dist = DirectedHausdorffDistance::distance(*g1, *g2);
        ensure(std::fabs(dist - expectedDistance) <= TOLERANCE);
    }

    void checkDistance(const std::string& wkt1, const std::string& wkt2,
                       double tolerance, const std::string& wktExpected) const
    {
        auto g1 = read(wkt1);
        auto g2 = read(wkt2);

        auto pts = DirectedHausdorffDistance::distancePoints(*g1, *g2, tolerance);
        ensure(pts.has_value());

        auto seq = std::make_unique<geos::geom::CoordinateSequence>(2);
        seq->setAt(pts.value()[0], 0);
        seq->setAt(pts.value()[1], 1);

        auto result = g1->getFactory()->createLineString(std::move(seq));
        ensure_equals_exact_geometry(static_cast<const Geometry*>(result.get()), read(wktExpected).get(), TOLERANCE);
    }

    void checkDistanceStartPtLen(const std::string wkt1, const std::string& wkt2,
                                 const std::string& wktExpected, double resultTolerance) const
    {
        auto g1 = read(wkt1);
        auto g2 = read(wkt2);

        auto pts = DirectedHausdorffDistance::distancePoints(*g1, *g2);
        ensure(pts.has_value());

        auto seq = std::make_unique<geos::geom::CoordinateSequence>(2);
        seq->setAt(pts.value()[0], 0);
        seq->setAt(pts.value()[1], 1);

        auto result = g1->getFactory()->createLineString(std::move(seq));
        auto expected = reader.read<LineString>(wktExpected);

        const CoordinateXY& resultPt = result->getCoordinatesRO()->getAt<CoordinateXY>(0);
        const CoordinateXY& expectedPt = expected->getCoordinatesRO()->getAt<CoordinateXY>(0);

        ensure_equals_xy(expectedPt, resultPt, resultTolerance);

        double distResult = result->getLength();
        double distExpected = expected->getLength();
        ensure_equals("distance is not within tolerance of expected", distExpected, distResult, resultTolerance);
    }

    void
    checkDistance(const std::string& wkt1, const std::string& wkt2,
                  double tolerance, double expectedDistance) const
    {
        auto g1 = read(wkt1);
        auto g2 = read(wkt2);
        double dist = DirectedHausdorffDistance::distance(*g1, *g2, tolerance);
        ensure(std::fabs(dist - expectedDistance) <= TOLERANCE);
    }

    void
    checkDistance(const std::string& wkt1, const std::string& wkt2,
                  const std::string& wktExpected) const
    {
        auto g1 = read(wkt1);
        auto g2 = read(wkt2);
        auto pts = DirectedHausdorffDistance::distancePoints(*g1, *g2);
        ensure(pts.has_value());

        auto seq = std::make_unique<geos::geom::CoordinateSequence>(2);
        seq->setAt(pts.value()[0], 0);
        seq->setAt(pts.value()[1], 1);

        auto result = g1->getFactory()->createLineString(std::move(seq));
        ensure_equals_exact_geometry(static_cast<const Geometry*>(result.get()), read(wktExpected).get(), TOLERANCE);

        ensure_equals_exact_geometry(static_cast<const Geometry*>(result.get()), read(wktExpected).get(), TOLERANCE);
    }

    void
    checkHausdorff(const std::string& wkt1, const std::string& wkt2,
                   const std::string& wktExpected) const
    {
        std::unique_ptr<Geometry> g1 = read(wkt1);
        std::unique_ptr<Geometry> g2 = read(wkt2);

        auto pts = DirectedHausdorffDistance::hausdorffDistancePoints(*g1, *g2);
        ensure(pts.has_value());

        auto seq = std::make_unique<geos::geom::CoordinateSequence>(2);
        seq->setAt(pts.value()[0], 0);
        seq->setAt(pts.value()[1], 1);

        auto result = g1->getFactory()->createLineString(std::move(seq));
        ensure_equals_exact_geometry(static_cast<const Geometry*>(result.get()), read(wktExpected).get(), TOLERANCE);
    }

    void
    checkDistanceEmpty(const std::string& a, const std::string& b) const
    {
        std::unique_ptr<Geometry> g1 = read(a);
        std::unique_ptr<Geometry> g2 = read(b);
        ensure(!DirectedHausdorffDistance::distancePoints(*g1, *g2).has_value());
        ensure(std::isnan(DirectedHausdorffDistance::distance(*g1, *g2)));
        ensure(std::isnan(DirectedHausdorffDistance::hausdorffDistance(*g1, *g2)));
    }

    void
    checkFullyWithinDistance(const std::string& a, const std::string& b,
                             double distance, bool expected) const
    {
        std::unique_ptr<Geometry> g1 = read(a);
        std::unique_ptr<Geometry> g2 = read(b);
        bool result = DirectedHausdorffDistance::isFullyWithinDistance(*g1, *g2, distance);
        ensure_equals(result, expected);
    }

    void checkFullyWithinDistanceEmpty(const std::string& a, const std::string& b) const
    {
        checkFullyWithinDistance(a, b, 0, false);
        checkFullyWithinDistance(b, a, 0, false);
        checkFullyWithinDistance(a, b, 1, false);
        checkFullyWithinDistance(b, a, 1, false);
        checkFullyWithinDistance(a, b, 1000, false);
        checkFullyWithinDistance(b, a, 1000, false);
    }

    PrecisionModel pm;
    GeometryFactory::Ptr gf;
    geos::io::WKTReader reader;
};

typedef test_group<test_directedhausdorffdistance_data> group;
typedef group::object object;

group test_DirectedHausdorffDistance_group(
    "geos::algorithm::distance::DirectedHausdorffDistance");

template<>
template<>
void object::test<1>()
{
    set_test_name("testEmpty");

    checkDistanceEmpty("POINT EMPTY", "POINT (1 1)");
    checkDistanceEmpty("LINESTRING EMPTY", "LINESTRING (0 0, 2 1)");
    checkDistanceEmpty("POLYGON EMPTY", "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))");
}

template<>
template<>
void object::test<2>()
{
    set_test_name("testZeroTolerancePoint");

    checkDistance("POINT (5 5)", "LINESTRING (5 1, 9 5)",
                  0,
                  "LINESTRING (5 5, 7 3)");
}

template<>
template<>
void object::test<3>()
{
    set_test_name("testZeroToleranceLine");

    checkDistance("LINESTRING (1 5, 5 5)", "LINESTRING (5 1, 9 5)",
                  0,
                  "LINESTRING (1 5, 5 1)");
}

template<>
template<>
void object::test<4>()
{
    set_test_name("testZeroToleranceZeroLengthLineQuery");

    checkDistance("LINESTRING (5 5, 5 5)", "LINESTRING (5 1, 9 5)",
                  0,
                  "LINESTRING (5 5, 7 3)");
}

template<>
template<>
void object::test<5>()
{
    set_test_name("testZeroLengthLineQuery");

    checkDistance("LINESTRING (5 5, 5 5)", "LINESTRING (5 1, 9 5)",
                  "LINESTRING (5 5, 7 3)");
}

template<>
template<>
void object::test<6>()
{
    set_test_name("testZeroLengthPolygonQuery");

    checkDistance("POLYGON ((5 5, 5 5, 5 5, 5 5))", "LINESTRING (5 1, 9 5)",
                  "LINESTRING (5 5, 7 3)");
}

template<>
template<>
void object::test<7>()
{
    set_test_name("testZeroLengthLineTarget");

    checkDistance("POINT (5 5)", "LINESTRING (5 1, 5 1)",
                  "LINESTRING (5 5, 5 1)");
}

template<>
template<>
void object::test<8>()
{
    set_test_name("testNegativeTolerancePoint");

    auto g1 = read("POINT (5 5)");
    auto g2 = read("LINESTRING (5 1, 9 5)");

    ensure_THROW(DirectedHausdorffDistance::distance(*g1, *g2, -1.0), geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<9>()
{
    set_test_name("testNegativeToleranceLine");

    auto g1 = read("LINESTRING (1 1, 5 5)");
    auto g2 = read("LINESTRING (5 1, 9 5)");

    ensure_THROW(DirectedHausdorffDistance::distance(*g1, *g2, -1.0), geos::util::IllegalArgumentException);
}

template<>
template<>
void object::test<10>()
{
    set_test_name("testPointPoint");

    checkHausdorff("POINT (0 0)", "POINT (1 1)", "LINESTRING (0 0, 1 1)");
}

template<>
template<>
void object::test<11>()
{
    set_test_name("testPointPoints");

    std::string a = "MULTIPOINT ((0 1), (2 3), (4 5), (6 6))";
    std::string b = "MULTIPOINT ((0.1 0), (1 0), (2 0), (3 0), (4 0), (5 0))";
    checkDistance(a, b, "LINESTRING (6 6, 5 0)");
    checkDistance(b, a, "LINESTRING (5 0, 2 3)");
    checkHausdorff(a, b, "LINESTRING (6 6, 5 0)");
}

template<>
template<>
void object::test<12>()
{
    set_test_name("testPointPolygonInterior");

    checkDistance("POINT (3 4)", "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))",
                  0);
}

template<>
template<>
void object::test<13>()
{
    set_test_name("testPointsPolygon");

    checkDistance("MULTIPOINT ((4 3), (2 8), (8 5))", "POLYGON ((6 9, 6 4, 9 1, 1 1, 6 9))",
                  "LINESTRING (2 8, 4.426966292134832 6.48314606741573)");
}

template<>
template<>
void object::test<14>()
{
    set_test_name("testLineSegments");

    checkHausdorff("LINESTRING (0 0, 2 0)", "LINESTRING (0 0, 2 1)",
                   "LINESTRING (2 0, 2 1)");
}

template<>
template<>
void object::test<15>()
{
    set_test_name("testLineSegments2");

    checkHausdorff("LINESTRING (0 0, 2 0)", "LINESTRING (0 1, 1 2, 2 1)",
                   "LINESTRING (1 0, 1 2)");
}

template<>
template<>
void object::test<16>()
{
    set_test_name("testLinePoints");

    checkHausdorff("LINESTRING (0 0, 2 0)", "MULTIPOINT (0 2, 1 0, 2 1)",
                   "LINESTRING (0 0, 0 2)");
}

template<>
template<>
void object::test<17>()
{
    set_test_name("testLinesTopoEqual");

    checkDistance("MULTILINESTRING ((10 10, 10 90, 40 30), (40 30, 60 80, 90 30, 40 10))",
                  "LINESTRING (10 10, 10 90, 40 30, 60 80, 90 30, 40 10)",
                  0.0);
}

template<>
template<>
void object::test<18>()
{
    set_test_name("testLinesPolygon");

    checkHausdorff("MULTILINESTRING ((1 1, 2 7), (7 1, 9 9))",
                   "POLYGON ((3 7, 6 7, 6 4, 3 4, 3 7))",
                   "LINESTRING (9 9, 6 7)");
}

template<>
template<>
void object::test<19>()
{
    set_test_name("testLinesPolygon2");

    std::string a = "MULTILINESTRING ((2 3, 2 7), (9 1, 9 8, 4 9))";
    std::string b = "POLYGON ((3 7, 6 8, 8 2, 3 4, 3 7))";
    checkDistance(a, b, "LINESTRING (9 8, 6.3 7.1)");
    checkHausdorff(a, b, "LINESTRING (2 3, 5.5 3)");
}

template<>
template<>
void object::test<20>()
{
    set_test_name("testPolygonLineCrossingBoundaryResult");

    checkDistance("POLYGON ((2 8, 8 2, 2 1, 2 8))",
                  "LINESTRING (6 5, 4 7, 0 0, 8 4)",
                  "LINESTRING (2 8, 3.9384615384615387 6.892307692307693)");
}

template<>
template<>
void object::test<21>()
{
    set_test_name("testPolygonLineCrossingInteriorPoint");

    checkDistanceStartPtLen("POLYGON ((2 8, 8 2, 2 1, 2 8))",
                            "LINESTRING (6 5, 4 7, 0 0, 9 1)",
                            "LINESTRING (4.555 2.989, 4.828 0.536)", 0.01);
}

template<>
template<>
void object::test<22>()
{
    set_test_name("testPolygonPolygon");

    std::string a = "POLYGON ((2 18, 18 18, 17 3, 2 2, 2 18))";
    std::string b = "POLYGON ((1 19, 5 12, 5 3, 14 10, 11 19, 19 19, 20 0, 1 1, 1 19))";
    checkDistance(b, a, "LINESTRING (20 0, 17 3)");
    checkDistance(a, b, "LINESTRING (6.6796875 18, 11 19)");
    checkHausdorff(a, b, "LINESTRING (6.6796875 18, 11 19)");
}

template<>
template<>
void object::test<23>()
{
    set_test_name("testPolygonPolygonHolesNested");

    // B is contained in A
    std::string a = "POLYGON ((1 19, 19 19, 19 1, 1 1, 1 19), (6 8, 11 14, 15 7, 6 8))";
    std::string b = "POLYGON ((2 18, 18 18, 18 2, 2 2, 2 18), (10 17, 3 7, 17 5, 10 17))";
    checkDistance(a, b, "LINESTRING (9.817138671875 12.58056640625, 7.863620425230705 13.948029178901006)");
    checkDistance(b, a, 0.0);
}

template<>
template<>
void object::test<24>()
{
    set_test_name("testMultiPolygons");

    std::string a = "MULTIPOLYGON (((1 1, 1 10, 5 1, 1 1)), ((4 17, 9 15, 9 6, 4 17)))";
    std::string b = "MULTIPOLYGON (((1 12, 4 13, 8 10, 1 12)), ((3 8, 7 7, 6 2, 3 8)))";
    checkDistance(a, b, "LINESTRING (1 1, 5.4 3.2)");
    checkDistanceStartPtLen(b, a,
                            "LINESTRING (2.669921875 12.556640625, 5.446115154109589 13.818546660958905)",
                            0.01);
}

template<>
template<>
void object::test<25>()
{
    set_test_name("testLinePolygonCrossing");

    std::string wkt1 = "LINESTRING (2 5, 5 10, 6 4)";
    std::string wkt2 = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))";
    checkDistance(wkt1, wkt2, "LINESTRING (5 10, 5 9)");
}

template<>
template<>
void object::test<26>()
{
    set_test_name("testNonVertexResult");

    std::string wkt1 = "LINESTRING (1 1, 5 10, 9 1)";
    std::string wkt2 = "LINESTRING (0 10, 0 0, 10 0)";

    checkHausdorff(wkt1, wkt2, "LINESTRING (6.53857421875 6.5382080078125, 6.53857421875 0)");
    checkDistance(wkt1, wkt2, "LINESTRING (6.53857421875 6.5382080078125, 6.53857421875 0)");
}

template<>
template<>
void object::test<27>()
{
    set_test_name("testDirectedLines");

    std::string wkt1 = "LINESTRING (1 6, 3 5, 1 4)";
    std::string wkt2 = "LINESTRING (1 10, 9 5, 1 2)";
    checkDistance(wkt1, wkt2, "LINESTRING (1 6, 2.797752808988764 8.876404494382022)");
    checkDistance(wkt2, wkt1, "LINESTRING (9 5, 3 5)");
}

template<>
template<>
void object::test<28>()
{
    set_test_name("testDirectedLines2");

    std::string wkt1 = "LINESTRING (1 6, 3 5, 1 4)";
    std::string wkt2 = "LINESTRING (1 3, 1 9, 9 5, 1 1)";
    checkDistance(wkt1, wkt2, "LINESTRING (3 5, 1 5)");
    checkDistance(wkt2, wkt1, "LINESTRING (9 5, 3 5)");
}

/**
 * Tests that segments are detected as interior even for a large tolerance.
 */
template<>
template<>
void object::test<29>()
{
    set_test_name("testInteriorSegmentsLargeTol");

    std::string a = "POLYGON ((4 6, 5 6, 5 5, 4 5, 4 6))";
    std::string b = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))";
    checkDistance(a, b, 2.0, 0.0);
}

/**
 * Tests that segment endpoint nearest points
 * which are interior to B have distance 0
 */
template<>
template<>
void object::test<30>()
{
    set_test_name("testInteriorSegmentsSameExterior");
    std::string a = "POLYGON ((1 9, 3 9, 4 5, 5.05 9, 9 9, 9 1, 1 1, 1 9))";
    std::string b = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9))";
    checkDistance(a, b, 0.0);
}

//-----------------------------------------------------

template<>
template<>
void object::test<31>()
{
    set_test_name("testFullyWithinDistanceEmptyPoints");
    std::string a = "POINT EMPTY";
    std::string b = "MULTIPOINT ((1 1), (9 9))";
    checkFullyWithinDistanceEmpty(a, b);
}

template<>
template<>
void object::test<32>()
{
    set_test_name("testFullyWithinDistanceEmptyLine");
    std::string a = "LINESTRING EMPTY";
    std::string b = "LINESTRING (9 9, 1 1)";
    checkFullyWithinDistanceEmpty(a, b);
}

template<>
template<>
void object::test<33>()
{
    //-- shows withinDistance envelope check not triggering for disconnected A
    set_test_name("testFullyWithinDistancePoints");
    std::string a = "MULTIPOINT ((1 9), (9 1))";
    std::string b = "MULTIPOINT ((1 1), (9 9))";
    checkFullyWithinDistance(a, b, 1, false);
    checkFullyWithinDistance(a, b, 8.1, true);
}

template<>
template<>
void object::test<34>()
{
    set_test_name("testFullyWithinDistanceDisconnectedLines");
    std::string a = "MULTILINESTRING ((1 9, 2 9), (8 1, 9 1))";
    std::string b = "LINESTRING (9 9, 1 1)";
    checkFullyWithinDistance(a, b, 1, false);
    checkFullyWithinDistance(a, b, 6, true);
    checkFullyWithinDistance(b, a, 1, false);
    checkFullyWithinDistance(b, a, 7.1, true);
}

template<>
template<>
void object::test<35>()
{
    set_test_name("testFullyWithinDistanceDisconnectedPolygons");
    std::string a = "MULTIPOLYGON (((1 9, 2 9, 2 8, 1 8, 1 9)), ((8 2, 9 2, 9 1, 8 1, 8 2)))";
    std::string b = "POLYGON ((1 2, 9 9, 2 1, 1 2))";
    checkFullyWithinDistance(a, b, 1, false);
    checkFullyWithinDistance(a, b, 5.3, true);
    checkFullyWithinDistance(b, a, 1, false);
    checkFullyWithinDistance(b, a, 7.1, true);
}

template<>
template<>
void object::test<36>()
{
    set_test_name("testFullyWithinDistanceLines");
    std::string a = "MULTILINESTRING ((1 1, 3 3), (7 7, 9 9))";
    std::string b = "MULTILINESTRING ((1 9, 1 5), (6 4, 8 2))";
    checkFullyWithinDistance(a, b, 1, false);
    checkFullyWithinDistance(a, b, 4, false);
    checkFullyWithinDistance(a, b, 6, true);
}

template<>
template<>
void object::test<37>()
{
    set_test_name("testFullyWithinDistancePolygons");
    std::string a = "POLYGON ((1 4, 4 4, 4 1, 1 1, 1 4))";
    std::string b = "POLYGON ((10 10, 10 15, 15 15, 15 10, 10 10))";
    checkFullyWithinDistance(a, b, 5, false);
    checkFullyWithinDistance(a, b, 10, false);
    checkFullyWithinDistance(a, b, 20, true);
}

template<>
template<>
void object::test<38>()
{
    set_test_name("testFullyWithinDistancePolygonsNestedWithHole");
    std::string a = "POLYGON ((2 8, 8 8, 8 2, 2 2, 2 8))";
    std::string b = "POLYGON ((1 9, 9 9, 9 1, 1 1, 1 9), (3 7, 7 7, 7 3, 3 3, 3 7))";
    checkFullyWithinDistance(a, b, 1, false);
    checkFullyWithinDistance(a, b, 2, true);
    checkFullyWithinDistance(a, b, 3, true);
}

template<>
template<>
void object::test<39>()
{
    set_test_name("CircularString/LineString");

    auto g1 = reader.read("CIRCULARSTRING (-5 0, 0 5, 0 -5)");
    auto g2 = reader.read("LINESTRING (-8 -2, -2 -6)");

    auto pts = DirectedHausdorffDistance::distancePoints(*g1, *g2);
    ensure(pts.has_value());
    std::cout << "LINESTRING (" << pts->at(0) << ", " << pts->at(1) << ")" << std::endl;

    pts = DirectedHausdorffDistance::distancePoints(*g2, *g1);
    ensure(pts.has_value());
    std::cout << "LINESTRING (" << pts->at(0) << ", " << pts->at(1) << ")" << std::endl;
}

template<>
template<>
void object::test<40>()
{
    set_test_name("CircularString / linearized approximation");

    auto g1 = reader.read("CIRCULARSTRING (-5 4, 13 5, 6 7, 5 -20, -5 4)");
    auto g2 = reader.read("LINESTRING(-5 4,-7.612791972449398 1.807191851097221,-9.747584231894455 -0.853212531361182,-11.32233790900999 -3.878975298929304,-12.276536103116497 -7.153818085316818,-12.573509512734969 -10.551890518795839,-12.201845616278373 -13.942606581408736,-11.175827248856024 -17.195662956404778,-9.534883720930234 -20.186046511627904,-7.342075572027454 -22.798838484077304,-4.681671189569057 -24.933630743522357,-1.655908422000929 -26.5083844206379,1.618934364386585 -27.462582614744406,5.017006797865607 -27.759556024362876,8.407722860478511 -27.387892127906277,11.660779235474545 -26.361873760483935,14.651162790697672 -24.720930232558143,17.263954763147073 -22.52812208365536,19.398747022592126 -19.867717701196966,20.973500699707667 -16.841954933628838,21.927698893814174 -13.567112147241323,22.224672303432644 -10.169039713762302,21.853008406976045 -6.778323651149398,20.826990039553703 -3.525267276153363,19.18604651162791 -0.534883720930234,16.99323836272513 2.077908251519164,14.332833980266734 4.212700510964218,11.307071212698606 5.787454188079759,8.032228426311093 6.741652382186266,6 7,8.577976384455235 6.303885187286109,10.970612412583517 5.118207782512325,13.085960423629622 3.488532703408558,14.842728779829878 1.477487449332455,16.173405858462555 -0.837644638409017,17.0268544856481 -3.367894332024982,17.370277109263284 -6.016025554770154,17.19047619047619 -8.680272108843537,16.4943613777623 -11.25824849329877,15.308683972988517 -13.650884521427054,13.67900889388475 -15.766232532473158,11.667963639808645 -17.523000888673415,9.352831552067176 -18.853677967306094,6.82258185845121 -19.70712659449164,4.174450635706037 -20.05054921810682,1.510204081632655 -19.87074829931973,-1.06777230282258 -19.17463348660584,-3.460408330950862 -17.988956081832054,-5.575756341996966 -16.35928100272829,-7.332524698197224 -14.348235748652183,-8.663201776829903 -12.033103660910712,-9.516650404015449 -9.50285396729475,-9.86007302763063 -6.854722744549579,-9.680272108843537 -4.190476190476193,-8.984157296129645 -1.612499806020957,-7.798479891355864 0.780136222107323,-6.168804812252101 2.895484233153426,-5 4)");

    auto pts = DirectedHausdorffDistance::distancePoints(*g1, *g2);
    ensure(pts.has_value());
    std::cout << "LINESTRING (" << pts->at(0) << ", " << pts->at(1) << ")" << std::endl;

    pts = DirectedHausdorffDistance::distancePoints(*g2, *g1);
    ensure(pts.has_value());
    std::cout << "LINESTRING (" << pts->at(0) << ", " << pts->at(1) << ")" << std::endl;
}

} // namespace tut
