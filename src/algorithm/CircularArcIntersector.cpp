/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2024 ISciences, LLC
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/algorithm/Angle.h>
#include <geos/algorithm/CircularArcIntersector.h>
#include <geos/algorithm/LineIntersector.h>
#include <geos/geom/CoordinateSequences.h>

#include <iomanip>

using geos::geom::CoordinateSequence;
using geos::geom::CoordinateXY;
using geos::geom::CircularArc;

namespace geos::algorithm {

static double
nextAngleCCW(double from, double a, double b)
{
    if (Angle::normalizePositive(a - from) < Angle::normalizePositive(b - from)) {
        return a;
    }
    else {
        return b;
    }
}

static double
angleFractionCCW(double a, double b, double x)
{
    if (x < a) {
        x += 2*MATH_PI;
    }
    if (b < a) {
        b += 2*MATH_PI;
    }
    return (x - a) / (b - a);
}

static void interpolateZM(const CircularArc& arc,
                   const CoordinateXY& pt,
                   double& z, double& m)
{
    using geos::geom::Ordinate;

    const geom::CoordinateSequence& seq = *arc.getCoordinateSequence();
    std::size_t i0 = arc.getCoordinatePosition();

    if (arc.p1().equals2D(pt)) {
        seq.applyAt(i0 + 1, [&z, &m](const auto& arcPt) {
            z = arcPt.template get<Ordinate::Z>();
            m = arcPt.template get<Ordinate::M>();
        });
        return;
    }

    double z0, m0;
    seq.applyAt(i0, [&z0, &m0](const auto& arcPt) {
        z0 = arcPt.template get<Ordinate::Z>();
        m0 = arcPt.template get<Ordinate::M>();
    });
    if (arc.p0().equals2D(pt)) {
        z = z0;
        m = m0;
        return;
    }

    double z2, m2;
    seq.applyAt(i0 + 2, [&z2, &m2](const auto& arcPt) {
        z2 = arcPt.template get<Ordinate::Z>();
        m2 = arcPt.template get<Ordinate::M>();
    });

    if (arc.p2().equals2D(pt)) {
        z = z2;
        m = m2;
        return;
    }

    // FIXME should probably interpolate between pt0/pt1 or pt1/pt2 rather than pt0/pt2
    const double theta0 = arc.theta0();
    const double theta2 = arc.theta2();
    const double theta = CircularArcs::getAngle(pt, arc.getCenter());

    const double frac = arc.isCCW() ? angleFractionCCW(theta0, theta2, theta) : angleFractionCCW(theta2, theta0, theta);

    z = frac*(z0 + z2);
    m = frac*(m0 + m2);
}

static void interpolateZM(const CoordinateSequence& seq,
                              std::size_t ind0, std::size_t ind1,
                              geom::CoordinateXY& pt, double& z, double& m)
{
    seq.applyAt(ind0, [&seq, &pt, ind1, &z, &m](const auto& p0) {
        // FIXME should probably interpolate between pt0/pt1 or pt1/pt2 rather than pt0/pt2
        const auto& p1 = seq.getAt<std::decay_t<decltype(p0)>>(ind1);

        z = Interpolate::zGetOrInterpolate(pt, p0, p1);
        m = Interpolate::mGetOrInterpolate(pt, p0, p1);
    });
}


static void interpolateZM(const CircularArc& arc0,
                   const CircularArc& arc1,
                   geom::CoordinateXYZM& pt)
{
    double z0, m0;
    double z1, m1;
    interpolateZM(arc0, pt, z0, m0);
    interpolateZM(arc1, pt, z1, m1);

    // FIXME are these semantics (getOrAverage) what we use for linear intersections?
    pt.z = Interpolate::getOrAverage(z0, z1);
    pt.m = Interpolate::getOrAverage(m0, m1);
}

static void interpolateZM(const CircularArc& arc0,
                       const geom::CoordinateSequence& seq,
                        std::size_t ind0, std::size_t ind1,
                       geom::CoordinateXYZM& pt)
{
    double z0, m0;
    double z1, m1;
    interpolateZM(arc0, pt, z0, m0);
    interpolateZM(seq, ind0, ind1, pt, z1, m1);

    // FIXME are these semantics (getOrAverage) what we use for linear intersections?
    pt.z = Interpolate::getOrAverage(z0, z1);
    pt.m = Interpolate::getOrAverage(m0, m1);
}

#if 0
class DoIntersectArcArc {
public:
    DoIntersectArcArc(CircularArcIntersector& cai,
                const CircularArc& arc0,
                const CircularArc& arc1) :
        m_cai(cai),
        m_arc0(arc0),
        m_arc1(arc1)
    {}

    template<typename C1, typename C2>
    void operator()()
    {
        const CoordinateSequence& seq0 = *m_arc0.getCoordinateSequence();
        const auto pos0 = m_arc0.getCoordinatePosition();

        const CoordinateSequence& seq1 = *m_arc1.getCoordinateSequence();
        const auto pos1 = m_arc1.getCoordinatePosition();


        const C1& p0 = seq0.getAt<C1>(pos0);
        const C1& p2 = seq0.getAt<C1>(pos0 + 2);

        const C2& q0 = seq1.getAt<C2>(pos1);
        const C2& q2 = seq1.getAt<C2>(pos1 + 2);

        m_cai.computeIntersectionArcArc(p0, p2, m_arc0.getCenter(), m_arc0.getRadius(), m_arc0.getOrientation(),
                                        q0, q2, m_arc1.getCenter(), m_arc0.getRadius(), m_arc0.getOrientation());
    }
private:
    CircularArcIntersector& m_cai;
    const CircularArc& m_arc0;
    const CircularArc& m_arc1;
};

void computeIntersection(const CircularArc& arc0, const CircularArc& arc1)
{

    const geom::CoordinateSequence& seq0 = *arc0.getCoordinateSequence();
    const geom::CoordinateSequence& seq1 = *arc1.getCoordinateSequence();

    DoIntersectArcArc di;
    geom::CoordinateSequences::binaryDispatch(seq0, seq1, di);

}
#endif

int
CircularArcIntersector::circleIntersects(const CoordinateXY& center,
                                         double r,
                                         const CoordinateXY& p0,
                                         const CoordinateXY& p1,
                                         CoordinateXY& ret0,
                                         CoordinateXY& ret1)
{
    const double& x0 = center.x;
    const double& y0 = center.y;

    Envelope segEnv(p0, p1);

    CoordinateXY isect0, isect1;
    int n = 0;

    if (p1.x == p0.x) {
        // vertical line
        double x = p1.x;

        double A = 1;
        double B = 2*y0;
        double C = x*x - 2*x*x0 + x0*x0 + y0*y0 - r*r;

        double d = std::sqrt(B*B - 4*A*C);
        double Y1 = (-B + d)/(2*A);
        double Y2 = (-B - d)/(2*A);

        isect0 = {x, Y1};
        isect1 = {x, Y2};
    }
    else {
        double m = (p1.y - p0.y) / (p1.x - p0.x);
        double b = p1.y - p1.x*m;

        // Ax^2 + Bx + C = 0
        double A = 1 + m*m;
        double B = -2*x0 + 2*m*b - 2*m*y0;
        double C = x0*x0 + b*b - 2*b*y0 + y0*y0 - r*r;

        double d = std::sqrt(B*B - 4*A*C);
        double X1 = (-B + d)/(2*A);
        double X2 = (-B - d)/(2*A);

        // TODO use robust quadratic solver such as https://github.com/archermarx/quadratic ?
        // auto [X1, X2] = quadratic::solve(A, B, C);

        isect0 = {X1, m* X1 + b};
        isect1 = {X2, m* X2 + b};
    }

    if (segEnv.intersects(isect0)) {
        ret0 = isect0;
        if (segEnv.intersects(isect1) && !isect1.equals2D(isect0)) {
            ret1 = isect1;
            n = 2;
        } else {
            n = 1;
        }
    } else if (segEnv.intersects(isect1)) {
        ret0 = isect1;
        n = 1;
    }

    return n;
}

void
CircularArcIntersector::intersects(const CircularArc& arc, const CoordinateSequence& seq, std::size_t segPos0, std::size_t segPos1)
//CircularArcIntersector::intersects(const CircularArc& arc, const CoordinateXY& p0, const CoordinateXY& p1)
{
    if (arc.isLinear()) {
        seq.applyAt(segPos0, [this, &seq, &arc, segPos1](const auto& p0) {
            const auto& p1 = seq.getAt<std::decay_t<decltype(p0)>>(segPos1);
            intersects(arc.p0(), arc.p2(), p0, p1);
        });
        return;
    }

    reset();

    // TODO: envelope check?
    const CoordinateXY& c = arc.getCenter();
    const double r = arc.getRadius();

    CoordinateXYZM isect0, isect1;
    auto n = circleIntersects(c, r, seq.getAt<CoordinateXY>(segPos0), seq.getAt<CoordinateXY>(segPos1), isect0, isect1);

    if (n > 0 && arc.containsPointOnCircle(isect0)) {
        intPt[nPt++] = isect0;
        interpolateZM(arc,seq, segPos0, segPos1, intPt[nPt - 1]);
    }

    if (n > 1  && arc.containsPointOnCircle(isect1)) {
        intPt[nPt++] = isect1;
        interpolateZM(arc,seq, segPos0, segPos1, intPt[nPt - 1]);
    }

    switch (nPt) {
    case 2:
        result = TWO_POINT_INTERSECTION;
        break;
    case 1:
        result = ONE_POINT_INTERSECTION;
        break;
    default:
        result = NO_INTERSECTION;
    }
}

#if 0
void
CircularArcIntersector::intersects(const CoordinateXY& p0, const CoordinateXY& p1,
                                   const CoordinateXY& q0, const CoordinateXY& q1)
{
}
#endif

void CircularArcIntersector::intersects(const CircularArc& arc1, const CircularArc& arc2)
{
    // Handle cases where one or both arcs are degenerate
    if (arc1.isLinear()) {
        if (arc2.isLinear()) {
            // FIXME Z/M
            intersects(arc1.p0(), arc1.p2(), arc2.p0(), arc2.p2());
            return;
        } else {
            intersects(arc2, *arc1.getCoordinateSequence(), arc1.getCoordinatePosition(), arc1.getCoordinatePosition() + 2);
            return;
        }
    } else if (arc2.isLinear()) {
        intersects(arc1, *arc2.getCoordinateSequence(), arc2.getCoordinatePosition(), arc2.getCoordinatePosition() + 2);
        return;
    }

    reset();

    const auto& c1 = arc1.getCenter();
    const auto& c2 = arc2.getCenter();

    const auto r1 = arc1.getRadius();
    const auto r2 = arc2.getRadius();

    auto d = c1.distance(c2);

    if (d > r1 + r2) {
        // Circles are disjoint
        result = NO_INTERSECTION;
        return;
    }

    if (d < std::abs(r1-r2)) {
        // One circle contained within the other; arcs cannot intersect
        result = NO_INTERSECTION;
        return;
    }

    // a: the distance from c1 to the "radical line", which connects the two intersection points
    // The following expression was rewritten by
    double a = (d*d + r1*r1 - r2*r2) / (2*d);
    // Expression rewritten by Herbie, https://herbie.uwplse.org/demo/
    // double a = std::fma(r1-r2, (r1 + r2) / (d+d), d*0.5);

    // TODO because the circle center calculation is inexact we need some kind of tolerance here.
    // Take a PrecisionModel like LineIntersector?
    if (a == 0 || (d == 0 && r1 == r2)) {
        // COCIRCULAR

        double ap0 = arc1.theta0();
        double ap1 = arc1.theta2();
        double bp0 = arc2.theta0();
        double bp1 = arc2.theta2();

        bool resultArcIsCCW = true;

        if (arc1.getOrientation() != Orientation::COUNTERCLOCKWISE) {
            std::swap(ap0, ap1);
            resultArcIsCCW = false;
        }
        if (arc2.getOrientation() != Orientation::COUNTERCLOCKWISE) {
            std::swap(bp0, bp1);
        }
        ap0 = Angle::normalizePositive(ap0);
        ap1 = Angle::normalizePositive(ap1);
        bp0 = Angle::normalizePositive(bp0);
        bp1 = Angle::normalizePositive(bp1);

        bool checkBp1inA = true;

        // check start of B within A?
        if (Angle::isWithinCCW(bp0, ap0, ap1)) {
            double start = bp0;
            double end = nextAngleCCW(start, bp1, ap1);

            if (end == bp1) {
                checkBp1inA = false;
            }

            if (start == end) {
                intPt[nPt++] = CircularArcs::createPoint(c1, r1, start);
            }
            else {
                if (resultArcIsCCW) {
                    intArc[nArc++] = CircularArc(start, end, c1, r1, Orientation::COUNTERCLOCKWISE);
                }
                else {
                    intArc[nArc++] = CircularArc(end, start, c1, r1, Orientation::CLOCKWISE);
                }
            }
        }

        if (checkBp1inA && Angle::isWithinCCW(bp1, ap0, ap1)) {
            // end of B within A?
            double start = ap0;
            double end = bp1;
            if (start == end) {
                intPt[nPt++] = CircularArcs::createPoint(c1, r1, start);
            }
            else {
                if (resultArcIsCCW) {
                    intArc[nArc++] = CircularArc(start, end, c1, r1, Orientation::COUNTERCLOCKWISE);
                }
                else {
                    intArc[nArc++] = CircularArc(end, start, c1, r1, Orientation::CLOCKWISE);
                }
            }
        }
    } else {
        // NOT COCIRCULAR

        double dx = c2.x-c1.x;
        double dy = c2.y-c1.y;

#if 1
        // point where a line between the two circle center points intersects
        // the radical line
        CoordinateXY p{c1.x + a* dx/d, c1.y+a* dy/d};

        // distance from p to the intersection points
        double h = std::sqrt(r1*r1 - a*a);

        CoordinateXY isect0{p.x + h* dy/d, p.y - h* dx/d };
        CoordinateXY isect1{p.x - h* dy/d, p.y + h* dx/d };

        if (arc1.containsPointOnCircle(isect0) && arc2.containsPointOnCircle(isect0)) {
            intPt[nPt++] = isect0;
            interpolateZM(arc1, arc2, intPt[nPt-1]);
        }
        if (!isect1.equals2D(isect0) && arc1.containsPointOnCircle(isect1) && arc2.containsPointOnCircle(isect1)) {
            intPt[nPt++] = isect1;
            interpolateZM(arc1, arc2, intPt[nPt-1]);
        }
#else
        // Alternate formulation.
        // Instead of calculating the intersection points and determining if they fall on the arc,
        // calculate the angles of the intersection points. If they fall on the arc, create intersection points
        // at those angles.

        double centerPointAngle = std::atan2(dy, dx);

        double arc1IntPtAngleDeviation = std::acos(a / r1);

        double a11 = Angle::normalize(centerPointAngle - arc1IntPtAngleDeviation);
        double a12 = Angle::normalize(centerPointAngle + arc1IntPtAngleDeviation);

        double b = d - a;
        double arc2IntPtAngleDeviation = std::acos(b / r2);

        double a21 = Angle::normalize(centerPointAngle + MATH_PI + arc2IntPtAngleDeviation);
        double a22 = Angle::normalize(centerPointAngle + MATH_PI - arc2IntPtAngleDeviation);

        if (arc1.containsAngle(a11) && arc2.containsAngle(a21)) {
            intPt[nPt++] = CircularArcs::createPoint(arc1.getCenter(), arc1.getRadius(), a11);
        }
        if (arc1.containsAngle(a12) && arc2.containsAngle(a22)) {
            intPt[nPt++] = CircularArcs::createPoint(arc1.getCenter(), arc1.getRadius(), a12);
            if (nPt == 2 && intPt[0].equals(intPt[1])) {
                nPt = 1;
            }
        }
#endif

#if 1
        // Add endpoint intersections missed due to precision issues.
        // TODO: Add some logic to prevent double-counting of endpoints. Ideally, the endpoint test would happen before
        // computing intersection points, so if there is an endpoint intersection we get the exact intersection point
        // instead of a computed one.
        if (nPt < 2 && arc1.p0().equals2D(arc2.p0()) && (nPt == 0 || !intPt[0].equals2D(arc1.p0()))) {
            intPt[nPt++] = arc1.p0();
            interpolateZM(arc1, arc2, intPt[nPt-1]);
        }
        if (nPt < 2 && arc1.p0().equals2D(arc2.p2()) && (nPt == 0 || !intPt[0].equals2D(arc1.p0()))) {
            intPt[nPt++] = arc1.p0();
            interpolateZM(arc1, arc2, intPt[nPt-1]);
        }
        if (nPt < 2 && arc1.p2().equals2D(arc2.p0()) && (nPt == 0 || !intPt[0].equals2D(arc1.p2()))) {
            intPt[nPt++] = arc1.p2();
            interpolateZM(arc1, arc2, intPt[nPt-1]);
        }
        if (nPt < 2 && arc1.p2().equals2D(arc2.p2()) && (nPt == 0 || !intPt[0].equals2D(arc1.p2()))) {
            intPt[nPt++] = arc1.p2();
            interpolateZM(arc1, arc2, intPt[nPt-1]);
        }
#endif
    }

    if (nArc) {
        result = COCIRCULAR_INTERSECTION;
    }
    else {
        switch (nPt) {
        case 2:
            result = TWO_POINT_INTERSECTION;
            break;
        case 1:
            result = ONE_POINT_INTERSECTION;
            break;
        case 0:
            result = NO_INTERSECTION;
            break;
        }
    }
}

}
