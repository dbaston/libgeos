/**********************************************************************
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.osgeo.org
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifdef _MSC_VER
# if defined(GEOS_DEBUG_MSVC_USE_VLD) && !defined(GEOS_TEST_USE_STACKWALKER)
#  include <vld.h>
# else
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
#  include <crtdbg.h>
# endif
#pragma warning(disable : 4127)
#endif

#include <geos/geom/Point.h>
#include <geos/geom/LineString.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/IntersectionMatrix.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/geom/prep/PreparedGeometry.h>
#include <geos/geom/prep/PreparedGeometryFactory.h>
#include <geos/geom/util/Densifier.h>
#include <geos/operation/overlay/OverlayOp.h>
#include <geos/operation/overlay/snap/GeometrySnapper.h>
#include <geos/operation/overlayng/OverlayNG.h>
#include <geos/operation/overlayng/OverlayNGRobust.h>
#include <geos/operation/overlayng/UnaryUnionNG.h>
#include <geos/operation/buffer/BufferBuilder.h>
#include <geos/operation/buffer/BufferParameters.h>
#include <geos/operation/buffer/BufferOp.h>
#include <geos/operation/polygonize/BuildArea.h>
#include <geos/operation/valid/MakeValid.h>
#include <geos/precision/MinimumClearance.h>
#include <geos/util.h>
#include <geos/util/Interrupt.h>
#include <geos/io/WKBReader.h>
#include <geos/io/WKBWriter.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/operation/relate/RelateOp.h>
#include <geos/operation/polygonize/Polygonizer.h>
#include <geos/operation/linemerge/LineMerger.h>
#include <geos/profiler.h>
#include <geos/unload.h>
#include <geos/operation/valid/IsValidOp.h>
#include "XMLTester.h"
#include "BufferResultMatcher.h"
#include "SingleSidedBufferResultMatcher.h"

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <memory>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include <signal.h>

#if defined(_MSC_VER) && defined(GEOS_TEST_USE_STACKWALKER)
#include <windows.h>
#include "Stackwalker.h"
#endif

using namespace geos;
using namespace geos::operation::polygonize;
using namespace geos::operation::linemerge;
using namespace geos::geom::prep;
using std::runtime_error;
using geos::operation::overlayng::OverlayNG;
using geos::operation::overlayng::UnaryUnionNG;
using geos::operation::overlayng::OverlayNGRobust;
using operation::valid::TopologyValidationError;

namespace {

std::unique_ptr<const PreparedGeometry>
prepare(const geom::Geometry* g)
{
    return PreparedGeometryFactory::prepare(g);
}

// Asymmetric Rounding Algorithm  - equivalent to Java Math.round()
// Copy from geos/util/math.cpp
double
java_math_round(double val)
{
    double n;
    double f = std::fabs(std::modf(val, &n));

    if(val >= 0) {
        if(f < 0.5) {
            return std::floor(val);
        }
        else if(f > 0.5) {
            return std::ceil(val);
        }
        else {
            return (n + 1.0);
        }
    }
    else {
        if(f < 0.5) {
            return std::ceil(val);
        }
        else if(f > 0.5) {
            return std::floor(val);
        }
        else {
            return n;
        }
    }
} // java_math_round


#ifdef not_used
// a utility function defining a very simple method to indent a line of text
const char*
getIndent(unsigned int numIndents)
{
    static const char* pINDENT = "                                      + ";
    static const unsigned int LENGTH = static_cast<unsigned int>(strlen(pINDENT));

    if(numIndents > LENGTH) {
        numIndents = LENGTH;
    }

    return &pINDENT[ LENGTH - numIndents ];
}
#endif


// void dump_to_stdout( const tinyxml2::XMLNode * pParent, unsigned int indent = 0 )
// {
//     if ( !pParent ) return;

//     printf( "%s", getIndent( indent));

//     tinyxml2::XMLPrinter printer;
//     pParent->Accept(&printer);
// }


}

void
tolower(std::string& str)
{
    std::transform(
        str.begin(),
        str.end(),
        str.begin(),
        [](char c){ return (char)std::tolower(c); }
        );
}

std::string
normalize_filename(const std::string& str)
{
    std::string newstring;

    std::string::size_type last_slash = str.find_last_of('/', str.size());
    if(last_slash == std::string::npos) {
        newstring = str;
    }
    else {
        newstring = str.substr(last_slash + 1);
    }

    for(std::string::iterator i = newstring.begin(), e = newstring.end(); i != e; ++i) {
        if(*i == '.') {
            *i = '_';
        }
    }

    tolower(newstring);

    return newstring;
}

static int
checkOverlaySuccess(geom::Geometry const& gRes, geom::Geometry const& gRealRes)
{
    double tol = operation::overlay::snap::GeometrySnapper::computeSizeBasedSnapTolerance(gRes);
    if(gRes.equals(&gRealRes)) {
        return 1;
    }
    std::cerr << "Using an overlay tolerance of " << tol << std::endl;
    if(gRes.equalsExact(&gRealRes, tol)) {
        return 1;
    }
    return 0;
}

/* Could be an XMLTester class private but oh well.. */
static int
checkBufferSuccess(geom::Geometry const& gRes, geom::Geometry const& gRealRes, double dist)
{

    using geos::xmltester::BufferResultMatcher;

    int success = 1;
    do {

        if(gRes.getGeometryTypeId() != gRealRes.getGeometryTypeId()) {
            std::cerr << "Expected result is of type "
                      << gRes.getGeometryType()
                      << "; obtained result is of type "
                      << gRealRes.getGeometryType()
                      << std::endl;
            success = 0;
            break;
        }

        // Is a buffer always an area ?
        if(gRes.getDimension() != 2) {
            std::cerr << "Don't know how to validate "
                      << "result of buffer operation "
                      << "when expected result is not an "
                      << "areal type."
                      << std::endl;
        }

        if(!BufferResultMatcher::isBufferResultMatch(gRealRes, gRes, dist)) {
            std::cerr << "BufferResultMatcher FAILED" << std::endl;
            success = 0;
            break;
        }

    }
    while(0);

    return success;
}

static int
checkSingleSidedBufferSuccess(geom::Geometry& gRes,
                              geom::Geometry& gRealRes, double dist)
{
    int success = 1;
    do {

        if(gRes.getGeometryTypeId() != gRealRes.getGeometryTypeId()) {
            std::cerr << "Expected result is of type "
                      << gRes.getGeometryType()
                      << "; obtained result is of type "
                      << gRealRes.getGeometryType()
                      << std::endl;
            success = 0;
            break;
        }

        geos::xmltester::SingleSidedBufferResultMatcher matcher;
        if(! matcher.isBufferResultMatch(gRealRes,
                                         gRes,
                                         dist)) {
            std::cerr << "SingleSidedBufferResultMatcher FAILED" << std::endl;
            success = 0;
            break;
        }

    }
    while(0);

    return success;
}

XMLTester::XMLTester()
    :
    gA(nullptr),
    gB(nullptr),
    pm(nullptr),
    factory(nullptr),
    wktreader(nullptr),
    wktwriter(nullptr),
    wkbreader(nullptr),
    wkbwriter(nullptr),
    test_predicates(0),
    failed(0),
    succeeded(0),
    caseCount(0),
    testCount(0),
    testFileCount(0),
    totalTestCount(0),
    curr_file(nullptr),
    testValidOutput(false),
    testValidInput(false),
    sqlOutput(false),
    HEXWKB_output(true)
{
    setVerbosityLevel(0);
}

int
XMLTester::setVerbosityLevel(int value)
{
    int old_value = verbose;

    verbose = value;

    return old_value;
}

/*private*/
void
XMLTester::printTest(bool success, const std::string& expected_result, const std::string& actual_result,
                     const util::Profile& prof)
{
    if(sqlOutput) {
        std::cout << "INSERT INTO \"" << normalize_filename(*curr_file) << "\" VALUES ("
                  << caseCount << ", "
                  << testCount << ", "
                  << "'" << opSignature << "', "
                  << "'" << curr_case_desc << "', ";

        std::string geomOut;

        if(gA) {
            std::cout << "'" << printGeom(gA) << "', ";
        }
        else {
            std::cout << "NULL, ";
        }

        if(gB) {
            std::cout << "'" << printGeom(gB) << "', ";
        }
        else {
            std::cout << "NULL, ";
        }

        std::cout << "'" << expected_result << "', "
                  << "'" << actual_result << "', ";

        if(success) {
            std::cout << "'t'";
        }
        else {
            std::cout << "'f'";
        }

        std::cout << ");" << std::endl;
    }

    else {
        std::cout << *curr_file << ":";
        std::cout << " case" << caseCount << ":";
        std::cout << " test" << testCount << ": "
                  << opSignature;
        std::cout << ": " << (success ? "ok." : "failed.");
        std::cout << " (" << std::setprecision(15) << java_math_round(prof.getTot() / 1000) << " ms)" << std::endl;

        // print geometry on failure for -v
        // print geometry no matter what for -v -v and above
        if (verbose > 1 || (verbose == 1 && !success)) {
            std::cout << "\tDescription: " << curr_case_desc << std::endl;

            if(gA) {
                std::cout << "\tGeometry A: ";
                printGeom(std::cout, gA);
                std::cout << std::endl;
            }

            if(gB) {
                std::cout << "\tGeometry B: ";
                printGeom(std::cout, gB);
                std::cout << std::endl;
            }

            std::cout << "\tExpected result: " << expected_result << std::endl;
            std::cout << "\tObtained result: " << actual_result << std::endl;
            std::cout << std::endl;
        }
    }
}

void
XMLTester::run(const std::string& source)
{
    curr_file = &source;

    if(sqlOutput) {
        std::cout << "CREATE TABLE \"" << normalize_filename(*curr_file) << "\""
                  << "( caseno integer, testno integer, "
                  << " operation varchar, description varchar, "
                  << " a geometry, b geometry, expected geometry, "
                  << " obtained geometry, result bool )"

                  // NOTE: qgis 0.7.4 require oids for proper operations.
                  //       The 'WITH OIDS' parameter is supported back to
                  //       PostgreSQL 7.2, so if you run an older version
                  //       rebuild with the next line commented out.
                  //<< " WITH OIDS"

                  << ";" << std::endl;
    }

    ++testFileCount;

    caseCount = 0;

    tinyxml2::XMLError e = xml.LoadFile(source.c_str());
    if(e) {
        std::stringstream err;
        err << "Could not load " << source << ": " << e << std::endl;
        xml.~XMLDocument(); // Deallocates various internal pools
        throw runtime_error(err.str());
    }

    //dump_to_stdout(&xml);

    const tinyxml2::XMLNode* node = xml.FirstChildElement("run");

    if(! node) {
        throw(runtime_error("Document has no childs"));
    }

    parseRun(node);

}

void
XMLTester::resultSummary(std::ostream& os) const
{
    os << "Files: " << testFileCount << std::endl;
    os << "Tests: " << totalTestCount << std::endl;
    os << "Failed: " << failed << std::endl;
    os << "Succeeded: " << succeeded << std::endl;
}

void
XMLTester::resetCounters()
{
    testFileCount = totalTestCount = failed = succeeded = 0;
}

void
XMLTester::parseRun(const tinyxml2::XMLNode* node)
{
    using geos::geom::PrecisionModel;

    assert(node);

    //dump_to_stdout(node);

    // Look for precisionModel element
    const tinyxml2::XMLElement* el = node->FirstChildElement("precisionModel");
    if(el) {
        parsePrecisionModel(el);
    }
    else {
        pm.reset(new PrecisionModel());
    }

    // Look for geometryOperation, if any
    usePrepared = false;
    el = node->FirstChildElement("geometryOperation");
    if(el) {
        const tinyxml2::XMLNode* txt = el->FirstChild();
        if(txt) {
            std::string op = trimBlanks(txt->Value());
            if(op.find("PreparedGeometryOperation")) {
                usePrepared = true;
            }
            else {
                std::cerr << *curr_file
                          << ": WARNING: unknown geometryOperation: "
                          << op << std::endl;
            }
        }
    }

    if(verbose > 1) {
        std::cerr << *curr_file << ": run: Precision Model: " << pm->toString();
        if(usePrepared) {
            std::cerr << " (prepared)";
        }
        std::cerr << std::endl;
    }


    factory = geom::GeometryFactory::create(pm.get());
    wktreader.reset(new io::WKTReader(factory.get()));
    wktwriter.reset(new io::WKTWriter());
    wktwriter->setTrim(true);
    wkbreader.reset(new io::WKBReader(*factory));
    wkbwriter.reset(new io::WKBWriter());

    const tinyxml2::XMLNode* casenode;
    for(casenode = node->FirstChildElement("case");
            casenode;
            casenode = casenode->NextSiblingElement("case")) {
        try {
            parseCase(casenode);
        }
        catch(const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
        }
    }

}

void
XMLTester::parsePrecisionModel(const tinyxml2::XMLElement* el)
{
    using geos::geom::PrecisionModel;

    //dump_to_stdout(el);

    /* This does not seem to work... */
    std::string type;
    const char* typeStr = el->Attribute("type");
    if(typeStr) {
        type = typeStr;
    }

    const char* scaleStr = el->Attribute("scale");

    if(! scaleStr) {
        if(type == "FLOATING_SINGLE") {
            pm.reset(new PrecisionModel(PrecisionModel::FLOATING_SINGLE));
        }
        else {
            pm.reset(new PrecisionModel());
        }
    }
    else {

        char* stopstring;

        double scale = std::strtod(scaleStr, &stopstring);
        double offsetX = 0;
        double offsetY = 2;

        if(! el->QueryDoubleAttribute("offsetx", &offsetX))
        {} // std::cerr << "No offsetx" << std::endl;

        if(! el->QueryDoubleAttribute("offsety", &offsetY))
        {} // std::cerr << "No offsety" << std::endl;

        // NOTE: PrecisionModel discards offsets anyway...
        pm.reset(new PrecisionModel(scale, offsetX, offsetY));
    }
}

bool
XMLTester::testValid(const geom::Geometry* g, const std::string& label)
{
    operation::valid::IsValidOp ivo(g);
    bool valid = ivo.isValid();
    if(! valid) {
        const TopologyValidationError* err = ivo.getValidationError();
        std::cerr << *curr_file << ":"
                  << " case" << caseCount << ":"
                  << " test" << testCount << ": "
                  << opSignature << ": "
                  << " invalid geometry (" << label
                  << "): " << err->toString() << std::endl;
    }
    return valid;
}

/**
 * Parse WKT or HEXWKB
 */
geom::Geometry*
XMLTester::parseGeometry(const std::string& in, const char* label)
{
    if((! wkbreader.get()) || (! wktreader.get())) {
        throw(runtime_error("No precision model specified"));
    }

    std::stringstream is(in, std::ios_base::in);
    char first_char;

    // Remove leading spaces
    while(is.get(first_char) && std::isspace(first_char));
    is.unget();

    std::unique_ptr<geom::Geometry> ret;

    switch(first_char) {
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    case 'A':
    case 'B':
    case 'C':
    case 'D':
    case 'E':
    case 'F':
        ret = wkbreader->readHEX(is);
        break;
    default:
        ret = wktreader->read(in);
        break;
    }

    if(testValidInput) {
        testValid(ret.get(), std::string(label));
    }

    //ret->normalize();

    return ret.release();
}

std::string
XMLTester::trimBlanks(const std::string& in)
{
    std::string out;
    std::string::size_type pos = in.find_first_not_of(" \t\n\r");
    if(pos != std::string::npos) {
        out = in.substr(pos);
    }
    pos = out.find_last_not_of(" \t\n\r");
    if(pos != std::string::npos) {
        out = out.substr(0, pos + 1);
    }
    return out;
}

void
XMLTester::parseCase(const tinyxml2::XMLNode* node)
{
    assert(node);

    std::string geomAin;
    std::string geomBin;
    std::string thrownException;

    gA = nullptr;
    gB = nullptr;


    //dump_to_stdout(node);

    curr_case_desc.clear();
    const tinyxml2::XMLNode* txt = node->FirstChildElement("desc");
    if(txt) {
        txt = txt->FirstChild();
        if(txt) {
            curr_case_desc = trimBlanks(txt->Value());
        }
    }

    //std::cerr << "Desc: " << curr_case_desc << std::endl;


    try {
        const tinyxml2::XMLNode* el = node->FirstChildElement("a");
        geomAin = el->FirstChild()->Value();
        geomAin = trimBlanks(geomAin);
        gA = parseGeometry(geomAin, "Geometry A");

        if(nullptr != (el = node->FirstChildElement("b"))) {
            geomBin = el->FirstChild()->Value();
            geomBin = trimBlanks(geomBin);
            gB = parseGeometry(geomBin, "Geometry B");
        }
    }
    catch(const std::exception& e) {
        thrownException = e.what();
    }
    catch(...) {
        thrownException = "Unknown exception";
    }

//std::cerr << "A: " << geomAin << std::endl;
//std::cerr << "B: " << geomBin << std::endl;


    if(thrownException != "") {
        std::cout << *curr_file << ":";
        std::cout << " case" << caseCount << ":";
        std::cout << " skipped (" << thrownException << ")." << std::endl;
        return;
    }

    ++caseCount;
    testCount = 0;

    const tinyxml2::XMLNode* testnode;
    for(testnode = node->FirstChildElement("test");
            testnode;
            testnode = testnode->NextSiblingElement("test")) {
        parseTest(testnode);
    }

    totalTestCount += testCount;

    delete gA;
    delete gB;
}

/*private*/
void
XMLTester::printGeom(std::ostream& os, const geom::Geometry* g)
{
    os << printGeom(g);
}

std::string
XMLTester::printGeom(const geom::Geometry* g)
{
    if(HEXWKB_output) {
        std::stringstream s(std::ios_base::binary | std::ios_base::in | std::ios_base::out);
        wkbwriter->write(*g, s);
        std::stringstream s2;
        wkbreader->printHEX(s, s2);
        return s2.str();
    }
    else {
        wktwriter->setRoundingPrecision(16);
        return wktwriter->write(g);
    }
}

/**
* Computes the maximum area delta value
* resulting from identity equations over the overlay operations.
* The delta value is normalized to the total area of the geometries.
* If the overlay operations are computed correctly
* the area delta is expected to be very small (e.g. < 1e-6).
*/
double
XMLTester::areaDelta(const geom::Geometry* a, const geom::Geometry* b, std::string& rsltMaxDiffOp, double maxDiff, std::stringstream& ss)
{
    double areaA = a == nullptr ? 0 : a->getArea();
    double areaB = b == nullptr ? 0 : b->getArea();

    // if an input is non-polygonal delta is 0
    if (areaA == 0 || areaB == 0)
      return 0;

    std::unique_ptr<geom::Geometry> geomU = OverlayNGRobust::Union(a, b);
    std::unique_ptr<geom::Geometry> geomI = OverlayNGRobust::Intersection(a, b);
    std::unique_ptr<geom::Geometry> geomDab = OverlayNGRobust::Difference(a, b);
    std::unique_ptr<geom::Geometry> geomDba = OverlayNGRobust::Difference(b, a);
    std::unique_ptr<geom::Geometry> geomSD = OverlayNGRobust::SymDifference(a, b);

    double areaU   = geomU->getArea();
    double areaI   = geomI->getArea();
    double areaDab = geomDab->getArea();
    double areaDba = geomDba->getArea();
    double areaSD  = geomSD->getArea();


    double maxDelta = 0;

    // & : intersection
    // - : difference
    // + : union
    // ^ : symdifference

    double delta = std::abs(areaA - areaI - areaDab);
    if (delta > maxDelta) {
        rsltMaxDiffOp = "A = ( A & B ) + ( A - B )";
        maxDelta = delta;
    }

    delta = std::abs(areaB - areaI - areaDba);
    if (delta > maxDelta) {
        rsltMaxDiffOp = "B = ( A & B ) + ( B - A )";
        maxDelta = delta;
    }

    delta = std::abs(areaDab + areaDba - areaSD);
    if (delta > maxDelta) {
        maxDelta = delta;
        rsltMaxDiffOp = "( A ^ B ) = ( A - B ) + ( B - A )";
    }

    delta = std::abs(areaI + areaSD - areaU);
    if (delta > maxDelta) {
        maxDelta = delta;
        rsltMaxDiffOp = "( A + B ) = ( A & B ) + ( A ^ B )";
    }

    delta = std::abs(areaU - areaI - areaDab - areaDba);
    if (delta > maxDelta) {
        maxDelta = delta;
        rsltMaxDiffOp = "( A + B ) = ( A & B ) + ( A - B ) + ( A - B )";
    }

    // normalize the area delta value
    double diffScore = maxDelta / (areaA + areaB);

    if (diffScore > maxDiff) {
        ss << std::endl << "A" << std::endl;
        ss << *a;
        ss << std::endl << "B" << std::endl;
        ss << *b;
        ss << std::endl << "geomU" << std::endl;
        ss << *geomU;
        ss << std::endl << "geomI" << std::endl;
        ss << *geomI;
        ss << std::endl << "geomDab" << std::endl;
        ss << *geomDab;
        ss << std::endl << "geomDba" << std::endl;
        ss << *geomDba;
        ss << std::endl << "geomSD" << std::endl;
        ss << *geomSD;
        ss << std::endl;
    }

    return diffScore;
}


class Value {
public:
    explicit Value(bool value) : m_type(BOOL), m_bool(value) {}
    explicit Value(double value) : m_type(DOUBLE), m_double(value) {}
    explicit Value(const std::string& value) : m_type(STRING), m_string(value) {}
    explicit Value(std::unique_ptr<geom::Geometry> value) : m_type(GEOMETRY), m_geometry(std::move(value)) {}

    enum Type {
        BOOL,
        DOUBLE,
        INT,
        STRING,
        GEOMETRY
    };

    Type getType() const {
        return m_type;
    }

    bool getBool() const {
        checkType(BOOL);
        return m_bool;
    }

    double getDouble() const {
        checkType(DOUBLE);
        return m_double;
    }

    double getInt() const {
        checkType(INT);
        return m_int;
    }

    const std::string& getString() const {
        checkType(STRING);
        return m_string;
    }

    const geom::Geometry* getGeometry() const {
        checkType(GEOMETRY);
        return m_geometry.get();
    }

    bool operator==(const Value& other) const {
        if (getType() != other.getType()) {
            return false;
        }
        switch (getType()) {
            case BOOL: return getBool() == other.getBool();
            case DOUBLE: return getDouble() == other.getDouble();
            case INT: return getInt() == other.getInt();
            case STRING: return getString() == other.getString();
            case GEOMETRY: return getGeometry()->equals(other.getGeometry());
        }
    }

private:

    void checkType(Type expected) const {
        if (m_type != expected) {
            throw std::runtime_error("Incorrect type access");
        }
    }

    Type m_type;
    std::string m_string;
    std::unique_ptr<geom::Geometry> m_geometry;
    double m_double;
    int m_int;
    bool m_bool;
};

using Result = Value;

class Args {

public:
    const std::string& get(std::size_t index) const {
        return m_args[index - 1];
    }

    void set(std::size_t index, const std::string& value) {
        m_args[index - 1] = value;
    }

    bool has(std::size_t index) const {
        return get(index) != "";
    }

    int getInt(std::size_t index) const {
        return std::atoi(get(index).c_str());
    }

    int getIntOr(std::size_t index, int default_value) const {
        return has(index) ? getInt(index) : default_value;
    }

    double getDouble(std::size_t index) const {
        return std::atof(get(index).c_str());
    }

    double getDoubleOr(std::size_t index, double default_value) const {
        return has(index) ? getDouble(index) : default_value;
    }

    const geom::Geometry* A() const {
        return swapGeomArgs() ? m_geoms[1] : m_geoms[0];
    }

    const geom::Geometry* B() const {
        return swapGeomArgs() ? m_geoms[0] : m_geoms[1];
    }

    const geom::prep::PreparedGeometry* pA() const {
        return swapGeomArgs() ? m_prep[1].get() : m_prep[0].get();
    }

    const geom::prep::PreparedGeometry* pB() const {
        return swapGeomArgs() ? m_prep[0].get() : m_prep[1].get();
    }

    bool usePrepared() const {
        return m_useprepared;
    }

    void setUsePrepared(bool value) {
        m_useprepared = value;
    }

private:
    bool swapGeomArgs() const {
        return (get(1) == "B" || get(1) == "b") && m_geoms[1] != nullptr;
    }

    void prepareGeometries() const {
        for (std::size_t i = 0; i < m_prep.size(); i++) {
            if (m_prep[i] == nullptr) {
                m_prep[i] = PreparedGeometryFactory::prepare(m_geoms[i]);
            }
        }
    }

    static constexpr std::size_t MAX_ARGS = 4;

    std::array<std::string, MAX_ARGS> m_args;
    std::array<const geom::Geometry*, 2> m_geoms;
    mutable std::array<std::unique_ptr<geom::prep::PreparedGeometry>, 2> m_prep;
    bool m_useprepared;
};

class Test {

public:
    const std::string& getName() const {
        return m_name;
    }

    void setName(const std::string& value) {
        m_name = tolower(XMLTester::trimBlanks(value));
    }

    void setExpected(const std::string& value) {
        m_expected = XMLTester::trimBlanks(value);
    }

    void setResult(bool value) {
        if (value) {
            m_actual = "true";
        } else {
            m_actual = "false";
        }
    }


    std::string getSignature() const {
        std::string opSig;

        for (const auto& opArg : m_args) {
            if (opArg != "") {
                if (opSig == "") {
                    opSig += ", ";
                }
                opSig += opArg;
            }
        }

        return getName() + "(" + opSig + ")";
    }

    bool isSuccess() const {
        return m_actual == m_expected;
    }

    const Args& getArgs() const {
        return m_args;
    }

    bool setUsePrepared(bool value) {
        m_args.setUsePrepared(value);
    }

private:
    std::string m_name;

    std::string m_expected;
    std::string m_actual;

    Args m_args;
};

class Fizz {
public:

    Fizz(std::function<Value(const Args&)> f) :
        m_run(f),
        m_test([](const Value& a, const Value& b) {
            return a == b;
        })
    {}

    Fizz(std::function<Value(const Args&)> f,
         std::function<bool(const Value&, const Value&)> t) :
        m_run(f),
        m_test(t)
    {}

    template<typename F>
    Fizz& operator=(F&& f) {
        m_run = std::move(f);
        return *this;
    }

private:
    std::function<Value(const Args&)> m_run;
    std::function<bool(const Value&, const Value&)> m_test;
};

std::map<std::string, std::function<Result(const Args&)>> getFunctions() {
    using operation::buffer::BufferBuilder;
    using operation::buffer::BufferOp;
    using operation::buffer::BufferParameters;

    //std::map<std::string, std::function<Result(const Args&)>> functions;
    std::map<std::string, Fizz> functions;

    functions["isvalid"] = [](const Args& args) {
        const auto* toTest = args.A();
        if ((args.get(1) == "B" || args.get(1) == "b") && args.B()) {
            toTest = args.B();
        }
        return Result(toTest->isValid());
    };

    functions["isSimple"] = [](const Args& args) {
        return Result(args.A()->isSimple());
    };

    // Overlay

    functions["intersection"] = [](const Args& args) {
        return Result(args.A()->intersection(args.B()));
    };

    functions["intersectionng"] = [](const Args& args) {
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::INTERSECTION));
    };

    functions["unionng"] = [](const Args& args) {
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::UNION));
    };

    functions["differenceng"] = [](const Args& args) {
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::DIFFERENCE));
    };

    functions["symdifferenceng"] = [](const Args& args) {
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::SYMDIFFERENCE));
    };

    functions["intersectionsr"] = [](const Args& args) {
        auto precision = args.getDoubleOr(3, 0.0);
        geom::PrecisionModel precMod(precision);
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::INTERSECTION, &precMod));
    };

    functions["intersectionsin"] = [](const Args& args) {
        auto precision = args.getDoubleOr(3, 0.0);
        geom::PrecisionModel precMod(precision); // never used?
        return Result(OverlayNGRobust::Intersection(args.A(), args.B()));
    };

    functions["unionsr"] = [](const Args& args) {
        auto precision = args.getDoubleOr(3, 0.0);
        geom::PrecisionModel precMod(precision);
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::UNION, &precMod));
    };

    functions["differencesr"] = [](const Args& args) {
        auto precision = args.getDoubleOr(3, 0.0);
        geom::PrecisionModel precMod(precision);
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::DIFFERENCE, &precMod));
    };

    functions["symdifferencesr"] = [](const Args& args) {
        auto precision = args.getDoubleOr(3, 0.0);
        geom::PrecisionModel precMod(precision);
        return Result(OverlayNG::overlay(args.A(), args.B(), OverlayNG::SYMDIFFERENCE, &precMod));
    };

    functions["union"] = [](const Args& args) {
        return Result(args.A()->Union(args.B()));
    };

    functions["difference"] = [](const Args& args) {
        return Result(args.A()->difference(args.B()));
    };

    functions["symdifference"] = [](const Args& args) {
        return Result(args.A()->symDifference(args.B()));
    };

    // Relate, predicates

    functions["relate"] = [](const Args& args) {
        auto im = args.A()->relate(args.B());
        assert(im.get());

        return Result(im->matches(args.get(3)));
    };

    functions["relatestring"] = [](const Args& args) {
        return Result(args.A()->relate(args.B())->toString());
    };


    functions["intersects"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->intersects(args.B()));
        }
        return Result(args.A()->intersects(args.B()));
    };

    functions["contains"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->contains(args.B()));
        }
        return Result(args.A()->contains(args.B()));
    };

    functions["overlaps"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->overlaps(args.B()));
        }
        return Result(args.A()->overlaps(args.B()));
    };

    functions["within"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->within(args.B()));
        }
        return Result(args.A()->within(args.B()));
    };

    functions["touches"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->touches(args.B()));
        }
        return Result(args.A()->touches(args.B()));
    };

    functions["crosses"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->crosses(args.B()));
        }
        return Result(args.A()->crosses(args.B()));
    };

    functions["disjoint"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->disjoint(args.B()));
        }
        return Result(args.A()->disjoint(args.B()));
    };

    functions["covers"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->covers(args.B()));
        }
        return Result(args.A()->covers(args.B()));
    };

    functions["coveredby"] = [](const Args& args) {
        if (args.usePrepared()) {
            return Result(args.pA()->coveredBy(args.B()));
        }
        return Result(args.A()->coveredBy(args.B()));
    };

    functions["equalstopo"] = [](const Args& args) {
        return Result(args.A()->equals(args.B()));
    };

    functions["equalsexact"] = [](const Args& args) {
        return Result(args.A()->equalsExact(args.B()));
    };

    functions["equalsnorm"] = [](const Args& args) {
        auto g1 = args.A()->clone();
        auto g2 = args.B()->clone();
        g1->normalize();
        g2->normalize();

        return Result(g1->equalsExact(g2.get()));
    };

    functions["iswithindistance"] = [](const Args& args) {
        double dist = args.getDouble(3);

        if (args.usePrepared()) {
            return Result(args.pA()->isWithinDistance(args.B(), dist));
        }
        return Result(args.A()->isWithinDistance(args.B(), dist));
    };

    functions["distance"] = [](const Args& args) {
        return Result(args.A()->distance(args.B()));
    };

    functions["minclearance"] = [](const Args& args) {
        precision::MinimumClearance mc(args.A());
    return Result(
                // Hack for Inf/1.7976931348623157E308 comparison
                std::min(
                    mc.getDistance(),
                    1.7976931348623157E308
                ));
    };

    functions["minclearanceline"] = [](const Args& args) {
        precision::MinimumClearance mc(args.A());
        return Result(mc.getLine());
    };

    // Construction

    functions["getboundary"] = [](const Args& args) {
        return Result(args.A()->getBoundary());
    };

    functions["getcentroid"] = [](const Args& args) {
        return Result(args.A()->getCentroid());
    };

    functions["densify"] = [](const Args& args) {
        geom::util::Densifier den(args.A());
        double distanceTolerance = args.getDouble(2);
        den.setDistanceTolerance(distanceTolerance);
        return Result(den.getResultGeometry());
    };

    functions["convexhull"] = [](const Args& args) {
        return Result(args.A()->convexHull());
    };

    functions["buffer"] = [](const Args& args) {
        double dist = args.getDouble(2);
        BufferParameters params;
        if (args.has(3)) {
            params.setQuadrantSegments(args.getInt(3));
        }

        BufferOp op(args.A(), params);
        return Result(op.getResultGeometry(dist));
    };

    functions["buffersinglesided"] = [](const Args& args) {
        double dist = args.getDouble(2);

        BufferParameters params;
        params.setJoinStyle(BufferParameters::JOIN_ROUND);
        if (args.has(3)) {
            params.setQuadrantSegments(args.getInt(3));
        }
        bool leftSide = args.get(4) != "right";

        BufferBuilder bufBuilder(params);

        return Result(bufBuilder.bufferLineSingleSided(args.A(), dist, leftSide));
    };

    functions["buffermitredjoin"] = [](const Args& args) {
        double dist = args.getDouble(2);

        BufferParameters params;
        params.setJoinStyle(BufferParameters::JOIN_MITRE);

        if (args.has(3)) {
            params.setQuadrantSegments(args.getInt(3));
        }

        BufferOp op(args.A(), params);
        return Result(op.getResultGeometry(dist));
    };

    functions["getinteriorpoint"] = [](const Args& args) {
        return Result(args.A()->getInteriorPoint());
    };

    functions["polygonize"] = [](const Args& args) {
        Polygonizer p;
        p.add(args.A());

        auto polys = p.getPolygons();
        return Result(args.A()->getFactory()->createGeometryCollection(std::move(polys)));
    };

    functions["buildarea"] = [](const Args& args) {
        return Result(BuildArea().build(args.A()));
    };

    functions["linemerge"] = [](const Args& args) {
        LineMerger merger;
        merger.add(args.A());

        auto lines = merger.getMergedLineStrings();
        return Result(args.A()->getFactory()->createGeometryCollection(std::move(lines)));
    };


    functions["makevalid"] = [](const Args& args) {
        return Result(operation::valid::MakeValid().build(args.A()));
    };

    functions["overlayareatest"] = [](const Args& args) {
        std::string maxDiffOp;
        std::stringstream ss;
        double maxDiff = 1e-6;

        double areaDiff = XMLTester::areaDelta(args.A(), args.B(), maxDiffOp, maxDiff, ss);

        return Result(areaDiff);
    };

    functions["unionlength"] = [](const Args& args) {
        auto unionResult = OverlayNGRobust::Union(args.A());
        return Result(unionResult->getLength());
    };

    functions["unionarea"] = Fizz(
                [](const Args& args) {
                    auto unionResult = OverlayNGRobust::Union(args.A());
                    return Result(unionResult->getArea());
    },
                [](const Value& actual, const Value& expected) {

    });

    functions["areatest"] = [](const Args& args) {
        double areaA = args.A()->getArea();
        double areaB = args.B()->getArea();
        double areaI = args.A()->intersection(args.B())->getArea();
        double areaDab = args.A()->difference(args.B())->getArea();
        double areaDba = args.B()->difference(args.A())->getArea();
        double areaSD = args.A()->symDifference(args.B())->getArea();
        double areaU = args.A()->Union(args.B())->getArea();

        double maxdiff = 0;
        std::string maxdiffop;
        // @ : symdifference
        // - : difference
        // + : union
        // ^ : intersection

        // A == ( A ^ B ) + ( A - B )
        double diff = std::fabs(areaA - areaI - areaDab);
        if(diff > maxdiff) {
            maxdiffop = "A == ( A ^ B ) + ( A - B )";
            maxdiff = diff;
        }

        // B == ( A ^ B ) + ( B - A )
        diff = std::fabs(areaB - areaI - areaDba);
        if(diff > maxdiff) {
            maxdiffop = "B == ( A ^ B ) + ( B - A )";
            maxdiff = diff;
        }

        //  ( A @ B ) == ( A - B ) + ( B - A )
        diff = std::fabs(areaDab + areaDba - areaSD);
        if(diff > maxdiff) {
            maxdiffop = "( A @ B ) == ( A - B ) + ( B - A )";
            maxdiff = diff;
        }

        //  ( A u B ) == ( A ^ B ) + ( A @ B )
        diff = std::fabs(areaI + areaSD - areaU);
        if(diff > maxdiff) {
            maxdiffop = "( A u B ) == ( A ^ B ) + ( A @ B )";
            maxdiff = diff;
        }

        return Result(maxdiff);
    };

    return functions;
}

void
XMLTester::parseTest(const tinyxml2::XMLNode* node)
{
    using namespace operation::overlay;

    typedef std::unique_ptr< geom::Geometry > GeomPtr;

    int success = 0; // no success by default

    Test op;

    ++testCount;

    const tinyxml2::XMLNode* opnode = node->FirstChildElement("op");
    if(! opnode) {
        throw(runtime_error("case has no op"));
    }

    //dump_to_stdout(opnode);

    const tinyxml2::XMLElement* opel = opnode->ToElement();

    op.setName(opel->Attribute("name"));
    op.setArg(1, opel->Attribute("arg1"));
    op.setArg(2, opel->Attribute("arg2"));
    op.setArg(3, opel->Attribute("arg3"));
    op.setArg(4, opel->Attribute("arg4"));

    op.setUsePrepared(usePrepared);

    const tinyxml2::XMLNode* resnode = opnode->FirstChild();
    if(! resnode) {
        std::stringstream p_tmp;
        p_tmp << "op of test " << testCount
              << " of case " << caseCount
              << " has no expected result child";
        throw(runtime_error(p_tmp.str()));
    }
    op.setExpected(resnode->Value());

    util::Profile profile("op");

    // TODO overlayareatest
    // TODO unionlength
    // TODO unionarea
    // TODO areatest




    auto functions = getFunctions();
    Result r;

    try {
        r = functions[op.getName()](op.getArgs());
    } catch(const std::exception& e) {
        if (expected_result == "exception") {
            success = true;
            actual_result = "exception";
        }
<<<<<<< Updated upstream

        else if(opName == "polygonize") {

            GeomPtr gRes(wktreader->read(opRes));
            gRes->normalize();

            Polygonizer plgnzr;
            plgnzr.add(gA);


            auto polys = plgnzr.getPolygons();
            GeomPtr gRealRes(factory->createGeometryCollection(std::move(polys)));
            gRealRes->normalize();


            if(gRes->compareTo(gRealRes.get()) == 0) {
                success = 1;
            }

            actual_result = printGeom(gRealRes.get());
            expected_result = printGeom(gRes.get());

            if(testValidOutput) {
                success &= int(testValid(gRealRes.get(), "result"));
            }
=======
        else {
            std::cerr << "EXCEPTION on case " << caseCount
                      << " test " << testCount << ": " << e.what()
                      << std::endl;
            actual_result = e.what();
>>>>>>> Stashed changes
        }
    }
    catch(...) {
        std::cerr << "Unknown EXEPTION on case "
                  << caseCount
                  << std::endl;
        actual_result = "Unknown exception thrown";
    }

    if(success) {
        ++succeeded;
    }
    else {
        ++failed;
    }

    if((!success && verbose) || verbose > 0) {
        printTest(!!success, expected_result, actual_result, profile);
    }

    if(test_predicates && gB && gA) {
        runPredicates(gA, gB);
    }


        else if(opName == "overlayareatest") {

            std::string maxDiffOp;
            std::stringstream p_tmp;
            double maxDiff = 1e-6;
            double areaDiff = areaDelta(gA, gB, maxDiffOp, maxDiff, p_tmp);

            // Debug output of actual geometries returned
            if (areaDiff < maxDiff && false) {
                std::cout << p_tmp.str();
            }

            p_tmp.str("");
            p_tmp << maxDiffOp << ": " << areaDiff;
            actual_result = p_tmp.str();
            p_tmp.str("");
            p_tmp << maxDiff;
            expected_result = p_tmp.str();

            if (areaDiff < maxDiff)
                success = 1;
        }

        else if(opName == "unionlength") {

            char* rest;
            GeomPtr result = OverlayNGRobust::Union(gA);
            double resultLength = result->getLength();
            double expectedLength = std::strtod(opRes.c_str(), &rest);
            if(rest == opRes.c_str()) {
                throw std::runtime_error("malformed testcase: missing expected length 'unionlength' op");
            }

            std::stringstream ss;
            ss << resultLength;
            actual_result = ss.str();

            if (std::abs(expectedLength-resultLength) / expectedLength < 1e-3) {
                success = 1;
            }

        }

        else if(opName == "unionarea") {

            char* rest;
            GeomPtr result = OverlayNGRobust::Union(gA);
            double resultArea  = result->getArea();
            double expectedArea = std::strtod(opRes.c_str(), &rest);
            if(rest == opRes.c_str()) {
                throw std::runtime_error("malformed testcase: missing expected area 'unionarea' op");
            }

            std::stringstream ss;
            ss << resultArea;
            actual_result = ss.str();

            if (std::abs(expectedArea-resultArea) / expectedArea < 1e-3) {
                success = 1;
            }

        }

        else if(opName == "areatest") {
            char* rest;
            double toleratedDiff = std::strtod(opRes.c_str(), &rest);
            int validOut = 1;

            if(rest == opRes.c_str()) {
                throw std::runtime_error("malformed testcase: missing tolerated area difference in 'areatest' op");
            }

            if(verbose > 1) {
                std::cerr << "Running intersection for areatest" << std::endl;
            }
            GeomPtr gI(gA->intersection(gB));

            if(testValidOutput) {
                validOut &= int(testValid(gI.get(), "areatest intersection"));
            }

            if(verbose > 1) {
                std::cerr << "Running difference(A,B) for areatest" << std::endl;
            }

            GeomPtr gDab(gA->difference(gB));

            if(testValidOutput) {
                validOut &= int(testValid(gI.get(), "areatest difference(a,b)"));
            }

            if(verbose > 1) {
                std::cerr << "Running difference(B,A) for areatest" << std::endl;
            }

            GeomPtr gDba(gB->difference(gA));

            if(testValidOutput) {
                validOut &= int(testValid(gI.get(), "areatest difference(b,a)"));
            }

            if(verbose > 1) {
                std::cerr << "Running symdifference for areatest" << std::endl;
            }

            GeomPtr gSD(gA->symDifference(gB));

            if(testValidOutput) {
                validOut &= int(testValid(gI.get(), "areatest symdifference"));
            }

            if(verbose > 1) {
                std::cerr << "Running union for areatest" << std::endl;
            }

            GeomPtr gU(gA->Union(gB));

            double areaA = gA->getArea();
            double areaB = gB->getArea();
            double areaI = gI->getArea();
            double areaDab = gDab->getArea();
            double areaDba = gDba->getArea();
            double areaSD = gSD->getArea();
            double areaU = gU->getArea();

            double maxdiff = 0;
            std::string maxdiffop;

            // @ : symdifference
            // - : difference
            // + : union
            // ^ : intersection

            // A == ( A ^ B ) + ( A - B )
            double diff = std::fabs(areaA - areaI - areaDab);
            if(diff > maxdiff) {
                maxdiffop = "A == ( A ^ B ) + ( A - B )";
                maxdiff = diff;
            }

            // B == ( A ^ B ) + ( B - A )
            diff = std::fabs(areaB - areaI - areaDba);
            if(diff > maxdiff) {
                maxdiffop = "B == ( A ^ B ) + ( B - A )";
                maxdiff = diff;
            }

            //  ( A @ B ) == ( A - B ) + ( B - A )
            diff = std::fabs(areaDab + areaDba - areaSD);
            if(diff > maxdiff) {
                maxdiffop = "( A @ B ) == ( A - B ) + ( B - A )";
                maxdiff = diff;
            }

            //  ( A u B ) == ( A ^ B ) + ( A @ B )
            diff = std::fabs(areaI + areaSD - areaU);
            if(diff > maxdiff) {
                maxdiffop = "( A u B ) == ( A ^ B ) + ( A @ B )";
                maxdiff = diff;
            }

            if(maxdiff <= toleratedDiff) {
                success = 1 && validOut;
            }

            std::stringstream p_tmp;
            p_tmp << maxdiffop << ": " << maxdiff;
            actual_result = p_tmp.str();
            expected_result = opRes;

        }


        else {
            std::cerr << *curr_file << ":";
            std::cerr << " case" << caseCount << ":";
            std::cerr << " test" << testCount << ": "
                      << opName << "(" << opSig << ")";
            std::cerr << ": skipped (unrecognized)." << std::endl;
            return;
        }

    }
    catch(const std::exception& e) {
        if (expected_result == "exception") {
            success = true;
            actual_result = "exception";
        }
        else {
            std::cerr << "EXCEPTION on case " << caseCount
                      << " test " << testCount << ": " << e.what()
                      << std::endl;
            actual_result = e.what();
        }
    }
    catch(...) {
        std::cerr << "Unknown EXEPTION on case "
                  << caseCount
                  << std::endl;
        actual_result = "Unknown exception thrown";
    }

    if(success) {
        ++succeeded;
    }
    else {
        ++failed;
    }

    if((!success && verbose) || verbose > 0) {
        printTest(!!success, expected_result, actual_result, profile);
    }

    if(test_predicates && gB && gA) {
        runPredicates(gA, gB);
    }

}

void
XMLTester::runPredicates(const geom::Geometry* p_gA, const geom::Geometry* p_gB)
{
    std::cout << "\t    Equals:\tAB=" << (p_gA->equals(p_gB) ? "T" : "F") << ", BA=" << (p_gB->equals(
                  p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t  Disjoint:\tAB=" << (p_gA->disjoint(p_gB) ? "T" : "F") << ", BA=" << (p_gB->disjoint(
                  p_gA) ? "T" : "F") << std::endl;
    std::cout << "\tIntersects:\tAB=" << (p_gA->intersects(p_gB) ? "T" : "F") << ", BA=" << (p_gB->intersects(
                  p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t   Touches:\tAB=" << (p_gA->touches(p_gB) ? "T" : "F") << ", BA=" << (p_gB->touches(
                  p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t   Crosses:\tAB=" << (p_gA->crosses(p_gB) ? "T" : "F") << ", BA=" << (p_gB->crosses(
                  p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t    Within:\tAB=" << (p_gA->within(p_gB) ? "T" : "F") << ", BA=" << (p_gB->within(
                  p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t  Contains:\tAB=" << (p_gA->contains(p_gB) ? "T" : "F") << ", BA=" << (p_gB->contains(
                  p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t  Overlaps:\tAB=" << (p_gA->overlaps(p_gB) ? "T" : "F") << ", BA=" << (p_gB->overlaps(
                  p_gA) ? "T" : "F") << std::endl;

    std::cout << "\t  Prepared Disjoint:\tAB=" << (prepare(p_gA)->disjoint(p_gB) ? "T" : "F") << ", BA=" << (prepare(
                  p_gB)->disjoint(p_gA) ? "T" : "F") << std::endl;
    std::cout << "\tPrepared Intersects:\tAB=" << (prepare(p_gA)->intersects(p_gB) ? "T" : "F") << ", BA=" << (prepare(
                  p_gB)->intersects(p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t   Prepared Touches:\tAB=" << (prepare(p_gA)->touches(p_gB) ? "T" : "F") << ", BA=" << (prepare(
                  p_gB)->touches(p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t   Prepared Crosses:\tAB=" << (prepare(p_gA)->crosses(p_gB) ? "T" : "F") << ", BA=" << (prepare(
                  p_gB)->crosses(p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t    Prepared Within:\tAB=" << (prepare(p_gA)->within(p_gB) ? "T" : "F") << ", BA=" << (prepare(
                  p_gB)->within(p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t  Prepared Contains:\tAB=" << (prepare(p_gA)->contains(p_gB) ? "T" : "F") << ", BA=" << (prepare(
                  p_gB)->contains(p_gA) ? "T" : "F") << std::endl;
    std::cout << "\t Prepared Overlaps:\tAB=" << (prepare(p_gA)->overlaps(p_gB) ? "T" : "F") << ", BA=" << (prepare(
                  p_gB)->overlaps(p_gA) ? "T" : "F") << std::endl;
}

XMLTester::~XMLTester()
{
}


static void
usage(char* me, int exitcode, std::ostream& os)
{
    os << "Usage: " << me << " [options] <test> [<test> ...]" << std::endl;
    os << "Options: " << std::endl;
    os << " -v                  Verbose mode "
       << "(multiple -v increment verbosity)" << std::endl
       << "--test-valid-output  Test output validity" << std::endl
       << "--test-valid-input   Test input validity" << std::endl
       << "--sql-output         Produce SQL output" << std::endl
       << "--wkb-output         Print Geometries as HEXWKB" << std::endl;

    std::exit(exitcode);
}

void
request_interrupt(int)
{
    geos::util::Interrupt::request();
}

int
main(int argC, char* argV[])
{
    int verbose = 0;
    bool sql_output = false;

#if defined(_MSC_VER) && defined(GEOS_TEST_USE_STACKWALKER)
    InitAllocCheck();
    {
#endif

        if(argC < 2) {
            usage(argV[0], 1, std::cerr);
        }

        signal(15, request_interrupt);

        XMLTester tester;
        tester.setVerbosityLevel(verbose);

        for(int i = 1; i < argC; ++i) {
            // increment verbosity level
            if(! std::strcmp(argV[i], "-v")) {
                ++verbose;
                tester.setVerbosityLevel(verbose);
                continue;
            }
            if(! std::strcmp(argV[i], "--test-valid-output")) {
                tester.testOutputValidity(true);
                continue;
            }
            if(! std::strcmp(argV[i], "--sql-output")) {
                sql_output = true;
                tester.setSQLOutput(sql_output);
                continue;
            }
            if(! std::strcmp(argV[i], "--wkb-output")) {
                sql_output = true;
                tester.setHEXWKBOutput(sql_output);
                continue;
            }
            if(! std::strcmp(argV[i], "--test-valid-input")) {
                tester.testInputValidity(true);
                continue;
            }

            std::string source = argV[i];

            try {
                tester.run(source);
            }
            catch(const std::exception& exc) {
                std::cerr << exc.what() << std::endl;
            }
        }

        if(! sql_output) {
            tester.resultSummary(std::cout);
        }
        else {
            tester.resultSummary(std::cerr);
        }

        io::Unload::Release();

        return tester.getFailuresCount();

#if defined(_MSC_VER) && defined(GEOS_TEST_USE_STACKWALKER)
    }
    DeInitAllocCheck();
#endif

}
