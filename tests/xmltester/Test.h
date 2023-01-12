#pragma once

#include "Args.h"
#include "Value.h"
#include "XMLTesterUtil.h"

namespace geos {
namespace xmltester {
class Test {

public:
    const std::string& getName() const {
        return m_name;
    }

    void setName(const std::string& value) {
        m_name = XMLTesterUtil::trimBlanks(value);
        XMLTesterUtil::tolower(m_name);
    }

    void setExpected(const char* value) {
        setExpected(std::string(value));
    }

    void setExpected(const std::string& value) {
        m_expected = Value(XMLTesterUtil::trimBlanks(value));
    }

    void setExpected(Value value) {
        m_expected = std::move(value);
    }

    const Value& getExpected() const {
        return m_expected;
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
            if (!opArg.isNull()) {
                if (opSig == "") {
                    opSig += ", ";
                }
                opSig += opArg.toString();
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

    void setGeomA(std::shared_ptr<const geom::Geometry> g) {
        m_args.setGeomA(g);
    }

    void setGeomB(std::shared_ptr<const geom::Geometry> g) {
        m_args.setGeomB(g);
    }

    void setArg(std::size_t i, const char* s) {
        if (s)
            setArg(i, std::string(s));
    }

    void setArg(std::size_t i, const std::string& s) {
        if (s != "")
            m_args.set(i, s);
    }

    void setUsePrepared(bool value) {
        m_args.setUsePrepared(value);
    }

private:
    std::string m_name;

    Value m_expected;
    Value m_actual;

    Args m_args;
};

}
}
