#pragma once

#include <geos/geom/Geometry.h>

#include <memory>
#include <sstream>

namespace geos {
namespace xmltester {

class TypeError : public std::exception {
public:
    TypeError() {}
};

class Value {
public:
    Value() : m_type(UNDEFINED) {}

    Value(bool value) : m_type(BOOL), m_bool(value) {}
    Value(int value) : m_type(INT), m_int(value) {}
    Value(double value) : m_type(DOUBLE), m_double(value) {}
    Value(const std::string& value) : m_type(STRING), m_string(value) {}
    Value(const char* value) : m_type(STRING), m_string(value) {}
    Value(std::unique_ptr<geom::Geometry> value) : m_type(GEOMETRY), m_geometry(std::move(value)) {}

    enum Type {
        BOOL,
        DOUBLE,
        INT,
        STRING,
        GEOMETRY,
        UNDEFINED
    };

    std::string toString() const {
        std::stringstream ss;

        switch(getType()) {
            case BOOL: ss << (getBool() ? "true" : "false"); break;
            case DOUBLE: ss << getDouble(); break;
            case INT: ss << getInt(); break;
            case STRING: ss << getString(); break;
            case GEOMETRY: ss << getGeometry()->toString(); break;
            case UNDEFINED: break;
        }

        return ss.str();
    }

    Type getType() const {
        return m_type;
    }

    bool isGeometry() const {
        return getType() == GEOMETRY;
    }

    bool isNull() const {
        return getType() == UNDEFINED;
    }

    bool getBool() const {
        if (m_type == BOOL) {
            return m_bool;
        }

        if (m_type == STRING) {
            std::string ustring = m_string;
            for (auto& c : ustring) {
                c = static_cast<char>(toupper(c));
            }

            if (ustring == "TRUE") {
                return true;
            }
            if (ustring == "FALSE") {
                return false;
            }
        }

        throw TypeError();
    }

    double getDouble() const {
        if (m_type == DOUBLE) {
            return m_double;
        }

        if (m_type == INT) {
            return static_cast<double>(m_int);
        }

        if (m_type == STRING) {
            std::size_t pos;
            double d = std::stod(m_string, &pos);
            if (pos == m_string.size()) {
                return d;
            }
        }

        throw TypeError();
    }

    int getInt() const {
        if (m_type == INT) {
            return m_int;
        }

        if (m_type == STRING) {
            std::size_t pos;
            int i = std::stoi(m_string, &pos);
            if (pos == m_string.size()) {
                return i;
            }
        }

        throw TypeError();
    }

    double getIntOr(int default_value) const {
        if (getType() == UNDEFINED) {
            return default_value;
        }
        return getInt();
    }

    double getDoubleOr(double default_value) const {
        if (getType() == UNDEFINED) {
            return default_value;
        }
        return getDouble();
    }

    const std::string& getString() const {
        checkType(STRING);
        return m_string;
    }

    const geom::Geometry* getGeometry() const {
        checkType(GEOMETRY);
        return m_geometry.get();
    }

    bool operator!=(const Value& other) const {
        return !(*this == other);
    }

    bool operator==(const Value& other) const {
        try {
            switch (getType()) {
                case BOOL: return getBool() == other.getBool();
                case DOUBLE: return getDouble() == other.getDouble();
                case INT: return getInt() == other.getInt();
                case STRING: return getString() == other.getString();
                case GEOMETRY: return getGeometry()->equals(other.getGeometry());
                default: return other.getType() == UNDEFINED;
            }
        } catch (const TypeError&) {
            return false;
        }
    }

private:
    void typeError() const {
        throw std::runtime_error("Incorrect type access");
    }

    void checkType(Type expected) const {
        if (m_type == STRING &&
                ( expected == INT ||
                  expected == DOUBLE )) {
                return;
        }

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

}
}
