#pragma once

#include <geos/geom/Geometry.h>
#include <geos/geom/prep/PreparedGeometry.h>
#include <geos/geom/prep/PreparedGeometryFactory.h>
#include "Value.h"

#include <array>
#include <memory>

namespace geos {
namespace xmltester {

class Args {

public:

    void set(std::size_t index, const std::string& value) {
        m_args[index - 1] = Value(value);
    }

    void setGeomA(std::shared_ptr<const geom::Geometry> g) {
        m_geoms[0] = g;
    }

    void setGeomB(std::shared_ptr<const geom::Geometry> g) {
        m_geoms[1] = g;
    }

    bool has(std::size_t index) const {
        const auto& val = get(index);
        return !val.isNull();
    }

    const geom::Geometry* A() const {
        return swapGeomArgs() ? m_geoms[1].get() : m_geoms[0].get();
    }

    const geom::Geometry* B() const {
        return swapGeomArgs() ? m_geoms[0].get() : m_geoms[1].get();
    }

    const geom::prep::PreparedGeometry* pA() const {
        prepareGeometries();
        return swapGeomArgs() ? m_prep[1].get() : m_prep[0].get();
    }

    const geom::prep::PreparedGeometry* pB() const {
        prepareGeometries();
        return swapGeomArgs() ? m_prep[0].get() : m_prep[1].get();
    }

    bool usePrepared() const {
        return m_useprepared;
    }

    void setUsePrepared(bool value) {
        m_useprepared = value;
    }

    const Value& get(std::size_t index) const {
        return m_args[index - 1];
    }

    const Value& operator[](std::size_t index) const {
        return get(index);
    }

private:
    bool swapGeomArgs() const {
        return (get(1) == "B" || get(1) == "b") && m_geoms[1] != nullptr;
    }

    void prepareGeometries() const {
        for (std::size_t i = 0; i < m_prep.size(); i++) {
            if (m_prep[i] == nullptr) {
                m_prep[i] = geom::prep::PreparedGeometryFactory::prepare(m_geoms[i].get());
            }
        }
    }

    static constexpr std::size_t MAX_ARGS = 4;
    std::array<Value, MAX_ARGS> m_args;
    std::array<std::shared_ptr<const geom::Geometry>, 2> m_geoms;
    mutable std::array<std::unique_ptr<geom::prep::PreparedGeometry>, 2> m_prep;
    bool m_useprepared;

    using const_iterator = decltype(m_args.cbegin());
public:
    const_iterator begin() const {
        return m_args.begin();
    }

    const_iterator end() const {
        return m_args.end();
    }
};

}
}

