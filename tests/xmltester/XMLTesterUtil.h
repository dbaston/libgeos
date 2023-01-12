#pragma once

#include <string>

namespace geos {
namespace xmltester {

class XMLTesterUtil {
public:
    static std::string trimBlanks(const std::string& in);
    static void tolower(std::string& str);
};

}
}
