#include "XMLTesterUtil.h"

#include <algorithm>

namespace geos {
namespace xmltester {

void
XMLTesterUtil::tolower(std::string& str)
{
    std::transform(
        str.begin(),
        str.end(),
        str.begin(),
        [](char c){ return (char)std::tolower(c); }
        );
}

std::string
XMLTesterUtil::trimBlanks(const std::string& in)
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

}
}
