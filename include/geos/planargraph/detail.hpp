#include <vector>
#include <algorithm>

namespace {

template <typename T>
void
find_and_erase(T *what, std::vector<T*> &where) {
	if (!what) return;
	where.erase(std::remove(where.begin(), where.end(), what), where.end());
}

}
