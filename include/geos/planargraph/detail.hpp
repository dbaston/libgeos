#include <vector>
#include <algorithm>

namespace {

template <typename T>
void
find_and_erase(T *what, std::vector<T*> &where) {
	if (!what) return;
	where.erase(std::remove(where.begin(), where.end(), what), where.end());
}


/** \brief
 * Sets the Marked state for the values of each map
 * container element, from start to end iterator.
 *
 * @param start the start element
 * @param end one past the last element
 * @param marked the state to set the visited flag to
 */
template <typename T>
void setMarkedMap(T &p_map, bool marked) {
    for(auto &e : p_map) {
        e.second->setMarked(marked);
    }
}

/** \brief
 * Sets the Visited state for the values of each map
 * container element, from start to end iterator.
 *
 * @param start the start element
 * @param end one past the last element
 * @param visited the state to set the visited flag to
 */
template <typename T>
void setVisitedMap(T &p_map, bool visited) {
    for(auto &e : p_map) {
        e.second->setVisited(visited);
    }
}

}
