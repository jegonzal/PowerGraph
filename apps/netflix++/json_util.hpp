#include <boost/format.hpp>
#include <boost/algorithm/string/find_format.hpp>
#include <boost/algorithm/string/join.hpp>
template<typename Map>
std::string map2json(const Map& map) {
  typedef typename Map::value_type kv_type;
  typedef typename Map::const_iterator iter_type;
  std::stringstream ss;
  ss << "{\n";
  iter_type iter = map.begin();
  if (iter != map.end()) {
    ss << "\"" << (iter->first) << "\": "
       << "\"" << (iter->second) << "\"";
    ++iter;
    while (iter != map.end()) {
      ss << ",\n"
         << "\"" << (iter->first) << "\": "
         << "\"" << (iter->second) << "\"";
      ++iter;
    }
  }
  ss << "}\n";
  return ss.str();
}


template<typename Array>
std::string arr2json(const Array& arr) {
  std::stringstream ss;
  ss << "[\n";
  for (size_t i = 0; i < arr.size(); ++i) {
    if (i) 
      ss << ", ";
    ss << "\"" << (arr[i]) << "\"";
  }
  ss << "]\n";
  return ss.str();
}
