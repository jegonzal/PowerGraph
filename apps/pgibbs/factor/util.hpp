#ifndef PGIBBS_UTIL_HPP
#define PGIBBS_UTIL_HPP

#include <string>
#include <sstream>
#include <iostream>



std::string make_filename(const std::string& base,
                          const std::string& suffix,
                          size_t number) {
  std::stringstream strm;
  strm << base
       << std::setw(10) << std::setfill('0')
       << number
       << suffix;
  std::cout << strm.str() << std::endl;
  return strm.str();
}





#endif
