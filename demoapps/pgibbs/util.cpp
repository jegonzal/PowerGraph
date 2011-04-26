#include <ctime>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>


#include "util.hpp"



size_t file_line_count(const std::string& experiment_file) {
  std::ifstream fin(experiment_file.c_str());
  size_t lines = 0;
  std::string line;
  while(getline(fin, line)) lines++;
  fin.close();
  return lines;
}


std::string make_filename(const std::string& base,
                          const std::string& suffix,
                          const size_t number) {
  std::stringstream strm;
  strm << base
       << std::setw(10) << std::setfill('0')
       << number
       << suffix;
  std::cout << strm.str() << std::endl;
  return strm.str();
}
