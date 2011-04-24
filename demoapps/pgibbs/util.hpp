#ifndef PGIBBS_UTIL_HPP
#define PGIBBS_UTIL_HPP

#include <string>

//! Get the number of lines in the file
size_t file_line_count(const std::string& experiment_file);


//! make a filename from base sufix and number
std::string make_filename(const std::string& base,
                          const std::string& suffix,
                          const size_t number);




#endif
