#include <boost/version.hpp>
#include <boost/filesystem.hpp>

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>




#include <graphlab/util/fs_util.hpp>



void graphlab::fs_util::
list_files_with_suffix(const std::string& pathname,
                       const std::string& suffix,
                       std::vector<std::string>& files) {
  namespace fs = boost::filesystem;
  fs::path path(pathname);
  assert(fs::exists(path));
  for(fs::directory_iterator iter( path ), end_iter; 
      iter != end_iter; ++iter) {
    if( ! fs::is_directory(iter->status()) ) {

#if BOOST_VERSION >=104601
      std::string filename(iter->path().filename().string());
#else
      std::string filename(iter->path().filename());
#endif
      size_t pos = 
        filename.size() >= suffix.size()?
        filename.size() - suffix.size() : 0;
      std::string ending(filename.substr(pos));
      if(ending == suffix) {
#if BOOST_VERSION >=104601
        files.push_back(iter->path().filename().string());
#else
        files.push_back(iter->path().filename());
#endif
      }
    }
  }
  std::sort(files.begin(), files.end());
} // end of list files with suffix  





std::string graphlab::fs_util::
change_suffix(const std::string& fname,
              const std::string& new_suffix) {             
  size_t pos = fname.rfind('.');
  assert(pos != std::string::npos); 
  const std::string new_base(fname.substr(0, pos));
  return new_base + new_suffix;
} // end of change_suffix


