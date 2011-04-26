#ifndef GRAPHLAB_FS_UTIL
#define GRAPHLAB_FS_UTIL

#include <string>
#include <vector>


namespace graphlab {

  namespace fs_util {

    /**
     * List all the files with the given suffix at the pathname
     * location
     * \ingroup util_internal
     */
    void list_files_with_suffix(const std::string& pathname,
                                const std::string& suffix,
                                std::vector<std::string>& files);

    /// \ingroup util_internal
    std::string change_suffix(const std::string& fname,
                                     const std::string& new_suffix);

  }; // end of fs_utils


}; // end of graphlab
#endif
