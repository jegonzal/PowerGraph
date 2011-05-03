/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

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
