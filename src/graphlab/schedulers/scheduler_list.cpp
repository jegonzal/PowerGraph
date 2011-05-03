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

#include <algorithm>
#include <graphlab/schedulers/scheduler_list.hpp>
#include <graphlab/util/stl_util.hpp>


namespace graphlab {
  
  std::vector<std::string> get_scheduler_names() {
    std::vector<std::string> ret;
#define __APPEND_TO_RET__(r_unused, data_unused, i,  elem)      \
    ret.push_back(BOOST_PP_TUPLE_ELEM(3,0,elem));
    BOOST_PP_SEQ_FOR_EACH_I(__APPEND_TO_RET__, _, __SCHEDULER_LIST__)
#undef __APPEND_TO_RET__
      return ret;
  }


  std::string get_scheduler_names_str() {
    std::string ret;
    std::vector<std::string> schednames;
    schednames = get_scheduler_names();
    for (size_t i = 0; i < schednames.size(); ++i) {
      if (i > 0) {
        ret = ret + ", ";
      }
      ret = ret + schednames[i];
    }
    return ret;
  }

  std::string add_line_breaks(const std::string &s, size_t numcols) {
    size_t pos = 0;
    std::string ret;
    while(pos < s.length() - 1) {
      size_t oldpos = pos;
      pos = std::min(pos + numcols, s.length());
    
      size_t newpos = pos;
      // search backward for a space if we are not at the end of the
      // string
      if (pos < s.length()) {
        newpos = s.rfind(" ", pos);
      }
      // if we get back to the old position, or we fail to find a
      // space, force the break
      if (newpos  == std::string::npos || newpos == oldpos) {
        newpos = pos;
      }
      // break
      ret = ret + trim(s.substr(oldpos, newpos - oldpos)) + "\n";
      pos = newpos;
    }
    return ret;
  }

  void print_scheduler_info(std::string s, std::ostream &out) {
    // this is annoying... I need to instantiate the graph<char, char> type to
    // even call the scheduler
#define __GENERATE_SCHEDULER_HELP__(r_unused, data_unused, i,  elem)    \
    BOOST_PP_EXPR_IF(i, else) if (s == BOOST_PP_TUPLE_ELEM(3,0,elem)) { \
      out << "\n";                                                      \
      out << BOOST_PP_TUPLE_ELEM(3,0,elem) << " scheduler\n";           \
      out << std::string(50, '-') << std::endl;                         \
      out << add_line_breaks(BOOST_PP_TUPLE_ELEM(3,2,elem), 50) << "\n"; \
      out << "Options: \n";                                             \
      BOOST_PP_TUPLE_ELEM(3,1,elem)<graph<char,char> >::print_options_help(out); \
    }
    /*
     * if (scheduler == "sweep") {
     *   sweep_scheduler<graph<char,char> >::print_options_help(out);
     * }
     * ...
     */
    // generate the construction calls
    BOOST_PP_SEQ_FOR_EACH_I(__GENERATE_SCHEDULER_HELP__, _, __SCHEDULER_LIST__)
    else {
      out << "Scheduler " << s << " not found" << "\n";
    }
#undef __GENERATE_SCHEDULER_HELP__
  }
}



