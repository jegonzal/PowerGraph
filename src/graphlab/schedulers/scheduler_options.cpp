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

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <graphlab/schedulers/scheduler_options.hpp>

namespace graphlab {
  
  std::string 
  scheduler_options::parse_scheduler_string(std::string scheduler_raw) {
    std::pair<std::string, scheduler_options> ret;
    // Break the scheduler string appart
    size_t first_paren = scheduler_raw.find_first_of('(');
    size_t last_paren = scheduler_raw.find_last_of(')');
    std::string scheduler = scheduler_raw.substr(0, first_paren);
    std::string arguments;
    // Fill in the arguments if such are possibe
    if(first_paren != std::string::npos &&
       last_paren != std::string::npos) {
      arguments = scheduler_raw.substr(first_paren + 1,
                                       last_paren - first_paren - 1 );
    }
    if(!arguments.empty()) {
      std::replace(arguments.begin(), arguments.end(), ',', ' ');
      std::replace(arguments.begin(), arguments.end(), ';', ' ');        
      std::stringstream arg_strm(arguments);
      parse_options(arg_strm);
    }     
    return scheduler;
  }



std::ostream& operator<<(std::ostream& out,
                         const graphlab::scheduler_options& opts) {
  // save the format flags
  std::ios_base::fmtflags fmt = out.flags();
  
  std::map<std::string,
    graphlab::scheduler_options::scheduler_option_values>::const_iterator
    i = opts.options.begin();
  while(i != opts.options.end()) {
    
    //out.setf(std::ios::left);
    out << std::setw(18) << std::left << i->first;
    
    out << std::setw(2) << "= ";
    
    //out.setf(std::ios::right);
    i->second.anyval.print(out);
    out << std::endl;
    ++i;
  }
  // reset the format flags
  out.flags(fmt);
  out << std::endl;
  return out;
}

}
