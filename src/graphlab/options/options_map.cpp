/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <graphlab/options/options_map.hpp>

namespace graphlab {
  
  std::string 
  options_map::parse_string(std::string options_raw) {
    std::pair<std::string, options_map> ret;
    // Break the string appart
    size_t first_paren = options_raw.find_first_of('(');
    size_t last_paren = options_raw.find_last_of(')');
    std::string fun_name = options_raw.substr(0, first_paren);
    std::string arguments;
    // Fill in the arguments if such are possibe
    if(first_paren != std::string::npos &&
       last_paren != std::string::npos) {
      arguments = options_raw.substr(first_paren + 1,
                                     last_paren - first_paren - 1 );
    }
    if(!arguments.empty()) {
      std::replace(arguments.begin(), arguments.end(), ',', ' ');
      std::replace(arguments.begin(), arguments.end(), ';', ' ');        
      std::stringstream arg_strm(arguments);
      parse_options(arg_strm);
    }     
    return fun_name;
  }



  std::ostream& operator<<(std::ostream& out,
                           const graphlab::options_map& opts) {
    // save the format flags
    std::ios_base::fmtflags fmt = out.flags();
  
    std::map<std::string,
             graphlab::options_map::option_values>::const_iterator
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

}; // end of namespace graphlab

