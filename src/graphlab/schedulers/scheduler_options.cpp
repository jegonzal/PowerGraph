#include <string>
#include <iostream>
#include <graphlab/schedulers/scheduler_options.hpp>

namespace graphlab {
  
std::ostream& operator<<(std::ostream& out,
                                const graphlab::scheduler_options& opts) {
  // save the format flags
  std::ios_base::fmtflags fmt = out.flags();

  std::map<std::string,
          graphlab::scheduler_options::scheduler_option_values>::const_iterator
                            i = opts.options.begin();
  while(i != opts.options.end()) {
    out.width(16);
    out.setf(std::ios::left);
    out << i->first;
    
    out.width(1);
    out << "=";
    
    out.width(8);
    out.setf(std::ios::right);
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