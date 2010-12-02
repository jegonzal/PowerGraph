#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
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


std::pair<std::string, scheduler_options> parse_scheduler_string(std::string scheduler_raw) {
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
  scheduler_options sched_opts;  
  if(!arguments.empty()) {
    std::replace(arguments.begin(), arguments.end(), ',', ' ');
    std::replace(arguments.begin(), arguments.end(), ';', ' ');        
    std::stringstream arg_strm(arguments);
    sched_opts.parse_options(arg_strm);
  }     
    
  ret.first = arguments;
  ret.second = sched_opts;
  return ret;

}
}