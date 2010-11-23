#ifndef SCHEDULER_OPTION_CACHE
#define SCHEDULER_OPTION_CACHE
#include <sstream>
#include <graphlab/schedulers/ischeduler.hpp>
class scheduler_option_cache {
 private:
   std::stringstream cmdoptions;
   std::vector<std::pair<graphlab::scheduler_options::options_enum, void*> > options;
 public:
  inline void parse_options(std::stringstream &strm) {
    cmdoptions.str(strm.str());
  }
  
  inline void set_option(graphlab::scheduler_options::options_enum opt, void* value) {
    options.push_back(std::make_pair(opt,value));
  }

  template <typename Scheduler> 
  void apply_options(Scheduler &sched) {
    if (cmdoptions.str().length() > 0) {
      sched.parse_options(cmdoptions);
    }
    for (size_t i = 0; i < options.size(); ++i) {
      sched.set_option(options[i].first, options[i].second);
    }
  }

  scheduler_option_cache& operator=(const scheduler_option_cache &cache) {
    std::string cachestr = cache.cmdoptions.str();
    if (cachestr.length() != 0) {
      cmdoptions.str(cachestr);
    }

    if (cache.options.size() > 0) {
      options = cache.options;
    }
    return *this;
  }
};
#endif