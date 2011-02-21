
#ifndef GRAPHLAB_IMETRICS_REPORTER
#define GRAPHLAB_IMETRICS_REPORTER

#include <map>

#include <graphlab/metrics/metrics.hpp>

namespace graphlab {
   class imetrics_reporter {
   
        public:
            virtual void do_report(std::string name, std::string id, std::map<std::string, metrics_entry> &  entries) = 0;
    
   };

};

#endif


