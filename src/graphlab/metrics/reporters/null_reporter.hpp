
#ifndef GRAPHLAB_NULL_REPORTER
#define GRAPHLAB_NULL_REPORTER

#include <graphlab/metrics/imetrics_reporter.hpp>

/**
 * Simple metrics reporter that dumps metrics to
 * standard output.
 */

namespace graphlab {

  class null_reporter : public imetrics_reporter {

  public:
    virtual void do_report(std::string name, std::string ident, std::map<std::string, metrics_entry> & entries) {
    }

  };

}



#endif

