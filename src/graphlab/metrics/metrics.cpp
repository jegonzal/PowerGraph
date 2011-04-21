#include <graphlab/metrics/metrics.hpp>

namespace graphlab {
  
void metrics::report(imetrics_reporter & reporter) {
  if (name != "") {
    reporter.do_report(name, ident, entries);
  }
}

}
