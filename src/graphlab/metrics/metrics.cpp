
#include <graphlab/metrics/metrics.hpp>
#include <graphlab/metrics/imetrics_reporter.hpp>


using namespace graphlab;

static std::map<std::string, metrics *> metricsmap; 

/* USE THIS TO CREATE A NEW METRICS INSTANCE */

/**
  * Creates a new named metrics instance, if it does not exist yet.
  * If it exist, and flag create_always_new is FALSE, previously created
  * instance is returned. Otherwise, a new instance with suffix ":n", where
  * n is the number of such named instances is created.
  */
metrics & metrics::create_metrics_instance(std::string name, bool create_always_new = false) {
    static mutex mlock;
    metrics * m;
    mlock.lock();
    if (metricsmap.count(name) == 0) {
        m = new metrics(name, name);
        metricsmap[name] = m;
    } else if (create_always_new) {
        int n = 1;
        // Count how many metrics instances with same name are created
        std::map<std::string, metrics *>::iterator miter;
        for(miter = metricsmap.begin(); miter != metricsmap.end(); ++miter) {
            if (miter->second->ident == name) n++;
        }
        char cc[255];
        sprintf(cc, "%s:%d", name.c_str(), n);
        std::string ident = std::string(cc);
        m = new metrics(name, ident);
        metricsmap[ident] = m;
    } else {
        m = metricsmap[name];
    }
    mlock.unlock();
    return *m;
}

metrics & metrics::create_metrics_instance(std::string name) {
    return metrics::create_metrics_instance(name, false);
}

/**
  * Iterates all metrics instances and asks
  * them to report.
  */
void metrics::report_all(imetrics_reporter & reporter) {
   std::map<std::string, metrics *>::iterator miter;
   for(miter = metricsmap.begin(); miter != metricsmap.end(); ++miter) {
      miter->second->report(reporter);
   }
}

void metrics::report(imetrics_reporter & reporter) {
    reporter.do_report(this->name, this->ident, entries);
}
        