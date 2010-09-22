/**
 * Simple metrics consodilation tool.
 * written by Aapo Kyrola
 **/
#ifndef DISTR_METRICS_HPP
#define DISTR_METRICS_HPP

#include <queue>
#include <cmath>
#include <cassert>
#include <vector>
#include <climits>


#include <logger/logger.hpp> 
#include <graphlab/parallel/pthread_tools.hpp> 
#include <graphlab/distributed/distributed_control.hpp> 


#include <graphlab/macros_def.hpp>

namespace graphlab {

  
class distributed_metrics {

  private:
 
 
 distributed_metrics(distributed_control * _dc) :
            stats(_dc->numprocs()) {
    dc = _dc;
 }
 
 public:
  

 
 static distributed_metrics * instance(distributed_control * dc);
 
 // Remote call
 static void remote_set_value(distributed_control& _dc, size_t source, void *ptr, size_t len,
        std::string key, double value);
        
 void set_value(std::string key, double value) {
    set_value(dc->procid(), key, value);
 }
 
 void set_value(size_t procid, std::string key, double value) {
    if (dc->procid() == 0) {
        synclock.lock();
        stats[procid][key] = value;
        synclock.unlock();
    } else {
        dc->remote_callxs(0, distributed_metrics::remote_set_value,
            NULL, 0, key, value);
    }
 }
 
 double get_value(int procid, std::string key) {
    std::map<std::string, double>::iterator iter = stats[procid].find(key);
    if (iter == stats[procid].end()) {
        ASSERT_MSG("A key %s was not found for procid %d!", key.c_str(), procid);
        return -1; // Not reached, prevents warning.
    } else return iter->second;
 }
 
 void report() {
    ASSERT_EQ(dc->procid(), 0);
    
    char fname[255];
    sprintf(fname, "stats_%d.txt", dc->numprocs());
    
    FILE * statsf = fopen(fname,"w");
    FILE * benchf = fopen(".runstats.R", "w");
    int np = stats.size();
    
    printf("=============== REPORT (see also stats.txt) ============== \n");
    printf("key,");
    fprintf(statsf, "key,");
    for(int i=0; i<np; i++) {
        printf("%d,", i);
        fprintf(statsf, "%d,", i);
    }
    printf("avg,sum\n");
    fprintf(statsf, "avg,sum\n");

    std::map<std::string, double>::iterator iter = stats[0].begin();
    for(; iter != stats[0].end(); ++iter) {
        std::string key = iter->first;
        printf("%s,", key.c_str());
        fprintf(statsf, "%s,", key.c_str());
        double sum = 0.0;
        
        for(int p=0; p<np; p++) {
            double v = get_value(p, key);
            sum += v;
            printf("%lf,", v);
            fprintf(statsf, "%lf,", v);
        }
        
        // Average and sum
        printf("%lf,%lf\n",  sum/np, sum);
        fprintf(statsf, "%lf,%lf\n",  sum/np, sum);

        // Export for benchmarking system
        if (key == "taskcount" || key == "ncpus") {
            fprintf(benchf, "%s=%ld\n", key.c_str(), (long int) sum);
        } else if (key == "residual") {
 	        fprintf(benchf, "%s=%lf\n", key.c_str(),  sum);

        } else {
            fprintf(benchf, "%s=%lf\n", key.c_str(), sum/np);
        }    
    }
    printf("=============== END REPORT ============== \n");
    
    fclose(statsf);
    fclose(benchf);
 }
     
     
 private:
    distributed_control * dc;
    mutex synclock;
    std::vector< std::map<std::string, double> > stats;
     
};
    

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
