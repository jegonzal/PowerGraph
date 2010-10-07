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
#include <map>


#include <graphlab/logger/logger.hpp> 
#include <graphlab/parallel/pthread_tools.hpp> 
#include <graphlab/distributed/distributed_control.hpp> 
#include <graphlab/util/timer.hpp> 


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
 
 // Remote calls
 static void remote_set_value(distributed_control& _dc, size_t source, void *ptr, size_t len,
        std::string key, double value);
        
 static void remote_start_bandwidth_test(distributed_control& _dc, size_t source, void *ptr, size_t len);
 
 static void receive_bandwidth_test_payload(distributed_control& _dc, size_t source, void *ptr, size_t len);
 
 //// BANDWITH TEST
 
 /**
    * Starts a bandwidth test which sends 30 megabyte chunk to all
    * other nodes. Each node computes the time it took to send the data.
    */
    
 std::string bandwidth_key(int from, int to) {
    std::ostringstream os;
    os << "bandwidth_" << from;
    std::string buffer(os.str());
    return buffer;
 }
    
 #define BANDWIDTH_TEST_SIZE (30*1024*256)
 
 mutex condlock;
 conditional cond;
 int received_payloads;
 
 void execute_bandwith_test() {
    if (dc->numprocs() < 2) return;
    int * payload = (int *) malloc(sizeof(char) * BANDWIDTH_TEST_SIZE);
    // Fill with random data to avoid any compression (I think there isnt, but just in case)
    for(unsigned int i=0; i<BANDWIDTH_TEST_SIZE/sizeof(int); i++) {
        payload[i] = (i*17 + i); // Silly
    }
    received_payloads=0;
    // set self
    set_value(bandwidth_key(dc->procid(), dc->procid()), 0.0);

    for(int i=1; i<dc->numprocs(); i++) {
      size_t targetthisround = (dc->procid() + i) % dc->numprocs();
      // the cycle I am in is (dc->procid() + k * i) % dc->numprocs();
      // what is the smallest id in this cycle?
      // (probably related to the gcd of dc->numprocs() and i), but being lazy...
      size_t lowestid = dc->numprocs();
      size_t cyclelength = 0;

      
      for (size_t k = 0;k <= dc->numprocs(); ++k) {
        size_t id_at_pos_k = (dc->procid() + k * i) % dc->numprocs();
        lowestid = std::min(lowestid, id_at_pos_k);
        // if this is not k = 0, and cycle length has not been set
        // and I see my current id again, the cycle length must be k
        // note that the <= is important in the for loop
        if (k > 0 && id_at_pos_k == dc->procid() && cyclelength == 0) cyclelength = k;
      }
      
      // ok then, now what is my position in the cycle?
      size_t cyclepos = 0;
      for (;cyclepos < dc->numprocs(); ++cyclepos) {
        if ((lowestid + cyclepos * i) % dc->numprocs() == dc->procid()) break;
      }
      //std::cout << "stepsize: " << i << std::endl;
      //std::cout << "cycle pos: " << cyclepos << " / " << cyclelength << std::endl;
      // 3 phases are sufficient for any cycle length.
      // if even length cycle, we first go even->odd and odd -> even
      // if odd length cycle, we first go even->odd(excluding the last element in the cycle which will be odd->odd)
      //                      then even->even
      //                      then odd->even
      
      // I am send first if I my cycle position is odd and I am not the last element in the cycle
      if (cyclepos % 2 == 0 && cyclepos != cyclelength - 1) {
          dc->remote_callxs(targetthisround, distributed_metrics::remote_start_bandwidth_test, NULL, 0);
          dc->remote_callxs(targetthisround, distributed_metrics::receive_bandwidth_test_payload, payload, BANDWIDTH_TEST_SIZE);
          sleep(2);
      }
      dc->barrier();
      // send if I my cycle position is odd and I am the last element in the cycle
      if (cyclepos % 2 == 0 && cyclepos == cyclelength - 1) {
          dc->remote_callxs(targetthisround, distributed_metrics::remote_start_bandwidth_test, NULL, 0);
          dc->remote_callxs(targetthisround, distributed_metrics::receive_bandwidth_test_payload, payload, BANDWIDTH_TEST_SIZE);
          sleep(2);
      }
      dc->barrier();
      if (cyclepos % 2 == 1) {
          dc->remote_callxs(targetthisround, distributed_metrics::remote_start_bandwidth_test, NULL, 0);
          dc->remote_callxs(targetthisround, distributed_metrics::receive_bandwidth_test_payload, payload, BANDWIDTH_TEST_SIZE);
          sleep(2);
      }
      dc->barrier();
    }
    free(payload);
    

 }
 
 std::map<size_t, timer> start_times;
 
 void start_bandwidth_test(size_t from_proc) {
    start_times[from_proc] = timer();
    start_times[from_proc].start();
    std::cout << "Start bandwidth test from " << from_proc << " to " << dc->procid() << std::endl;
 }
 
 void end_bandwith_test(size_t from_proc, size_t len) {
    double time = start_times[from_proc].current_time();
    double bandwidth = (double(len)/1024.0/1024.0)/time;
    set_value(bandwidth_key(from_proc, dc->procid()), bandwidth);
    received_payloads++;
    std::cout << "Finished " << received_payloads << ". bandwidth test from " << from_proc << " to " << dc->procid() << ": " << bandwidth << " MB/sec" << std::endl; 
 }
 
 ///  METRICS

 
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
