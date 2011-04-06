
/**
 * Metrics tool for Graphlab.
 */
  
#ifndef GRAPHLAB_METRICS_BASE
#define GRAPHLAB_METRICS_BASE

#include <cassert>
#include <cstring>
#include <map>
#include <vector>

#include <graphlab/util/timer.hpp> 
#include <graphlab/parallel/pthread_tools.hpp> 


namespace graphlab {

  class imetrics_reporter;
    
  enum metrictype {REAL, INTEGER, TIME, STRING, VECTOR};
    
  // Data structure for storing metric entries
  // NOTE: This data structure is not very optimal, should
  // of course use inheritance. But for this purpose,
  // it works fine as the number of metrics entry is small.
  struct metrics_entry {
    int count;
    double value;
    double minvalue;
    double cumvalue;
    double maxvalue;
    metrictype valtype;
    std::string stringval;
    std::vector<double> v;
    timer _timer;
        
    metrics_entry() {} 
        
    metrics_entry(double firstvalue, metrictype _valtype) {
      minvalue = firstvalue;
      maxvalue = firstvalue;
      value = firstvalue;
      valtype = _valtype;
      cumvalue = value;
      count = 1;
      if (valtype == VECTOR) v.push_back(firstvalue);
    };
    metrics_entry(std::string svalue) {
      valtype = STRING;
      stringval = svalue;
    }
    metrics_entry(metrictype _valtype) {
      valtype = _valtype;
      count = 0;
      cumvalue = 0;
    }
    void adj(double v) {
      if (count == 0) {
        minvalue = v;
        maxvalue = v;
      } else {
        minvalue = std::min(v,minvalue);
        maxvalue = std::max(v,maxvalue);
      }
    }
    
    
    void add(double x) {
      adj(x);
      value += x;
      cumvalue += x;
      count++;
      if (valtype == VECTOR) v.push_back(x);
    }
    void set(double v) {
      adj(v);
      value = v;
      cumvalue += v;
      count++;
    }
    void set(std::string s) {
      stringval = s;
    }
    void timer_start() {
      _timer.start();
    }
    void timer_stop() {
      add(_timer.current_time());
    }
  };

  // Metrics instance. Create using method create_metrics_instance()
  class metrics {
    
    std::string name, ident;
    std::map<std::string, metrics_entry> entries;
    mutex safelock;
        
  protected:
    metrics(std::string _name, std::string _id) : name(_name), ident (_id) {
    }
        
  private:
    metrics(const metrics&);          // No definition
    void operator=(const metrics &); // No definition
        
  public: 
            
    /**
     * Add to an existing value or create new.
     */
    void add(std::string key, double value, metrictype type = REAL) {
      safelock.lock();
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(value, type);
      } else {
        entries[key].add(value);
      }
      safelock.unlock();
    }
    
    void add_vector(std::string key, double value) {
       safelock.lock();
       if (entries.count(key) == 0) {
         entries[key] = metrics_entry(value, VECTOR);
       } else {
         entries[key].add(value);
       }
       safelock.unlock(); 
    }
        
    void set(std::string key, size_t value) {
      set(key, value, INTEGER);
    }
        
    void set(std::string key, double value, metrictype type = REAL) {
      safelock.lock();
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(value, type);
      } else {
        entries[key].set(value);
      }
      safelock.unlock(); 
    }
    void set(std::string key, std::string s) {
      safelock.lock();
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(s);
      } else {
        entries[key].set(s);
      }
      safelock.unlock(); 
    }
        
    void start_time(std::string key) {
      safelock.lock();
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(TIME);
      } 
      entries[key].timer_start();
      safelock.unlock(); 
    }
        
    void stop_time(std::string key) {
      safelock.lock();
      assert(entries.count(key) > 0);
      entries[key].timer_stop();
      safelock.unlock(); 
    }
        
    metrics_entry get(std::string key) {
      assert(entries.count(key) > 0);
      return entries[key];
    }
        
    void report(imetrics_reporter & reporter);
        

  public:
        
    /* USE THIS TO CREATE A NEW METRICS INSTANCE */
        
    /**
     * Creates a new named metrics instance, if it does not exist yet.
     * If it exist, and flag create_always_new is FALSE, previously created
     * instance is returned. Otherwise, a new instance with suffix ":n", where
     * n is the number of such named instances is created.
     */
    static metrics & create_metrics_instance(std::string name, bool create_always_new);
    static metrics & create_metrics_instance(std::string name);
    static void clear_all_metrics_instances();
    static void report_all(imetrics_reporter & reporter);
  };


    



};






#endif
