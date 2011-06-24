/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */



/**
 * Metrics tool for Graphlab.
 */
  
#ifndef GRAPHLAB_METRICS_BASE
#define GRAPHLAB_METRICS_BASE

#include <cstring>
#include <map>
#include <vector>
#include <limits>
#include <graphlab/logger/assertions.hpp>

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
    size_t count;
    double value;
    double minvalue;
    double cumvalue;
    double maxvalue;
    metrictype valtype;
    std::string stringval;
    std::vector<double> v;
    timer _timer;
        
    metrics_entry() {} 
        
    inline metrics_entry(double firstvalue, metrictype _valtype) {
      minvalue = firstvalue;
      maxvalue = firstvalue;
      value = firstvalue;
      valtype = _valtype;
      cumvalue = value;
      count = 1;
      if (valtype == VECTOR) v.push_back(firstvalue);
    };
    inline metrics_entry(std::string svalue) {
      valtype = STRING;
      stringval = svalue;
    }
    inline metrics_entry(metrictype _valtype) {
      valtype = _valtype;
      count = 0;
      cumvalue = 0;
      value = 0;
      minvalue = std::numeric_limits<double>::max();
      maxvalue = std::numeric_limits<double>::min();
    }
    inline void adj(double v) {
      if (count == 0) {
        minvalue = v;
        maxvalue = v;
      } else {
        minvalue = std::min(v,minvalue);
        maxvalue = std::max(v,maxvalue);
      }
    }
    
    
    inline void add(double x) {
      adj(x);
      value += x;
      cumvalue += x;
      if (valtype == VECTOR) {
        ++count;
        v.push_back(x);
      }
    }


    inline void set(double v) {
      adj(v);
      value = v;
      cumvalue += v;
    }
    inline void set(std::string s) {
      stringval = s;
    }

    inline void add_vector_entry(size_t i, double x) {
      ASSERT_EQ(valtype, VECTOR);
      if (v.size() < i + 1) v.resize(i + 1);
      count = v.size();
      value += x;
      cumvalue += x;
      v[i] += x;
      adj(v[i]);
    }
    
    inline void set_vector_entry(size_t i, double x) {
      ASSERT_EQ(valtype, VECTOR);
      if (v.size() < i + 1) v.resize(i + 1);
      count = v.size();
      value = value - v[i] + x;
      cumvalue = cumvalue - v[i] + x;
      v[i] = x;

      minvalue = x; maxvalue = x;
      for (size_t i = 0; i < v.size(); ++i) {
        adj(v[i]);
      }
    }
    
    inline void timer_start() {
      _timer.start();
    }
    inline void timer_stop() {
      add(_timer.current_time());
    }
  };

  /**
   * Metrics instance for logging metrics of a single object type.
   * Name of the metrics instance is set on construction.
   */
  class metrics {
    
    std::string name, ident;
    std::map<std::string, metrics_entry> entries;
        
  public: 
    inline metrics(std::string _name = "", std::string _id = "") : name(_name), ident (_id) {
    }

    inline void clear() {
      entries.clear();
    }
    /**
     * Add to an existing value or create new.
     */
    inline void add(std::string key, double value, metrictype type = REAL) {
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(value, type);
      } else {
        entries[key].add(value);
      }
    }
    
    inline void add_to_vector(std::string key, double value) {
       if (entries.count(key) == 0) {
         entries[key] = metrics_entry(value, VECTOR);
       } else {
         entries[key].add(value);
       }
    }

    inline void add_vector_entry(std::string key, size_t idx, double value) {
       if (entries.count(key) == 0) {
         entries[key] = metrics_entry(VECTOR);
       }
       entries[key].add_vector_entry(idx, value);
    }
    
    inline void set(std::string key, size_t value) {
      set(key, (double)value, INTEGER);
    }
        
    inline void set(std::string key, double value, metrictype type = REAL) {
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(value, type);
      } else {
        entries[key].set(value);
      }
    }
    
    inline void set_integer(std::string key, size_t value) {
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry((double)value, INTEGER);
      } else {
        entries[key].set((double)value);
      }
    }
    
    inline void set(std::string key, std::string s) {
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(s);
      } else {
        entries[key].set(s);
      }
    }

    inline void set_vector_entry_integer(std::string key, size_t idx, size_t value) {
      set_vector_entry(key, idx, (double)(value));
    }
    
    inline void set_vector_entry(std::string key, size_t idx, double value) {
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(VECTOR);
      } 
      entries[key].set_vector_entry(idx, value);
    }
    
    inline void start_time(std::string key) {
      if (entries.count(key) == 0) {
        entries[key] = metrics_entry(TIME);
      } 
      entries[key].timer_start();
    }
        
    inline void stop_time(std::string key) {
      ASSERT_GE(entries.count(key), 0);
      entries[key].timer_stop();
    }
        
    inline metrics_entry get(std::string key) {
      ASSERT_GE(entries.count(key), 0);
      return entries[key];
    }
        
    void report(imetrics_reporter & reporter);
  };



};



#include <graphlab/metrics/imetrics_reporter.hpp>


#endif

