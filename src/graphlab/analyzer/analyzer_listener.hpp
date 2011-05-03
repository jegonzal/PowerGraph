/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef analyzer_listener_HPP
#define analyzer_listener_HPP

#include <vector>
#include <string>

#include <graphlab/app_support/appstats.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

namespace graphlab {
  
  // Sorry about the naming, but I like it.
  typedef std::vector<double>(*global_dumper)(blob_graph  &g);
  
  class analyzer_listener : 
    public imonitor<blob_graph> {
    
    typedef imonitor<blob_graph>::update_task_type
    update_task_type;

    typedef imonitor<blob_graph>::iengine_type iengine_type;

    blob_graph* g;
    std::vector<vertex_stats> stats;
    timer t;
    
    std::string appname;
    
    /* Locks */
    mutex startlock;
    conditional startcond;
    
    bool pause;
    
    global_dumper dumperfunc;
    int dumpfreq;
    std::vector<std::string> dumpheaders;
    FILE * dumpfile;
    
    atomic<size_t> taskcount;
    
    
  public:
        
    analyzer_listener(blob_graph* _g, std::string _appname) : 
      stats(_g->num_vertices()),taskcount(0) {      
      g = _g;
      appname = _appname;
      for(unsigned int i=0; i<g->num_vertices(); i++) {
        stats[i].vertexid = i;
        stats[i].priority = 0.0f;
        stats[i].last_update_time = 0.0f;
        stats[i].value = 0.0f;
        stats[i].updatecount = 0;
      }
      pause = false;
      dumperfunc = NULL;
    }
    
    ~analyzer_listener() {
    } 
    
    virtual void set_global_dumper(std::vector<std::string> headers, global_dumper dumpfunc, int dumpfreq_updates=5000) {
      dumperfunc = dumpfunc;
      dumpfreq = dumpfreq_updates;
      dumpheaders = headers;
      
      /* Open file */
      std::string filename = appname + "_dump.dat";
      dumpfile = fopen(filename.c_str(), "w");
      
      /* Write headers */
      fprintf(dumpfile, "updates\t");
      for(unsigned int i=0; i<headers.size(); i++) {
        fprintf(dumpfile, "\t");
        fprintf(dumpfile, "%s", headers[i].c_str());
      }
      fprintf(dumpfile,"\n");
    }
    
    virtual void write_graph_datafile() {
      std::string filename = appname + "_graph.dat";
      FILE * f = fopen(filename.c_str(), "w");      
      fprintf(f, "indegree\toutdegree\n");
      for(vertex_id_t vid = 0; vid< g->num_vertices(); vid++) {
        fprintf(f, "%u\t%u\n", (unsigned int)g->num_in_neighbors(vid), 
                (unsigned int)g->num_out_neighbors(vid));
      }
      fclose(f);
      
      printf("Wrote %s \n", filename.c_str());
    }
    
    virtual void write_final_stats() {
      std::string filename = appname + "_updates.dat";
      FILE * f = fopen(filename.c_str(), "w");      
      fprintf(f, "num_updates\n");
      
      for(vertex_id_t vid = 0; vid< g->num_vertices(); vid++) {
        fprintf(f, "%u\n", (unsigned int)stats[vid].updatecount);
      }
      
      fclose(f);
    }
    
    
    virtual void dump(size_t tc) {
      std::vector<double> newdump = dumperfunc(*g);
      fprintf(dumpfile, "%u", (unsigned int) tc);
      for(unsigned int i=0; i<newdump.size(); i++) {
        fprintf(dumpfile, "\t%lf", newdump[i]); 
      }
      fprintf(dumpfile, "\n");
    }
    
    virtual void init(iengine_type* engine) {
      printf("====== ANALYZER SERVER LISTENER STARTED ===== \n");
      // Todo: start web server and wait for connection.
      
      write_graph_datafile();
    }
    
    /**
     * Pause & continue hacks
     */
    virtual void pauseExec() {
      startlock.lock();
      pause = true;
      startlock.unlock();
    }
    
    virtual void continueExec() {
      startlock.lock();
      pause = false;
      startcond.broadcast();
      startlock.unlock();
    }
    
    
    void scheduler_task_scheduled(update_task_type task, double current_max_priority) {
      vertex_id_t vid = task.vertex();
      stats[vid].updatecount++;
      
      /* Hack to pause - don't return control to scheduler/engine */
      while(pause) {

        startlock.lock();
        if (pause) {
          startcond.wait(startlock);
        }
        startlock.unlock();
      }
      
      /* Dump */
      if (dumperfunc != NULL) {
        size_t tc = taskcount.inc();
        if (tc % dumpfreq == 0) {
          /* Pause - there can be some odd updates, but maybe this is ok anyway */
          pauseExec();
          dump(tc);
          continueExec();
        }
      }
      
    }
    
    
    std::vector<vertex_stats> get_stats_snapshot() {
      return stats;
    }
    
    /* Scheduler calls */
    virtual void scheduler_task_added(update_task_type task, 
                                      double priority) { }
    
    virtual void scheduler_task_promoted(update_task_type task, 
                                         double diffpriority, 
                                         double totalpriority) {  }  
    
    
    
    /* Application calls */
    virtual void app_set_vertex_value(vertex_id_t vid, double value) { 
      stats[vid].value = (float)value;
    }
    
    /* Called by application to help visualizers to scale values properly */
    virtual void app_set_vertex_value_scale(double min, double max) { }
    
    virtual void engine_worker_dies(size_t cpuid, int tcc) { 
      if (cpuid == 0) {
        write_final_stats();
        size_t tc = taskcount.inc()-1;
        dump(tc);
        //DB fclose(dumpfile); //potential problem in a loop
      }
    }
    
    
  };
  
  
}

#endif
