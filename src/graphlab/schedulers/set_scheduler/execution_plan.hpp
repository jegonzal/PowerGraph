#ifndef EXECUTION_PLAN_HPP
#define EXECUTION_PLAN_HPP
#include <map>
#include <set>
#include <vector>
#include <iostream>

#include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
#include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/parallel/atomic.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  template<typename Graph>
  class execution_plan {
  public:
    typedef Graph graph_type;

    typedef update_task<Graph> update_task_type;

    typedef typename update_task_type::update_function_type
    update_function_type;
    
    typedef ivertex_set<Graph> ivertex_set_type;

  public:
    struct plan_vertex {
      unsigned short setid;
      vertex_id_t v;
    };

    struct plan_edge { };

    struct cpu_plan_node {
      int waitid;
      vertex_id_t v;
      update_function_type updatefunc;
      std::vector<size_t> decrementid;
    };
    

    typedef graph<plan_vertex, plan_edge> plan_graph_type;

    execution_plan() {}

    void execute(ivertex_set_type &srcset, update_function_type u) {
      planorder.push_back(&srcset);
      updates.push_back(u);
    }
  
    void generate_plan(const Graph& g, size_t ncpus) {
      // plan graph generated
      plan_graph_type plangraph;
      build_plan_graph(g, plangraph);
      //build_metis_plan(g,plangraph, ncpus);
      build_topsort_plan(g, plangraph, ncpus);
      atomics_reduction();
    }
  
  
    void build_topsort_plan(const Graph& g, 
                            plan_graph_type& plangraph, 
                            size_t ncpus) {
      std::vector<vertex_id_t> topsort;
      if (plangraph.topological_sort(topsort) == false) {
        assert(false);
      }
      assert(topsort.size() == plangraph.num_vertices());
      cpuplans.resize(ncpus);
      // construct the partitioning
      std::vector<uint32_t> partitionids;
      partitionids.reserve(topsort.size());
      for (size_t i = 0; i < topsort.size(); ++i) {
        partitionids[i] = i % ncpus;
      }
      build_plan_from_partition(g, plangraph,ncpus,partitionids);    
    }
  
    void build_metis_plan(const Graph& g, 
                          plan_graph_type& plangraph, 
                          size_t ncpus) {
      // partition to minimize edges.
      std::vector<uint32_t> partitionids;
      plangraph.metis_partition(ncpus, partitionids);
      build_plan_from_partition(g, plangraph,ncpus,partitionids);

    }
  
    void build_plan_from_partition(const Graph& g, 
                                   plan_graph_type& plangraph, 
                                   size_t ncpus,
                                   std::vector<uint32_t> &partitionids) {
      cpuplans.resize(ncpus);
      // plan graph creation necessitates that all edges
      // go from a lower id to a higher id
      std::vector<vertex_id_t> atomic_to_planvertex;
      std::map<vertex_id_t, size_t> planvertex_to_atomic;
    
      for (vertex_id_t curplanv = 0; 
           curplanv < plangraph.num_vertices(); ++curplanv) {
        cpu_plan_node plannode;
        // do I need to wait for anyone?
        if (planvertex_to_atomic.find(curplanv) != planvertex_to_atomic.end()) {
          plannode.waitid = planvertex_to_atomic[curplanv];
        } else {
          plannode.waitid = -1;
        }
        plannode.updatefunc = 
          updates[plangraph.vertex_data(curplanv).setid];
        plannode.v = plangraph.vertex_data(curplanv).v;
        // look at the forward edges
        foreach(edge_id_t eid, plangraph.out_edge_ids(curplanv)) {
          vertex_id_t destplanv = plangraph.target(eid);
          // if I cross a partition, I need to create an atomic
          // otherwise nothing need to be done
          if (partitionids[destplanv] !=partitionids[curplanv]) {
            size_t atomicid = 0;
            // there already exists an atomic, increase the initial value
            if (planvertex_to_atomic.find(destplanv) != planvertex_to_atomic.end()) {
              atomicid = planvertex_to_atomic[destplanv];
              atomics_initial_value[atomicid]++;
            } else {
              // I need to create a new atomic
              atomicid = atomics_initial_value.size();
              // update all the parallel array
              atomics_initial_value.push_back(1);     // initial value of 1
              atomic_to_planvertex.push_back(destplanv); 
              planvertex_to_atomic[destplanv] = atomicid; 
            }
            plannode.decrementid.push_back(atomicid);
          }
        }
        cpuplans[partitionids[curplanv]].push_back(plannode);
      }
    
      atomics.resize(atomics_initial_value.size());
      cpu_curloc.resize(ncpus);
    }
  
    void atomics_reduction() {
      // if I ever increment a particular atomic more than once,
      // I only need to do it the last time.
      for (size_t i = 0;i < cpuplans.size(); ++i) {
        // get the set of atomics I touch
        std::set<size_t> atomicstouched;
        for (int j = cpuplans[i].size() - 1; j >=0; --j) {
          for (int k = cpuplans[i][j].decrementid.size() - 1; k >= 0 ; --k) {
            size_t atomicid = cpuplans[i][j].decrementid[k];
            if (atomicstouched.find(atomicid) != atomicstouched.end()) {
              //erase the decrement
              cpuplans[i][j].decrementid[k] = cpuplans[i][j].decrementid[cpuplans[i][j].decrementid.size() - 1];
              cpuplans[i][j].decrementid.resize(cpuplans[i][j].decrementid.size() - 1);
              // edecrease the atomic strength
              atomics_initial_value[atomicid]--;
            }
            else {
              atomicstouched.insert(cpuplans[i][j].decrementid[k]);
            }
          }
        }
      }
    }
    void init_plan() {
      for (size_t i = 0;i < atomics_initial_value.size(); ++i) {
        atomics[i].value = atomics_initial_value[i];
      }
      for (size_t i = 0; i <cpu_curloc.size(); ++i) cpu_curloc[i] = 0;
    }
  
    bool get_next_task(size_t cpuid, update_task_type &u) {
      size_t &curloc = cpu_curloc[cpuid];
      // we might be done
      if (curloc >= cpuplans[cpuid].size()) return false;
      cpu_plan_node& plan = cpuplans[cpuid][curloc];
      //wait on what we need to wait. if we need to wait
      int waitid = plan.waitid;
      if (waitid >= 0) {
        volatile size_t* i = &(atomics[waitid].value);
        if ((*i) != 0) {
          while(true) {
            if ((*i) == 0) {
              break;
            }
          }
        }
      }
      u = update_task_type(plan.v, plan.updatefunc);
      return true;
    }
  
    void completed_task(size_t cpuid) {
      size_t &curloc = cpu_curloc[cpuid];
      cpu_plan_node& plan = cpuplans[cpuid][curloc]; 
      // trigger everything downstream
      for (size_t i = 0; i < plan.decrementid.size(); ++i) {
        atomics[plan.decrementid[i]].dec();
      }
      ++curloc;
    }
  
    bool done() {
      for (size_t i = 0;i < cpu_curloc.size(); ++i) {
        if (cpu_curloc[i] < cpuplans[i].size()) return false;
      }
      return true;
    }
  
    void print_plan(size_t cpuid) {
      for (size_t i = 0 ;i < cpuplans[cpuid].size(); ++i) {
        if (cpuplans[cpuid][i].waitid >= 0) {
          std::cout << "Wait(" << cpuplans[cpuid][i].waitid << ")\n";
        }
        std::cout << "Execute(" << cpuplans[cpuid][i].v << ")\n";
        for (size_t j = 0;j < cpuplans[cpuid][i].decrementid.size(); ++j) {
          std::cout << "Dec(" << cpuplans[cpuid][i].decrementid[j] << ")\n";
        }
      }
    }
  private:
    std::vector<ivertex_set_type*> planorder;
    std::vector<update_function_type> updates;
  
    std::vector<std::vector<cpu_plan_node> > cpuplans;
    std::vector<atomic<size_t> > atomics;
    std::vector<size_t> atomics_initial_value;
    std::vector<size_t> cpu_curloc;
  
    void build_plan_graph(const Graph& g,
                          plan_graph_type& plangraph) {
      size_t numv = g.num_vertices();
      // this vector contains the ID of the last set containing this vertex
      std::vector<short> last_set_with_vertex;
      last_set_with_vertex.assign(numv, -1);
      // thios vector contains the vertexID of the last planvertex which
      // contains this graph vertex
      std::vector<vertex_id_t> last_planvertex_of_vertex;
      last_planvertex_of_vertex.assign(numv, -1);
    
      plan_vertex p;
      for (int i = 0; i < (int)(planorder.size()); ++i) {
        p.setid = i;
        // current set of vertices to insert
        ss_set_type s = planorder[i]->get();
        ss_set_type_iterator j = begin(s);
        // iterate through all the vertices in this set
        while(j!=end(s)) {
          // insert the vertex into the plan graph
          vertex_id_t curv = *j;
          p.v = curv;
          vertex_id_t curplanv =  plangraph.add_vertex(p) ;
        
          // one special case. If I am already in the plan, and my 
          // last set is higher than all my neighbors, I just need to add edge to 
          // myself
          bool specialcase = true;
          foreach(edge_id_t eid, g.in_edge_ids(curv)) {
            vertex_id_t srcv = g.source(eid);
            if (last_set_with_vertex[srcv] >= last_set_with_vertex[curv]) {
              specialcase = false; 
              break;
            }
          }
          foreach(edge_id_t eid, g.out_edge_ids(curv)) {
            vertex_id_t destv = g.target(eid);
            if (last_set_with_vertex[destv] >= last_set_with_vertex[curv]) {
              specialcase = false; 
              break;
            }
          }

          // if special case hit, we just add an edge from the
          // last occurance of my vertex in the plan
          // otherwise we have to do the complete construction
          if (specialcase && last_set_with_vertex[curv] != -1
              && last_set_with_vertex[curv] < i) {
            vertex_id_t srcplanv = last_planvertex_of_vertex[curv];
            plangraph.add_edge(srcplanv, curplanv);
          } else {
            // look at the neighbors and build edges to them        
            std::set<vertex_id_t> visited;
            foreach(edge_id_t eid, g.in_edge_ids(curv)) {
              vertex_id_t srcv = g.source(eid);
              visited.insert(srcv);
              vertex_id_t srcplanv = last_planvertex_of_vertex[srcv] ;
              // add an edge if I have been inserted before
              // and I am in a previous set
              if (srcplanv != vertex_id_t(-1) && last_set_with_vertex[srcv] < i) {
                plangraph.add_edge(srcplanv, curplanv);
              }
            }
            foreach(edge_id_t eid, g.out_edge_ids(curv)) {
              vertex_id_t destv = g.target(eid);
              if (visited.count(destv) == 0) {
                vertex_id_t destplanv = last_planvertex_of_vertex[destv] ;
                // add an edge if I have been inserted before
                // and I am in a previous set
                if (destplanv != vertex_id_t(-1) && last_set_with_vertex[destv] < i) {
                  plangraph.add_edge(destplanv, curplanv);
                }
              }
            }
            if (last_planvertex_of_vertex[curv] != vertex_id_t(-1)) {
              plangraph.add_edge(last_planvertex_of_vertex[curv], curplanv);
            }
          }
          // update the plan data structure
          last_set_with_vertex[curv] = i;
          last_planvertex_of_vertex[curv] = curplanv;
          ++j;
        }
      }    
    }
  };

}
#include <graphlab/macros_undef.hpp>
#endif
