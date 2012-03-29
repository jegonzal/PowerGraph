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



#include <cmath>
#include <cstdio>
#include <limits>
#include <iostream>
#include "graphlab.hpp"
#include "graphlab/graph/graph3.hpp"
#include "graphlab/graph/multigraph.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include "../shared/stats.hpp"

#include <graphlab/macros_def.hpp>
using namespace graphlab;

bool debug = false;
int max_iter = 50;
size_t * active_nodes_num;
size_t * active_links_num;
int iiter = 0; //current iteration
int nodes = 0;
timer gt;

enum kcore_output_fields{
  KCORE_INDEX = 1
};

struct vertex_data {
  bool active;
  unsigned char kcore;
  uint cur_links;

  vertex_data() : active(true), kcore(-1), cur_links(0)  {}

  void add_self_edge(double value) { }

  void set_val(double value, int field_type) { 
  }
  //only one output for jacobi - solution x
  double get_output(int field_type){ 
    if (field_type == KCORE_INDEX)
      return kcore;
    return -1;
  }
}; // end of vertex_data

struct edge_data {
  edge_data()  { }
  //compatible with parser which have edge value (we don't need it)
  edge_data(double val)  { }
};

typedef graphlab::multigraph<vertex_data, edge_data> multigraph_type;
typedef graphlab::graph3<vertex_data, edge_data> graph_type;

typedef vertex_data vertex_data_type;
typedef edge_data edge_data_type;
typedef edge_list edge_list_type;

graph_type * pgraph;
multigraph_type * pmultigraph;

#ifndef USE_GRAPHLAB_ENGINE
struct dummy_context{
  uint id;
  dummy_context(uint _id){ id = _id; }
  uint vertex_id(){ return id; }
  vertex_data_type & vertex_data(){ return pmultigraph->get_vertex_data(id); }
  edge_list out_edges() { return pgraph->out_edges(id); }
  edge_list in_edges() { return pgraph->in_edges(id); }
  const vertex_data_type & const_vertex_data(uint id) { return pmultigraph->get_vertex_data(id); }
//#ifdef USE_GRAPH2
  edge_data_type &edge_data(const edge_type &i){ return pgraph->edge_data(i); }
//#endif
};

#endif
 


/*
void calc_initial_degree(multigraph_type * _g, bipartite_graph_descriptor & desc){
  int active = 0;
  for (int j=0; j< _g->num_graphs(); j++){
        _g->doload(j);
        graph_type * g = _g->graph(0);
#pragma omp parallel for
     for (int i=0; i< desc.total(); i++){
       vertex_data & data = g->vertex_data(i);
       data.degree = g->out_edges(i).size() + g->in_edges(i).size();
       data.active = data.degree > 0;
       if (data.active)
         active++;
     }
     _g->unload_all();
  }
  printf("Number of active nodes in round 0 is %d\n", active);
  printf("Number of active links in round 0 is %d\n", (int)_g->num_edges());

  active_nodes_num[0] = active;
  active_links_num[0] = _g->num_edges();
}
*/
vec  fill_output(multigraph_type * g, bipartite_graph_descriptor & matrix_info, int field_type){
  vec out = zeros(matrix_info.total());
  for (int i = 0; i < matrix_info.total(); i++){
    out[i] = g->get_vertex_data(i).get_output(field_type);
  }
  return out;
}




struct kcore_update :
  public graphlab::iupdate_functor<graph_type, kcore_update> {
  void operator()(icontext_type& context) {
  } 
};

size_t num_active = 0; 
size_t links = 0;


#ifdef USE_GRAPHLAB_ENGINE
class aggregator :
  public graphlab::iaggregator<graph_type, kcore_update, aggregator> {
private:
  int num_active;
  int links;
public:
  aggregator() : num_active(0), links(0) { }
  void operator()(icontext_type& context) {
#else
  void update_function_Axb(dummy_context &context){
#endif
    uint id = context.vertex_id();
    if (id % 14000000 == 0)
      logstream(LOG_INFO)<<"Now in node: " << id << " num active: " << num_active << " links: " << links << " at time: " << gt.current_time() << std::endl;

    vertex_data & vdata = context.vertex_data();
    if (debug)
      logstream(LOG_INFO)<<"Entering node: " << id << std::endl;

    if (!vdata.active)
      return;

    //vdata.cur_links = 0;
    int increasing_links = 0;
    
    edge_list_type edges = context.out_edges();
    //edge_list_type inedgeid = context.in_edges();
   for (uint * start = edges.start_ptr + edges.abs_offset; start < edges.start_ptr + edges.abs_offset + edges._size; start++){
    //for(size_t i = 0; i < outedgeid.size(); i++) {
      const vertex_data & other = context.const_vertex_data(*start);
        if (other.active){
	  vdata.cur_links++;
          increasing_links++;
        }
    }
    /*for (size_t i =0; i < inedgeid.size(); i++){
      const vertex_data & other = context.const_vertex_data(inedgeid[i].source());
        if (other.active)
          vdata.cur_links++;
    }*/
   links += increasing_links;
   if (debug)
       logstream(LOG_INFO)<<"Exting node: " << id << " with: " << vdata.cur_links << " status: " << vdata.active << std::endl;


    };
#ifdef USE_GRAPHLAB_ENGINE
  void operator+=(const aggregator& other) { 
    num_active += other.num_active;
    links += other.links;
  }

  void finalize(iglobal_context_type& context) {
   active_nodes_num[iiter] = num_active;
   if (num_active == 0)
	links = 0;
   printf("Number of active nodes in round %d is %ld, links: %ld at time %lg\n", iiter, num_active, links, gt.current_time());
   active_links_num[iiter] = links;

   if (num_active == 0){
     context.terminate(); 
     max_iter = iiter;
   }
 }
}; // end of  aggregator
#else
  void finalize() {
   active_nodes_num[iiter] = num_active;
   if (num_active == 0)
	links = 0;
   printf("Number of active nodes in round %d is %ld, links: %ld at time %lg\n", iiter, num_active, links, gt.current_time());
   active_links_num[iiter] = links;

   if (num_active == 0){
     max_iter = iiter;
   }
 }
#endif



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  gt.start();
  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile;
  std::string format = "matrixmarket";
  int unittest = 0;
  int lineformat = MATRIX_MARKET_4;
  bool gzip = true;
  bool stats = false;
  std::string filter = "day";
  std::string listdir = "/usr2/bickson/daily.sorted/";
  std::string dirpath = "/usr2/bickson/bin.graphs/";
 
  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("lineformat", &lineformat, lineformat, "matrix line format");
  
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=TBD");
  clopts.attach_option("max_iter", &max_iter, max_iter, "maximal number of cores");
  clopts.attach_option("nodes", &nodes, nodes, "number of nodes"); 
  clopts.attach_option("gzip", &gzip, gzip, "gzipped input file?");
  clopts.attach_option("stats", &stats, stats, "calculate graph stats and exit");
  clopts.attach_option("filter", & filter, filter, "Filter - parse files starting with prefix");
  clopts.attach_option("listdir", &listdir, listdir, "Directory with list of files");
  clopts.attach_option("dirpath", &dirpath, dirpath, "Directory path");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  omp_set_num_threads(clopts.get_ncpus());

  active_nodes_num = new size_t[max_iter+1];
  active_links_num = new size_t[max_iter+1];
  memset(active_links_num, 0, sizeof(size_t)*(max_iter+1));
  memset(active_nodes_num, 0, sizeof(size_t)*(max_iter+1));

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab Linear solver library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Gaussian Belief Propagation, "
    << "Jacobi method, Conjugate Gradient" << std::endl;



  // Create a core
  graphlab::core<graph_type, kcore_update> core;
  core.set_options(clopts); // Set the engine options

  //unit testing
  if (unittest == 1){
     datafile = "kcores_unittest1";
  }

  std::cout << "Load graph" << std::endl;
  bipartite_graph_descriptor matrix_info;

  //nodes = 123306178;
  //nodes=149747010;
  //nodes = 121408373;
  //matrix_info.rows = matrix_info.cols = nodes;
  //matrix_info.nonzeros = 1000000000;
  /* Rows:      95526
 * Cols:      3561
 * Nonzeros:  3298163
 */
  //matrix_info.rows = 95526;  matrix_info.cols = 3561; matrix_info.nonzeros = 3298163;
  //std::string dirpath="/mnt/bigbrofs/usr0/bickson/out_phone_calls/";
 //core.graph().set_undirected();
  core.set_scope_type("null");

    
    multigraph_type multigraph;
    multigraph.load(listdir, dirpath, filter, true);
    pmultigraph = &multigraph;
  if (stats){
    calc_multigraph_stats_and_exit<multigraph_type>(&multigraph, matrix_info);
  }
  graphlab::timer mytimer; mytimer.start();
  bool single_graph = (pmultigraph->num_graphs() == 1);
  bool first_time = true;

  int pass = 0;
  for (iiter=1; iiter< max_iter+1; iiter++){
    logstream(LOG_INFO)<<mytimer.current_time() << ") Going to run k-cores iteration " << iiter << " at time: " << gt.current_time() << std::endl;
    while(true){
      int prev_nodes = active_nodes_num[iiter];
      num_active = 0; links = 0;

      for (int i=0; i< multigraph.num_graphs(); i++){
        if (!single_graph || first_time)
           multigraph.doload(i, true, true, true);
       core.graph() = *multigraph.graph(0);
       matrix_info.nonzeros = core.graph().num_edges();
       matrix_info.rows = matrix_info.cols = core.graph().num_vertices();
       pgraph = multigraph.graph(0);
#ifdef USE_GRAPHLAB_ENGINE
       aggregator acum;
       core.add_aggregator("sync", acum, 1000);
       core.add_global("NUM_ACTIVE", int(0));
       glcore->aggregate_now("sync"); 
#else
#pragma omp parallel for
       for (int t=0; t< matrix_info.total(); t++){
            if (pmultigraph->get_vertex_data(t).active){
              dummy_context con(t);
              update_function_Axb(con);
            }
       }
       if (!single_graph)
          multigraph.unload_all(); 
       else first_time = false;
      }
      int t;
#pragma omp parallel for private(t) reduction(+: num_active)
     for (t=0; t< matrix_info.total(); t++){
         vertex_data& vdata = pmultigraph->get_vertex_data(t);
        if (vdata.cur_links <= iiter){
          vdata.active = false;
          vdata.kcore = iiter;
	//links -= (outedgeid.size() + inedgeid.size());
        }
       vdata.cur_links = 0;
       if (vdata.active)
         num_active = num_active + 1;

     }
 
      finalize();
#endif
      pass++;
      int cur_nodes = active_nodes_num[iiter];
      if (prev_nodes == cur_nodes)
        break; 
    }
    if (active_nodes_num[iiter] == 0)
	break;
  }
 
  std::cout << "KCORES finished in " << mytimer.current_time() << std::endl;
  std::cout << "Number of updates: " << pass*core.graph().num_vertices() << " pass: " << pass << std::endl;
  imat retmat = imat(max_iter+1, 4);
  memset((int*)data(retmat),0,sizeof(int)*retmat.size());

  std::cout<<active_nodes_num<<std::endl;
  std::cout<<active_links_num<<std::endl;
  active_nodes_num[0] = pmultigraph->num_vertices();
  active_links_num[0] = pmultigraph->num_edges();
  set_val(retmat, 0, 1, active_nodes_num[0]);
  set_val(retmat, 0, 2, active_nodes_num[0]);
  set_val(retmat, 0, 3, active_links_num[0]);

  for (int i=0; i <= max_iter; i++){
    set_val(retmat, i, 0, i);
    if (i >= 1){
      set_val(retmat, i, 1, active_nodes_num[i-1]-active_nodes_num[i]);
      set_val(retmat, i, 2, active_nodes_num[0]-active_nodes_num[i]);
      set_val(retmat, i, 3, active_links_num[i]);
    }
  } 
  //write_output_matrix(datafile + ".kcores.out", format, retmat);
  std::cout<<retmat<<std::endl;


  vec ret = fill_output(pmultigraph, matrix_info, KCORE_INDEX);
  write_output_vector(datafile + "x.out", format, ret,false);

  if (unittest == 1){
    imat sol = init_imat("0 0 0 0; 1 1 1 1; 2 4 5 7; 3 4 9 13", 4, 4);
    assert(sumsum(sol - retmat) == 0);
  }


   return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
