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

  vertex_data() : active(true), kcore(-1){}

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
    vertex_data & vdata = context.vertex_data();
    if (!vdata.active)
      return;

    int cur_iter = iiter;
    int cur_links = 0;
    
    const edge_list_type outedgeid = context.out_edges();
    const edge_list_type inedgeid = context.in_edges();

    for(size_t i = 0; i < outedgeid.size(); i++) {
      const vertex_data & other = context.const_vertex_data(outedgeid[i].target());
        if (other.active) {
	  cur_links++;
        }
    }
    for (size_t i =0; i < inedgeid.size(); i++){
      const vertex_data & other = context.const_vertex_data(inedgeid[i].source());
        if (other.active)
          cur_links++;
    }
    if (cur_links <= cur_iter){
        vdata.active = false;
        vdata.kcore = cur_iter;
        for(size_t i = 0; i < outedgeid.size(); i++) {
          const vertex_data & other = context.const_vertex_data(outedgeid[i].target());
           if (other.active){
             context.schedule(outedgeid[i].target(), kcore_update());
           }
        }
        for (size_t i =0; i < inedgeid.size(); i++){
          const vertex_data & other = context.const_vertex_data(inedgeid[i].source());
           if (other.active)
	     context.schedule(inedgeid[i].source(), kcore_update());    
        }
    }
  };

};

class aggregator :
  public graphlab::iaggregator<graph_type, kcore_update, aggregator> {
private:
  int num_active;
  int links;
public:
  aggregator() : num_active(0), links(0) { }

  void operator()(icontext_type& context) {
    vertex_data & vdata = context.vertex_data();
    if (!vdata.active) return;
    int increasing_links = 0;
    const edge_list_type outedgeid = context.out_edges();
    for(size_t i = 0; i < outedgeid.size(); i++) {
      const vertex_data & other = context.const_vertex_data(outedgeid[i].target());
        if (other.active){
          increasing_links++;
        }
    }
    links += increasing_links;
    num_active++;
  };

  void operator+=(const aggregator& other) { 
    num_active += other.num_active;
    links += other.links;
  }

  void finalize(iglobal_context_type& context) {
   active_nodes_num[iiter] = num_active;
   if (num_active == 0)
	links = 0;
   printf("Number of active nodes in round %d is %di, links: %d\n", iiter, num_active, links);
   active_links_num[iiter] = links;

   if (num_active == 0){
     context.terminate(); 
     max_iter = iiter;
   }
 }
}; // end of  aggregator





int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  gt.start();
  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile;
  std::string format = "matrixmarket";
  int unittest = 0;
  int lineformat = MATRIX_MARKET_3;
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
  
    multigraph_type multigraph;
    multigraph.load(listdir, dirpath, filter, true);
    pmultigraph = &multigraph;
  if (stats){
    calc_multigraph_stats_and_exit<multigraph_type>(&multigraph, matrix_info);
  }
 
  matrix_info.rows = matrix_info.cols = nodes;
  core.set_scope_type("none");
  assert(in_files.size() > 0);
  for (int i=0; i< std::min(max_files, (int)in_files.size()); i++){
      graphlab::timer mt; mt.start();
      core.graph().load_directed(dirpath + in_files[i], false, false);

    matrix_info.nonzeros = core.graph().num_edges();

    //DB: to be cleaned later
    if (datafile == "smallnetflix"){
      matrix_info.rows = 95526; matrix_info.cols = 3561; matrix_info.nonzeros= 3298163;
    }
    else if (datafile == "netflixe"){
      matrix_info.rows = 480189; matrix_info.cols = 17770; matrix_info.nonzeros= 99072112;
    }
      
    logstream(LOG_INFO)<<mt.current_time() << "Takes to load graph" << std::endl;
  } 


  core.schedule_all(kcore_update());
 
  aggregator acum;
  core.add_aggregator("sync", acum, 1000000000);
  core.add_global("NUM_ACTIVE", int(0));


    multigraph_type multigraph;
    multigraph.load(listdir, dirpath, filter, true);
    pmultigraph = &multigraph;
  if (stats){
    calc_multigraph_stats_and_exit<multigraph_type>(&multigraph, matrix_info);
  }
 
  graphlab::timer mytimer; mytimer.start();

  for (iiter=1; iiter< max_iter+1; iiter++){
    logstream(LOG_INFO)<<mytimer.current_time() << ") Going to run k-cores iteration " << iiter << std::endl;
      core.schedule_all(kcore_update());
      core.start();
      core.aggregate_now("sync");
      if (active_nodes_num[iiter] == 0)
	break;
  }
 
  std::cout << "KCORES finished in " << mytimer.current_time() << std::endl;
  std::cout << "Number of updates: " << core.last_update_count() << " per node: " << ((double)core.last_update_count())/core.graph().num_vertices() << std::endl;

  imat retmat = imat(max_iter+1, 4);
  memset((int*)data(retmat),0,sizeof(int)*retmat.size());

  std::cout<<active_nodes_num<<std::endl;
  std::cout<<active_links_num<<std::endl;

  for (int i=0; i <= max_iter; i++){
    set_val(retmat, i, 0, i);
    if (i >= 1){
      set_val(retmat, i, 1, active_nodes_num[i-1]-active_nodes_num[i]);
      set_val(retmat, i, 2, active_nodes_num[0]-active_nodes_num[i]);
      set_val(retmat, i, 3, core.graph().num_edges() - active_links_num[i]);
    }
  } 
  //write_output_matrix(datafile + ".kcores.out", format, retmat);
  std::cout<<retmat<<std::endl;


  vec ret = fill_output(&core.graph(), matrix_info, KCORE_INDEX);
  write_output_vector(datafile + "x.out", format, ret,false);


  if (unittest == 1){
    imat sol = init_imat("0 0 0 0; 1 1 1 1; 2 4 5 7; 3 4 9 13", 4, 4);
    assert(sumsum(sol - retmat) == 0);
  }

   return EXIT_SUCCESS;
}



