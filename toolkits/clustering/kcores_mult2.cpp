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
ivec active_nodes_num;
ivec active_links_num;
int iiter = 0; //current iteration
int nodes = 0;
uint * edge_count;
unsigned long long total_edges = 0;

enum kcore_output_fields{
  KCORE_INDEX = 1
};

struct vertex_data {
  bool active;
  int kcore, degree;

  vertex_data() : active(true), kcore(-1), degree(0)  {}

  void add_self_edge(double value) { }

  void set_val(double value, int field_type) { 
  }
  //only one output for jacobi - solution x
  double get_output(int field_type){ 
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
graph_type * reference_graph = NULL;

void calc_initial_degree(multigraph_type * _g, bipartite_graph_descriptor & desc){
  int active = 0;
  for (int j=0; j< _g->num_graphs(); j++){
    graph_type * g = _g->graph(j);
  for (int i=0; i< desc.total(); i++){
     vertex_data & data = g->vertex_data(i);
     data.degree = g->out_edges(i).size() + g->in_edges(i).size();
     data.active = data.degree > 0;
     if (data.active)
       active++;
  }
  }
  printf("Number of active nodes in round 0 is %d\n", active);
  printf("Number of active links in round 0 is %d\n", (int)_g->num_edges());

  active_nodes_num[0] = active;
  active_links_num[0] = _g->num_edges();
}




struct kcore_update :
  public graphlab::iupdate_functor<graph_type, kcore_update> {
  void operator()(icontext_type& context) {
  } 
};

class accumulator :
  public graphlab::iaccumulator<graph_type, kcore_update, accumulator> {
private:
  int num_active;
  int links;
public:
  accumulator() : num_active(0), links(0) { }

  void operator()(icontext_type& context) {
   

    vertex_data & vdata = context.vertex_data();
    if (debug)
      logstream(LOG_INFO)<<"Entering node: " << context.vertex_id() << std::endl;

    if (!vdata.active)
      return;

    int increasing_links = 0;
    
    edge_list_type outedgeid = context.out_edges();
    //edge_list_type inedgeid = context.in_edges();

    for(size_t i = 0; i < outedgeid.size(); i++) {
      const edge_type& edge = reference_graph->find(context.vertex_id(), outedgeid[i].target());
      if (edge.offset() != (uint)-1)
          edge_count[edge.offset()]++;
      total_edges++;
    }
    /*for (size_t i =0; i < inedgeid.size(); i++){
      const vertex_data & other = context.const_vertex_data(inedgeid[i].source());
        if (other.active)
          cur_links++;
    }*/
    //links += increasing_links;
    //if (vdata.active)
    //  num_active++;
  };

  void operator+=(const accumulator& other) { 
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
}; // end of  accumulator





int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile;
  std::string format = "matrixmarket";
  int unittest = 0;
  int lineformat = MATRIX_MARKET_4;
  bool gzip = true;
  bool stats = false;
  std::string filter = "";
  int reference = 0;
  int max_graph = 1000;
  std::string list_dir = "/usr2/bickson/daily.sorted/";
  std::string dir_path = "/usr2/bickson/bin.graphs/";
  std::string out_dir = "/usr0/bickson/"; 
 
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
  clopts.attach_option("references", &reference, reference, "reference - why day to compare to?"); 
  clopts.attach_option("max_graph", &max_graph, max_graph, "maximum number of graphs parsed");
  clopts.attach_option("list_dir", &list_dir, list_dir, "directory with a list of file names to parse");
  clopts.attach_option("dir_path", &dir_path, dir_path, "actual directory where files are found");
  clopts.attach_option("outdir", &out_dir, out_dir, "output dir");
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  active_nodes_num = ivec(max_iter+1);
  active_links_num = ivec(max_iter+1);


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

  nodes = 121408373;
  matrix_info.rows = matrix_info.cols = nodes;
  core.set_scope_type("vertex");

    
    multigraph_type multigraph;
    multigraph.load(list_dir, dir_path, filter, true);
    matrix_info.nonzeros = core.graph().num_edges();


  if (stats){
    calc_multigraph_stats_and_exit<multigraph_type>(&multigraph, matrix_info);
  }
  logstream(LOG_INFO)<<"Going to load reference graph: " << reference << std::endl;
  graphlab::timer mytimer; mytimer.start();
  multigraph.doload(reference);
  reference_graph = multigraph.graph(0);
  edge_count = new uint[multigraph.graph(0)->num_edges()];

  int pass = 0;
    logstream(LOG_INFO)<<mytimer.current_time() << ") Going to run k-cores iteration " << iiter << std::endl;

      for (int i=0; i< std::min(multigraph.num_graphs(),max_graph); i++){
       if (i != reference){
       accumulator acum;
       multigraph.doload(i);
       core.graph() = *multigraph.graph(1);
       assert(multigraph.get_node_vdata()->size() == nodes);
       core.add_sync("sync", acum, 1000);
       core.add_global("NUM_ACTIVE", int(0));
       core.sync_now("sync");
       logstream(LOG_INFO)<<mytimer.current_time()<<") Finished giong over graph number " << i << std::endl;
       multigraph.unload(1);
       } 
     }

  
  std::cout << "KCORES finished in " << mytimer.current_time() << std::endl;
  std::cout << "Number of updates: " << pass*core.graph().num_vertices() << " pass: " << pass << std::endl;
  for (int i=0; i< 1000; i++)
     std::cout<<i<<": "<<edge_count[i]<<std::endl;



    write_output_vector_binary(out_dir + boost::lexical_cast<std::string>(reference) + "edge_count.bin", edge_count, reference_graph->num_edges());

    return EXIT_SUCCESS;
    uint * hist = histogram(edge_count, reference_graph->num_edges(), 29);
     std::ofstream out_file(std::string(out_dir +  ".hist.gz").c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening output file " << out_dir  << ".hist.gz" << std::endl;
    boost::iostreams::filtering_stream<boost::iostreams::output> fout;
    fout.push(boost::iostreams::gzip_compressor());
    fout.push(out_file);
    assert(fout.good()); 

     

    //for (int i=0; i< 29; i++)
    //  fout << hist[i] << std::endl;
   boost::unordered_map<uint, std::string> nodeid2hash;
   nodeid2hash.rehash(nodes);
   std::ifstream ifs((out_dir + ".reverse.map").c_str());
   {
   graphlab::iarchive ia(ifs);
   ia >> nodeid2hash;
   }
  
   for (int i=0; i< reference_graph->num_vertices(); i++){
      edge_list edges = reference_graph->out_edges(i);
      for (int j=0; j < edges.size(); j++){
        if (edge_count[edges[j].offset()] == 28)
          fout << nodeid2hash[edges[j].source()] << " " << nodeid2hash[edges[j].target()] << endl;     
      }      
   }
   fout.pop(); fout.pop();
   out_file.close();


  //multigraph.unload_all();
   return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
