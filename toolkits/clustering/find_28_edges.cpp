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
int nodes = 0;
uint * edge_count;
unsigned long long total_edges = 0;
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



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile;
  std::string format = "matrixmarket";
  int unittest = 0;
  int lineformat = MATRIX_MARKET_4;
  bool gzip = true;
  std::string filter = "";
  int reference = 0;
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
  clopts.attach_option("nodes", &nodes, nodes, "number of nodes"); 
  clopts.attach_option("gzip", &gzip, gzip, "gzipped input file?");
  clopts.attach_option("filter", & filter, filter, "Filter - parse files starting with prefix");
  clopts.attach_option("references", &reference, reference, "reference - why day to compare to?"); 
  clopts.attach_option("list_dir", &list_dir, list_dir, "directory with a list of file names to parse");
  clopts.attach_option("dir_path", &dir_path, dir_path, "actual directory where files are found");
  clopts.attach_option("outdir", &out_dir, out_dir, "output dir");
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }


  std::cout << "Load graph" << std::endl;
  nodes = 121408372;
    
  multigraph_type multigraph;
  multigraph.load(list_dir, dir_path, filter, true);


  logstream(LOG_INFO)<<"Going to load reference graph: " << reference << std::endl;
  graphlab::timer mytimer; mytimer.start();
  multigraph.doload(reference);
  reference_graph = multigraph.graph(0);
  edge_count = new uint[multigraph.graph(0)->num_edges()];

    uint * edge_count = read_input_vector_binary<uint>(out_dir + boost::lexical_cast<std::string>(reference) + "edge_count.bin", (int)reference_graph->num_edges());

    uint * hist = histogram(edge_count, reference_graph->num_edges(), 29);
   logstream(LOG_INFO)<<"Historgram of 28 edges is: " << hist[28] << endl;
    gzip_out_file(out_dir + ".hist.gz");
     
   boost::unordered_map<uint, std::string> nodeid2hash;
   nodeid2hash.rehash(nodes);
   load_map_from_file(nodeid2hash, out_dir + ".reverse.map");
   logstream(LOG_INFO)<<"Loaded a map of size: " << nodeid2hash.size() << endl;
   assert(nodeid2hash.size() == (uint)nodes);

   boost::unordered_map<std::string, bool> edges_in_28;
   int cnt =0; 
   int found = 0;
   for (int i=0; i< (int)reference_graph->num_vertices(); i++){
      edge_list edges = reference_graph->out_edges(i);
      for (int j=0; j < (int)edges.size(); j++){
        if (edge_count[edges[j].offset()] == 28){
          //fout << no
          //deid2hash[edges[j].source()] << " " << nodeid2hash[edges[j].target()] << endl;    
          std::string srcid = nodeid2hash[edges[j].source()+1];
          std::string dstid = nodeid2hash[edges[j].target()+1];
	  assert(srcid != dstid);
          if (srcid.size() != 11 && srcid.size() != 19)
            logstream(LOG_WARNING)<<" Found an src edge with too short hash: " << srcid << " " << srcid.size() << " " << " id: " << edges[j].source() << " " << dstid << endl;
          if (dstid.size() != 11 && dstid.size() != 19)
            logstream(LOG_WARNING)<<" Found an dst edge with too short hash: " << dstid << " " << dstid.size() << " " << " id: " << edges[j].target() << " " << srcid << endl;
           //assert(srcid.size() == 11 || srcid.size() == 19);
          assert(dstid.size() == 11 || dstid.size() == 19);
          edges_in_28[srcid + " " + dstid] = true;
          found++;
        }
        cnt++;
      }     
   }
   logstream(LOG_INFO) << "Found " << found << " edges with 28 days" << endl;
   assert(edges_in_28.size() == (uint)found);
   save_map_to_file(edges_in_28,out_dir + ".28.edges");

  //multigraph.unload_all();
   return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
