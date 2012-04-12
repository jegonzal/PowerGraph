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


#include <cstdio>
#include <map>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/unordered_map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/unordered_map.hpp>
#define AVOID_PARALLEL_SORT 1
#define GRAPH2_NO_EDGE
#include <graphlab/graph/graph3.hpp>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include <graphlab/macros_def.hpp>
using namespace graphlab;
using namespace std;

bool debug = false;
bool quick = false;
boost::unordered_map<string,uint> hash2nodeid;
string datafile;
string prediction_file;
uint conseq_id;
graphlab::mutex mymutex;
graphlab::timer mytime;
unsigned long long total_lines = 0;
unsigned long long self_edges = 0;

struct vertex_data {
  string filename;
  string prediction_file; 
  vertex_data() { }
  vertex_data(string _filename, string _prediction_file) : filename(_filename), prediction_file(_prediction_file) { }
}; // end of vertex_data

struct edge_data {
  double val;
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data, int & dateret, int & timeret);

struct vertex_data2 {
  double A_ii;
  double value;
  //real_type y, Aii;
  //real_type pvec[JACOBI_X], pvec[JACOBI_REAL_X], pvec[JACOBI_PREV_X];
  vertex_data2(): A_ii(1) { //: y(0), Aii(1), pvec[JACOBI_X](0), pvec[JACOBI_REAL_X](0), 
                 // pvec[JACOBI_PREV_X](-1) 
  }
  void add_self_edge(double value) { A_ii = value; }

  void set_val(double value, int field_type) { 
  }  
  double get_output(int field_type){ return -1; }
}; // end of vertex_data


struct edge_data2 {
  double weight;
  edge_data2(double weight = 0) : weight(weight) { }
};


typedef graphlab::graph<vertex_data2, edge_data2> graph_type2;

struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {
   
typedef graphlab::graph<vertex_data2, edge_data2>::edge_list_type edge_list;
   int nodes = context.get_global<int>("NUM_NODES");

    vertex_data& vdata = context.vertex_data();
    FILE * pfile = open_file(vdata.filename + ".out", "w");
    fprintf(pfile, "id,clicks\n"); 
    graph_type2 _graph; 
    bipartite_graph_descriptor info;
    info.force_non_square = true;
    load_graph(vdata.filename, "matrixmarket", info, _graph, MATRIX_MARKET_4, 
       true, false, false);
    _graph.finalize();
   
    vec ret;
    bipartite_graph_descriptor vecinfo;
    load_matrix_market_vector(vdata.prediction_file, vecinfo, ret, false, true);
    ASSERT_EQ(ret.size(), _graph.num_edges()); 
 
    int cnt = 0;
    int users = 0;
    for (int i=0; i< nodes; i++){
      vec ratings = zeros(_graph.num_out_edges(i));
      edge_list edges  = _graph.out_edges(i);
      for (uint j=0; j< edges.size(); j++){
         ratings[j] = ret[cnt]; 
         cnt++;
      }
      if (edges.size() > 0 ){
        users++;
        ivec pos = sort_index(ratings);
        fprintf(pfile, "%d,", i+1);
        for (int j=0; j< std::min((int)ratings.size(),3); j++){
          fprintf(pfile, "%d ", edges[pos[ratings.size() - j - 1]].target()-nodes+1);
        }
        fprintf(pfile, "\n");
        if (i ==100019){ //sanity check to verify this user has the right item (there are only 2 items for her)
           assert((edges[pos[0]].target() -nodes+1)==1606574 || (edges[pos[1]].target() -nodes+1)==1606574);
           assert((edges[pos[0]].target() -nodes+1)==1774684 || (edges[pos[1]].target() -nodes+1)==1774684);
        }
      }
     } //for
     fclose(pfile);
    ASSERT_EQ(cnt, ret.size());
    logstream(LOG_INFO)<<"Finished going over " << cnt << " rating " << " for " << users << " unique users " << std::endl; 
    } //operator()    
}; //update_functor




int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  int numnodes = 2421057; // 1420949264;//121408373;
  
  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("prediction_file", &prediction_file, prediction_file, "prediction file name");
  clopts.add_positional("prediction_file");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("num_nodes", &numnodes, numnodes, "Number of nodes");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

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
  graphlab::core<graph_type, stringzipparser_update> core;
  core.set_options(clopts); // Set the engine options
  core.set_scope_type("vertex");
  mytime.start();

  vertex_data data(datafile, prediction_file);
  core.graph().add_vertex(0, data);

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringzipparser_update());
 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
  core.add_global("NUM_NODES", numnodes);


  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;
   return EXIT_SUCCESS;
}


#include <graphlab/macros_undef.hpp>

