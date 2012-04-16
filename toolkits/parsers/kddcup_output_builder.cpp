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
 *      http://graphlab.org
 *
 * Code by Danny Bickson, CMU
 *
 */


#include <cstdio>
#include <map>
#include <iostream>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;

string datafile;
string datafileB;
string prediction_file;
string prediction_fileB;
bool debug = false;

struct vertex_data {
  string filename; //file name with test data (either earlier or later)
  string prediction_file;  //file name with computed predictions for each test item
  vertex_data() { }
  vertex_data(string _filename, string _prediction_file) : filename(_filename), prediction_file(_prediction_file) { }
}; // end of vertex_data

struct edge_data {
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

struct kdd_output_builder :
   public graphlab::iupdate_functor<graph_type, kdd_output_builder>{
   void operator()(icontext_type& context) {
   
typedef graphlab::graph<vertex_data2, edge_data2>::edge_list_type edge_list;
   int nodes = context.get_global<int>("NUM_NODES");

    vertex_data& vdata = context.vertex_data();
    FILE * pfile = open_file(vdata.filename + ".csv", "w");
    if (vdata.filename == datafile)
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
    //go over each user in test data
    for (int i=0; i< nodes; i++){
      vec ratings = zeros(_graph.num_out_edges(i));
      edge_list edges  = _graph.out_edges(i);
      //find out the ratings computed by graphlab
      for (uint j=0; j< edges.size(); j++){
         ratings[j] = ret[cnt]; 
         cnt++;
      }
      if (edges.size() > 0 ){
        users++;
        //sort the rating
        ivec pos = sort_index(ratings);
        //print user id
        fprintf(pfile, "%d,", i+1);
        //print out the 3 items with the highest ratings (or only 2 if there are only two items)
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
    logstream(LOG_INFO)<<"Successfully created csvput file: " << vdata.filename << ".csv" << std::endl;
    } //operator()    
}; //update_functor




int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  int numnodes = 2421057; // 1420949264;//121408373;
  
  clopts.attach_option("dataA", &datafile, datafile,
                       "test A input file");
  clopts.add_positional("dataA");
  clopts.attach_option("prediction_fileA", &prediction_file, prediction_file, "prediction file name for testA");
  clopts.add_positional("prediction_fileA");
  clopts.attach_option("dataB", &datafileB, datafileB,
                       "test B input file");
  clopts.add_positional("dataB");
  clopts.attach_option("prediction_fileB", &prediction_fileB, prediction_fileB, "prediction file name for testB");
  clopts.add_positional("prediction_fileB");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab Parsers library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl;
 
  // Create a core
  graphlab::core<graph_type, kdd_output_builder> core;
  core.set_options(clopts); // Set the engine options

  vertex_data data(datafile, prediction_file);
  vertex_data dataB(datafileB, prediction_fileB);
  core.graph().add_vertex(0, data);
  core.graph().add_vertex(1, dataB);

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(kdd_output_builder());
 
  core.add_global("NUM_NODES", numnodes);
  double runtime= core.start();
  std::cout << "Finished in " << runtime << std::endl;

  logstream(LOG_INFO)<<"Merging the two files together"<<endl;
  int rc = system((string("cat ") + datafileB + ".csv >> " + datafile + ".csv").c_str());
  if (rc == -1){
     perror("failed cat");
     logstream(LOG_FATAL)<<"failed to concatanate the two test files"<<endl;
  }
     
  rc = remove((datafileB + ".csv").c_str());
  if (rc == -1){
     perror("failed delete");
     logstream(LOG_FATAL)<<"failed to remove temp file"<< datafileB << " .csv" << endl;
  }
 

  rc = system(("zip submission.zip " + datafile + ".csv").c_str());
  if (rc == -1){
     perror("failed zip");
     logstream(LOG_FATAL)<<"failed to zip submisison file"<< datafile << " .csv" << endl;
  }
  logstream(LOG_INFO)<<"Successfully created zip file: submission.zip" << endl;
  rc = remove((datafile + ".csv").c_str());
  if (rc == -1){
     perror("failed delete");
     logstream(LOG_FATAL)<<"failed to remove temp file"<< datafile << " .csv" << endl;
  }
  logstream(LOG_INFO)<<"removed temporary files" << endl;
   return EXIT_SUCCESS;
}



