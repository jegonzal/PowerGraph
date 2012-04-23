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
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;


bool debug = false;
bool quick = true;
bool gzip = false;
string user_data_file;
string item_data_file;
string training_data, validation_data;
string outfile;
int nodes = 2421057;
int split_training_time = 1320595199;
int MAX_FEATURE = 410;
int pos_offset = 0;
const int MATRIX_MARKET = 1;
const int VW = 2;
int output_format = 1;

struct vertex_data {
  vertex_data(){ 
  }
  void add_self_edge(double value) { }

  void set_val(double value, int field_type) { 
  }  
  double get_output(int field_type){ return NAN; }
}; 

struct edge_data2 {
  int value;
  int time;
  edge_data2(int value) : value(value) { time = 0; }
  edge_data2(int value, int time) : value(value), time(time) { }
  void set_field(int pos, int val){
     time = val;
  }
  int get_field(int pos){ if (pos == 0) return value; else if (pos == 1) return time; else assert(false); }
};

struct edge_data{
  double weight;
  edge_data(double weight): weight(weight) { }
  void set_field(int pos, double val){}
  double get_field(int pos) { return weight; }
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::graph<vertex_data, edge_data2> graph_type2;
typedef graph_type::edge_list_type edge_list;
typedef graph_type2::edge_list_type edge_list2;

int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Parsers  Library");

  clopts.attach_option("user_data", &user_data_file, user_data_file,
                       "user feature input file");
  clopts.add_positional("user_data");
  clopts.attach_option("item_data", &item_data_file, item_data_file, "item feature data file");
  clopts.add_positional("item_data");
  clopts.attach_option("training_data", &training_data, training_data, "training data");
  clopts.add_positional("training_data");
  clopts.attach_option("out_file", &outfile, outfile, "output file");
  clopts.add_positional("out_file");  
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("gzip", &gzip, gzip, "Gzipped input file?");
  clopts.attach_option("pos_offset", &pos_offset, pos_offset, "added offset to position values");
  clopts.attach_option("output_format", &output_format, output_format, "output format 1=Matrix market, 2=VW");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab parsers library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented parsers are: Call data records, document tokens "<< std::endl;

  timer mytimer; mytimer.start();
  graph_type2 training, validation;
  graph_type user_data, item_data;
  bipartite_graph_descriptor training_info, validation_info, user_data_info, item_data_info;
  training_info.force_non_square = true;
  validation_info.force_non_square = true; 
  load_matrixmarket_graph(training_data, training_info, 
                          training, MATRIX_MARKET_4, true);
  //we allow zero values of user/item features
  load_matrixmarket_graph(user_data_file, user_data_info, user_data, MATRIX_MARKET_3, true);
  load_matrixmarket_graph(item_data_file, item_data_info, item_data, MATRIX_MARKET_3, true);

  vec training_rating = zeros(training.num_edges());

  graph_type  out_graph;
  out_graph.resize(training.num_edges()+MAX_FEATURE);
  int training_instance = 0;
  int added_training = 0;
  int missing = 0;
  
  gzip_out_file fout(training_data + ".info", gzip);
  gzip_out_file fout2(training_data + ".data", gzip);
  MM_typecode out_typecode;
  mm_clear_typecode(&out_typecode);
  mm_set_real(&out_typecode); 
  mm_set_sparse(&out_typecode); 
  mm_set_matrix(&out_typecode);
  mm_write_cpp_banner(fout.get_sp(), out_typecode);
  
  for (int i=0; i< nodes; i++){
     //vertex_data & data = training.vertex_data(i);
     edge_list2 out_edges = training.out_edges(i);
     if (out_edges.size() == 0)
       continue;

     for (uint j=0; j < out_edges.size(); j++){
       training_instance++;
       ASSERT_LE(training_instance, training_rating.size());
       if (training_instance % 100000 == 0)
          logstream(LOG_INFO)<<"Handling training sample: " << training_instance<< endl;
       uint user = out_edges[j].source();
       uint item = out_edges[j].target() - nodes; 
       edge_data2 & edge = training.edge_data(out_edges[j]);
       int rating = edge.value;
       assert(rating == -1 || rating == 1 || rating == 0);
       edge_list user_features = user_data.out_edges(user);
       edge_list item_features = item_data.out_edges(item);
       if (!user_features.size() || !item_features.size()){
         missing++;
         continue;
       }
         
       training_rating[training_instance-1] = rating;
       if (output_format == VW)
         fout2.get_sp() << rating << " | ";
       
       for (uint k=0; k < user_features.size(); k++){
         assert(user_features[k].source() == user);
         int pos = user_features[k].target() - nodes;
         ASSERT_GE(pos , 0);
         ASSERT_LE(pos, MAX_FEATURE);
         edge_data edge2 = user_data.edge_data(user_features[k]);
         ASSERT_GE(edge2.weight , 0);
         //edge_data newedge(edge2.weight);
         if (output_format == MATRIX_MARKET)
             fout2.get_sp()<<training_instance+1<<" "<<pos+1<<" "<<edge2.weight<<endl;
         else 
             fout2.get_sp()<<pos+1<<" "<<edge2.weight<<" ";     
         
         //out_graph.add_edge(training_instance, training.num_edges() + pos, newedge);
         added_training++;
       }//for user features
       for (uint k=0; k < item_features.size(); k++){
         assert(item_features[k].source() == item);
         int pos = item_features[k].target() - nodes;
         ASSERT_GE(pos , 0);
         ASSERT_LE(pos, MAX_FEATURE);
         edge_data edge2 = item_data.edge_data(item_features[k]);
         ASSERT_GE(edge2.weight , 0);
         //edge_data newedge(edge2.weight);
         //out_graph.add_edge(training_instance, training.num_edges() + pos, newedge);
         if (output_format == VW)
           fout2.get_sp()<<pos+1<<" "<<edge2.weight<<" ";
         else
           fout2.get_sp()<<training_instance+1<<" "<<pos+1<<" "<<edge2.weight<<endl;
         added_training++;
        }//for item features
        if (output_format == VW)
          fout2.get_sp()<<std::endl;
      }//for out_edges
 } //for nodes
  mm_write_cpp_mtx_crd_size(fout.get_sp(), training_instance, MAX_FEATURE, added_training);
  //save_matrix_market_format_vector(training_data+".vec", training_rating, true, "%vector of ratings\n");
  if (output_format == MATRIX_MARKET)
    merge(training_data+".info", training_data+".data");
  std::cout << "Finished in " << mytimer.current_time() << " missing entries: " << missing << std::endl;
  return EXIT_SUCCESS;
} //main



