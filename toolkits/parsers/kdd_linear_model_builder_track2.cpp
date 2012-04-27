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
bool test = false;
string query_data_file, purchase_data_file, description_data_file, title_data_file, user_data_file;
string training_data_file;
int nodes =26243606;
int MAX_FEATURE = 410;
int pos_offset = 0;
const int MATRIX_MARKET = 1;
const int VW = 2;
int output_format = 1;
int missing = 0;

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

void expand_features(graph_type & query_data, size_t val, 
                     size_t & added_training, 
                     boost::iostreams::filtering_stream<boost::iostreams::output>& fout){

 edge_list query_features = query_data.out_edges(val);
 //ASSERT_GT(query_features.size() , 0);
 if (query_features.size() == 0)
   missing++;
 for (uint k=0; k < query_features.size(); k++){
 ASSERT_EQ(query_features[k].source(),val);
 int pos = query_features[k].target() - nodes;
 ASSERT_GE(pos , 0);
 ASSERT_LE(pos, MAX_FEATURE);
 edge_data edge2 = query_data.edge_data(query_features[k]);
 ASSERT_GE(edge2.weight , 0);
             //edge_data newedge(edge2.weight);
 fout<<pos<<" ";     
 added_training++;
 }
}  

int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Parsers  Library");

  clopts.attach_option("query_data", &query_data_file, query_data_file,
                       "query feature input file");
  clopts.add_positional("query_data");
  clopts.attach_option("purchase_data", &purchase_data_file, purchase_data_file, "item feature data file");
  clopts.add_positional("purchase_data");
  clopts.attach_option("training_data", &training_data_file, training_data_file, "training data");
  clopts.add_positional("training_data");
  clopts.attach_option("description_data", &description_data_file, description_data_file, "user keyword file");
  clopts.add_positional("description_data");
  clopts.attach_option("title_data", &title_data_file, title_data_file, "title data file");
  clopts.add_positional("title_data");
  clopts.attach_option("user_data", &user_data_file, user_data_file, "user data file");
  clopts.add_positional("user_data");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("gzip", &gzip, gzip, "Gzipped input file?");
  clopts.attach_option("pos_offset", &pos_offset, pos_offset, "added offset to position values");
  clopts.attach_option("output_format", &output_format, output_format, "output format 1=Matrix market, 2=VW");
  clopts.attach_option("max_feature", &MAX_FEATURE, MAX_FEATURE, "max number of feature");
  clopts.attach_option("test", &test, test, "test ?");
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
  graph_type query_data, purchase_data, description_data, user_data, title_data;
  bipartite_graph_descriptor training_info, validation_info, query_data_info, purchase_data_info, description_data_info, user_data_info, title_data_info;
  //we allow zero values of user/item features
  if (query_data_file != "")
    load_matrixmarket_graph(query_data_file, query_data_info, query_data, MATRIX_MARKET_3, true);
  if (purchase_data_file != "")
    load_matrixmarket_graph(purchase_data_file, purchase_data_info, purchase_data, MATRIX_MARKET_3, true);
  if (description_data_file != "")
     load_matrixmarket_graph(description_data_file, description_data_info, description_data, MATRIX_MARKET_3, true);
  if (title_data_file != "")
     load_matrixmarket_graph(title_data_file, title_data_info, title_data, MATRIX_MARKET_3, true);
  if (user_data_file != "")
     load_matrixmarket_graph(user_data_file, user_data_info, user_data, MATRIX_MARKET_3, true);
  vec training_rating = zeros(training.num_edges());

  graph_type  out_graph;
  out_graph.resize(training.num_edges()+MAX_FEATURE);
  int training_instance = 0;
  size_t added_training = 0;
  int missing = 0;
  int keyword_items = 0;

   gzip_in_file fin(training_data_file, gzip);
   gzip_out_file fout(training_data_file +".out", gzip);
    char linebuf[24000];
    char saveptr[1024];
    int items = 0; 
    int total_lines = 0;

    while(true){
      fin.get_sp().getline(linebuf, 24000);
      if (fin.get_sp().eof())
        break;

      char *pch = strtok_r(linebuf," \r\n\t:;",(char**)&saveptr);
      ASSERT_NE(pch, NULL);
      size_t val = 0;
      if (!test){
      int click = atoi(pch);
      ASSERT_GE(click, 0);
      pch = strtok_r(NULL, " \r\n\t:;",(char**)&saveptr);
      ASSERT_NE(pch, NULL);
      int impression = atoi(pch);
      ASSERT_GE(impression, 0);
      fout.get_sp() << (((double)click)/ ((double)impression)) << " | ";
      }
      else fout.get_sp() << "0 | " << atoll(pch) << " ";

      int field = 3;
      while(true){
         pch = strtok_r(NULL, " \r\n\t:;",(char**)&saveptr);
         if (pch == NULL && field == 13)
           break;
         if (pch == NULL && test && field == 12)
           break;
         val = atoll(pch);
         items++;
         if (field < 8)
           fout.get_sp() << val << " " ;
         else if (field == 8)
           expand_features(query_data, val, added_training, fout.get_sp());
	 else if (field == 9)
           expand_features(purchase_data, val, added_training, fout.get_sp());
         else if (field == 10)
           expand_features(title_data, val, added_training, fout.get_sp());
         else if (field == 11)
           expand_features(description_data, val,added_training, fout.get_sp());
         else if (field == 12)
           expand_features(user_data, val, added_training, fout.get_sp());
           
         field++;
         if (fin.get_sp().eof())
           break;

      }
      fout.get_sp() <<std::endl;
      total_lines++;

      if (debug && (total_lines % 1000000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << total_lines << endl;
  }

   logstream(LOG_INFO) <<"Finished parsing total of " << total_lines << " lines in file " << training_data_file << endl <<
	                 " wrote total items (nnz) " << items << std::endl << " missing: " << missing << endl;

  return EXIT_SUCCESS;
} //main



