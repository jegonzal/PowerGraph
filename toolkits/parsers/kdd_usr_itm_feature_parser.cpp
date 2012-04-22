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
std::string datafile;
unsigned long long total_lines = 0;
int pos_offset = 0;
int nodes = 2421057;
int MAX_FEATURE = 410;

struct vertex_data {
  string filename;
  vertex_data() { }
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data
struct vertex_data2 {
  vertex_data2(){ 
  }
  void add_self_edge(double value) { }

  void set_val(double value, int field_type) { 
  }  
  double get_output(int field_type){ return NAN; }
}; 

struct edge_data2{
  double weight;
  edge_data2(double weight): weight(weight) { }
  void set_field(int pos, double val){}
  double get_field(int pos) { return weight; }
};

struct edge_data{};
typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::graph<vertex_data2, edge_data2> graph_type2;


/***
* Line format is:
100044	39:1	2:17.91759	4:15	6:9	8:0	10:3	12:97	14:0	16:0	18:0	20:10.98612	22:34	24:0.05714286	26:19	28:0	30:0	32:6
100054	0:25	40:1	2:19.4591	4:4	6:0	8:0	10:0	12:1	14:0	16:0	18:0	20:0	22:0	24:0	26:0	28:0	30:0	32:68
100065	0:23	39:1	2:40.60443	4:12	6:0	8:2	10:0	12:19	14:0	16:0	18:0	20:6.931472	22:15	24:0.0625	26:9	28:0	30:0	32:10
100080	0:26	39:1	2:34.65736	4:1	6:10	8:0	10:0	12:3	14:0	16:0	18:2	20:0	22:13	24:0	26:13	28:1	30:30	32:0
100086	0:26	39:1	2:48.67535	4:5	6:0	8:0	10:1	12:1	14:0	16:1	18:0	20:6.931472	22:26	24:0.03703704	26:24	28:0	30:0	32:11
100097	0:31	39:1	2:43.30733	4:1	6:0	8:0	10:0	12:0	14:0	16:0	18:0	20:0	22:28	24:0	26:28	28:0	30:0	32:27
*/
struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {
    
   vertex_data& vdata = context.vertex_data();
   gzip_in_file fin((vdata.filename), gzip);
   graph_type2 out_graph;
   out_graph.resize(nodes+MAX_FEATURE);

    char linebuf[24000];
    char saveptr[1024];
    uint line = 1;
    uint lines = context.get_global<uint>("LINES");
    int pos = 0;
    double val = 0;
    int items = 0; 
    int min_pos = 99999999, max_pos = 0;

    while(true){
      fin.get_sp().getline(linebuf, 24000);
      if (fin.get_sp().eof())
        break;

      char *pch = strtok_r(linebuf," \r\n\t:",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line << "[" << linebuf << "]" << std::endl;
        return;
       }
      int from = atoi(pch) - 1;
      if (from <= 10000){
         logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line << " document ID is zero or less: " << from << std::endl;
         return;
      }
      ASSERT_LT(from, nodes);
  
      while(true){
         pch = strtok_r(NULL, " \r\n\t:",(char**)&saveptr);
         if (pch == NULL)
             break;
         pos = atoi(pch);
         ASSERT_GE(pos, 0);
         ASSERT_LT(pos, MAX_FEATURE);
         min_pos = std::min(pos, min_pos);
         max_pos = std::max(pos, max_pos);
         pch = strtok_r(NULL, " \r\n\t:",(char**)&saveptr);
         ASSERT_NE(pch, NULL);
         val = atof(pch);
         ASSERT_GE(val, 0);
         edge_data2 edge(val);
         ASSERT_NE(from, pos+nodes);
         out_graph.add_edge(from, pos+nodes, edge);
         items++;
         if (fin.get_sp().eof() || pch == NULL)
           break;
      }

      line++;
      total_lines++;
      if (lines && line>=lines)
	 break;

      if (debug && (line % 50000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << line << endl;
  }

  bipartite_graph_descriptor out_info;
  out_info.rows = nodes; out_info.cols = MAX_FEATURE;
  out_info.nonzeros = out_graph.num_edges();
  out_info.force_non_square = true;
  out_graph.finalize();
  save_matrixmarket_graph(vdata.filename +".out", out_info, out_graph);
   logstream(LOG_INFO) <<"Finished parsing total of " << line << " lines in file " << vdata.filename << endl <<
	                 " wrote total items (nnz) " << items << " min pos: " << min_pos << " max pos: " << max_pos << endl;


 }
};



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string format = "plain";
  std::string dir = "/mnt/bigbrofs/usr3/bickson/phone_calls/";
  std::string outdir = "/mnt/bigbrofs/usr3/bickson/out_phone_calls/";
  // int unittest = 0;
  uint lines = 0;
  clopts.attach_option("data", &datafile, datafile,
                       "feature input file");
  clopts.add_positional("data");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("lines", &lines, lines, "limit number of read lines to XX");
  clopts.attach_option("dir", &dir, dir, "path to files");
  clopts.attach_option("outdir", &outdir, outdir, "output directory");
  clopts.attach_option("gzip", &gzip, gzip, "Gzipped input file?");
  clopts.attach_option("pos_offset", &pos_offset, pos_offset, "added offset to position values");

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

  // Create a core
  graphlab::core<graph_type, stringzipparser_update> core;
  core.set_options(clopts); // Set the engine options
  core.set_scope_type("none");

 vertex_data data(datafile);
 core.graph().add_vertex(vertex_id_type(0), datafile);

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringzipparser_update());
 
  core.add_global("LINES", lines); 
  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;

   return EXIT_SUCCESS;
}



