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
typedef graph_type2::edge_type edge_type2;

/*
 * line format is:
 * 0	1	4298118681424644510	7686695	385	3	3	1601	5521	7709	576	490234
 * 0	1	4860571499428580850	21560664	37484	2	2	2255103	317	48989	44771	490234
 * 0	1	9704320783495875564	21748480	36759	3	3	4532751	60721	685038	29681	490234
 * 0	1	13677630321509009335	3517124	23778	3	1	1601	2155	1207	1422	490234
 * 0	1	3284760244799604489	20758093	34535	1	1	4532751	77819	266618	222223	490234
 *
*/
struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {
    
   vertex_data& vdata = context.vertex_data();
   gzip_in_file fin((vdata.filename), gzip);
   gzip_out_file fout(vdata.filename +".out", gzip);

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
      size_t val;
      int click = atoi(pch);
      ASSERT_GE(click, 0);
      pch = strtok_r(NULL, " \r\n\t:;",(char**)&saveptr);
      ASSERT_NE(pch, NULL);
      int impression = atoi(pch);
      ASSERT_GE(impression, 0);
      fout.get_sp() << ((double)click)/ impression << " | ";
      while(true){
         pch = strtok_r(NULL, " \r\n\t:;",(char**)&saveptr);
         if (pch == NULL)
             break;
         val = atoll(pch);
         items++;

         fout.get_sp() << val << " ";
         if (fin.get_sp().eof())
           break;
      }
      fout.get_sp() <<std::endl;
      total_lines++;

      if (debug && (total_lines % 1000000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << total_lines << endl;
  }

   logstream(LOG_INFO) <<"Finished parsing total of " << total_lines << " lines in file " << vdata.filename << endl <<
	                 " wrote total items (nnz) " << items << std::endl;

 }
};



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  clopts.attach_option("data", &datafile, datafile,
                       "feature input file");
  clopts.add_positional("data");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("gzip", &gzip, gzip, "Gzipped input file?");
  // Parse the command line arguments
  //
  
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
 
  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;

   return EXIT_SUCCESS;
}



