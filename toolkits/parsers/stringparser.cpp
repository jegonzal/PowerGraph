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
#include <map>
#include <iostream>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;


bool debug = false;

map<string,uint> hash2nodeid;
std::string datafile;
int conseq_id = 0;
struct vertex_data {
  string filename;
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

struct edge_data {
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data);



void find_ids(uint & from, uint & to, const string &buf1, const string& buf2){
   from = hash2nodeid[buf1];
   if (from == 0){
	hash2nodeid[buf1] = ++conseq_id;
        from = conseq_id;
   }
   to = hash2nodeid[buf2];
   if (to == 0){
       hash2nodeid[buf2] = ++conseq_id;
       to = conseq_id;
   }
        
   assert(from != to && from > 0 && to > 0);
}


/***
* Line format is: PnLaCsEnqei atslBvPNusB 050803 235959 590 
*/
struct stringparser_update :
  public graphlab::iupdate_functor<graph_type, stringparser_update> {
  void operator()(icontext_type& context) {
    /* GET current vertex data */
    vertex_data& vdata = context.vertex_data();
    FILE * f = open_file(vdata.filename.c_str(), "r");

    char buf1[256], buf2[256], buf3[256], buf4[256];
    int duration;
    int line = 1;

    while(true){
      int rc = fscanf(f, "%s %s %s %s %d\n", buf1, buf2, buf3, buf4,&duration);
      if (rc != 5){
        if (rc == EOF && errno == 0){
           logstream(LOG_INFO) << "Finishing reading " << (line-1) << " From file: " << vdata.filename << endl;
	   return;
        }
  	else logstream(LOG_FATAL)<< "Failed to parse line "<<line<<" in file: " << vdata.filename << endl;
      }
      char timebuf[256];
      sprintf(timebuf, "%s %s", buf3, buf4);
      unsigned long int timeret = datestr2uint64(std::string(timebuf));
      uint from, to;
      find_ids(from, to, buf1, buf2);
      if (debug && line <= 10)
         cout<<"Read line: " << line << " From: " << from << " To: " << to << " time: " << timeret << " val: " << duration << endl;
      line++;
    }
  }
}; // end of update_functor




/*
class accumulator :
  public graphlab::iaccumulator<graph_type, stringparser_update, accumulator> {
private:
  real_type real_norm, relative_norm;
public:
  accumulator() {}
  void operator()(icontext_type& context) {
  }
  void operator+=(const accumulator& other) { 
  }
  void finalize(iglobal_context_type& context) {
  }
}; // end of  accumulator
*/




int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string format = "matrixmarket";
  int unittest = 0;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
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
  graphlab::core<graph_type, stringparser_update> core;
  core.set_options(clopts); // Set the engine options

  //unit testing
  if (unittest == 1){
  }


  vertex_data vdata(datafile);
  core.graph().add_vertex(vdata);
  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringparser_update());

 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
 
  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;

  //vec ret = fill_output(&core.graph(), matrix_info, JACOBI_X);

  //write_output_vector(datafile + "x.out", format, ret);


  if (unittest == 1){
  }

   return EXIT_SUCCESS;
}



