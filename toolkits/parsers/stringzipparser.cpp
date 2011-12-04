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
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>





#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;


bool debug = false;

map<string,uint> hash2nodeid;
std::string datafile;
atomic<unsigned int> conseq_id;
graphlab::mutex mymutex;


struct vertex_data {
  string filename;
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

struct edge_data {
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data);



void assign_id(uint & outval, const string &name){
  mymutex.lock();
  outval = hash2nodeid[name];
  if (outval == 0){
      conseq_id.inc();
      hash2nodeid[name] = conseq_id.value;
      outval = conseq_id.value;
  }
  mymutex.unlock();
}


void find_ids(uint & from, uint & to, const string &buf1, const string& buf2){

   assign_id(from, buf1);
   assign_id(to, buf2);
   assert(from != to && from > 0 && to > 0);
}


/***
* Line format is: PnLaCsEnqei atslBvPNusB 050803 235959 590 
*/

struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {

    //open file
    vertex_data& vdata = context.vertex_data();
    std::ifstream in_file(vdata.filename.c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening input file: " << vdata.filename << std::endl;
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);  

    std::ofstream out_file(std::string(vdata.filename + ".out").c_str());

    char linebuf[256], buf1[256], buf2[256], buf3[256], buf4[256];
    char saveptr[1024];
    int duration;
    int line = 1;
    int lines = context.get_global<int>("LINES");
   
    while(true){
      fin.getline(linebuf, 128);
      char *pch = strtok_r(linebuf," \r\n\t",(char**)&saveptr);
      strncpy(buf1, pch, 18);
      pch = strtok_r(NULL, " \r\n\t;",(char**)&saveptr);
      strncpy(buf2, pch, 18);
      pch = strtok_r(NULL, " ",(char**)&saveptr);
      strncpy(buf3, pch, 6);
      buf3[6] = ' ';
      pch = strtok_r(NULL, " ",(char**)&saveptr);
      strncpy(buf3+7,pch,6);
      pch = strtok_r(NULL, " \r\n\t",(char**)&saveptr);
      duration = atoi(pch);
      unsigned long int timeret = datestr2uint64(std::string(buf3));
      uint from, to;
      find_ids(from, to, buf1, buf2);
      if (debug && line <= 10)
         cout<<"Read line: " << line << " From: " << from << " To: " << to << " time: " << timeret << " val: " << duration << endl;
 
      out_file << from << " " << to << " " << timeret << " " << duration << endl;
      fin.read(buf1,1); //go over \n
      line++;
      if (lines && line>=lines)
	 break;

      if (debug && (line % 50000 == 0))
        logstream(LOG_INFO) << "Parsed line: " << line << " map size is: " << hash2nodeid.size() << std::endl;
      if (fin.eof())
        break;
    } 

   logstream(LOG_INFO) <<"Finished parsing total of " << line << " lines in file " << vdata.filename << endl;

    // close file
    fin.pop(); fin.pop();
    in_file.close();
    out_file.close();
  }





  };


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

  std::string format = "plain";
  int unittest = 0;
  int lines = 0;
  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("lines", &lines, lines, "limit number of read lines to XX");
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

  //unit testing
  if (unittest == 1){
  }


  vertex_data vdata("HIDDEN.20050803.gz");
  core.graph().add_vertex(vdata);
  vertex_data vdata1("HIDDEN.20050810.gz");
  core.graph().add_vertex(vdata1);
  vertex_data vdata2("HIDDEN.20050817.gz");
  core.graph().add_vertex(vdata2);
  vertex_data vdata3("HIDDEN.20050824.gz");
  core.graph().add_vertex(vdata3);
     std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringzipparser_update());
 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
  core.add_global("LINES", lines); 
  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;

  //vec ret = fill_output(&core.graph(), matrix_info, JACOBI_X);

  //write_output_vector(datafile + "x.out", format, ret);


  if (unittest == 1){
  }

   return EXIT_SUCCESS;
}



