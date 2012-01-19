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
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;
using namespace std;


bool debug = false;
bool quick = false;
boost::unordered_map<string,uint> hash2nodeid;
std::string datafile;

graphlab::timer mytime;


int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string format = "plain";
  std::string dir = "/mnt/bigbrofs/usr0/bickson/phone_calls/";
  std::string outdir = "/usr2/bickson/filtered.hours/";
  std::string filter;
  int unittest = 0;
  int lines = 0;
  int nodes = 121408373;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("lines", &lines, lines, "limit number of read lines to XX");
  clopts.attach_option("quick", &quick, quick, "quick mode");
  clopts.attach_option("dir", &dir, dir, "path to files");
  clopts.attach_option("filter", &filter, filter, "select files starting with prefi [filter]");
  clopts.attach_option("outdir", &outdir, outdir, "output directory");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;

  mytime.start();
 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
    mytime.start();
    hash2nodeid.rehash(nodes);
    logstream(LOG_INFO)<<"Opening input file " << dir << datafile << ".map" << std::endl;
   std::ifstream ifs((dir + ".map").c_str());
   {
   graphlab::iarchive ia(ifs);
   ia >> hash2nodeid;
   }
   logstream(LOG_INFO)<<"Finished reading input file in " << mytime.current_time() << std::endl;
   

  boost::unordered_map<uint, std::string> nodeid2hash;
  boost::unordered_map<std::string, uint>::const_iterator i;
  for (i = hash2nodeid.begin(); i != hash2nodeid.end(); i++){ 
    nodeid2hash[i->second] = i->first;
  }
  std::cout << "Finished in " << mytime.current_time() << std::endl;
  assert(hash2nodeid.size() == nodeid2hash.size());
   std::ofstream ofs((outdir + ".reverse.map").c_str());
  graphlab::oarchive oa(ofs);
  {
  oa << nodeid2hash; 
  }
  logstream(LOG_INFO)<<"Finished saving reverse map in: " << mytime.current_time() << std::endl;
   return EXIT_SUCCESS;
}



