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
#include <graphlab/graph/graph2.hpp>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include <graphlab/macros_def.hpp>
using namespace graphlab;
using namespace std;

#define DUMMY 987654321
bool debug = false;
bool quick = false;
boost::unordered_map<std::string, bool> edges_in_28;
std::map<std::string, uint> edges_in_28_count;
boost::unordered_map<std::string, uint> hash2nodeid;
std::string datafile;
//atomic<unsigned int> conseq_id;
uint conseq_id;
graphlab::mutex mymutex;
graphlab::timer mytime;
unsigned long long total_lines = 0;
unsigned long long total_selected = 0;

struct vertex_data {
  string filename;
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

struct edge_data {
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data, int & dateret, int & timeret, int trhead_id);

struct vertex_data2 {
  bool active;
  int kcore, degree;

  vertex_data2() : active(true), kcore(-1), degree(0)  {}

  void add_self_edge(double value) { }

  void set_val(double value, int field_type) { 
  }
  //only one output for jacobi - solution x
  double get_output(int field_type){ 
    return -1; //TODO
  }
}; // end of vertex_data

struct edge_data2 {
  edge_data2()  { }
  //compatible with parser which have edge value (we don't need it)
  edge_data2(double val)  { }
};





/***
* Line format is: PnLaCsEnqei atslBvPNusB 050803 235959 590 
*/
/*
YVjAeZQjnVA IfrTTVlatui 050803 000000 156
GNgrmichxmG GNgriWokEhN 050803 000000 143
YnRdCKZkLao MHexzaXWCPL 050803 000000 0
RGNReqpKcZw RGNRSTDdqew 050803 000000 0
LPHSeuGhYkN ZFwbovKzAxY 050803 000000 1
sijmyRRfkwl XtqJaHYFEPqbZqNGPCr 050803 000000 68
*/
struct stringzipparser_update :
   public graphlab::iupdate_functor<graph_type, stringzipparser_update>{
   void operator()(icontext_type& context) {
    
   std::string dir = context.get_global<std::string>("PATH");
   std::string outdir = context.get_global<std::string>("OUTPATH");
   int mythreadid = thread::thread_id();

    //open file
    vertex_data& vdata = context.vertex_data();
    gzip_in_file fin((dir + vdata.filename));
    //gzip_out_file fout(outdir + vdata.filename + ".out.gz");

    edge_data2 edge; 

    char linebuf[256], buf1[256], buf2[256], buf3[256];
    char saveptr[1024];
    int duration, dateret, timeret;
    int line = 1;
    int lines = context.get_global<int>("LINES");
   
    while(true){
      fin.get_sp().getline(linebuf, 128);
      if (fin.get_sp().eof())
        break;

      if (linebuf[0] == '%') //skip comments
        continue;

      char *pch = strtok_r(linebuf," \r\n\t",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
        return;
       }
      strncpy(buf1, pch, 20);
      pch = strtok_r(NULL, " \r\n\t;",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
       strncpy(buf2, pch, 20);
       if (buf1[0] == '9' && buf2[0] == '9') //placeholder for matrix market size, to be done later
           continue;
       
      pch = strtok_r(NULL, " ",(char**)&saveptr);
      if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
         strncpy(buf3, pch, 6);
        buf3[6] = ' ';
        pch = strtok_r(NULL, " ",(char**)&saveptr);
        if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
         strncpy(buf3+7,pch,6);
        pch = strtok_r(NULL, " \r\n\t",(char**)&saveptr);
       if (!pch){
        logstream(LOG_ERROR) << "Error when parsing file: " << vdata.filename << ":" << line <<std::endl;
         return;
       }
        duration = atoi(pch);
        datestr2uint64(std::string(buf3), timeret, dateret, mythreadid);
   
        string hour = boost::lexical_cast<string>(timeret/3600);
        std::string srcid = std::string(buf1);
        std::string dstid = std::string(buf2);
        if (srcid.size() != 11 && srcid.size() != 19)
         logstream(LOG_WARNING)<<"Invalid string: " << srcid << " in line: " << line << endl;
       if (dstid.size() != 11 && dstid.size() != 19)
         logstream(LOG_WARNING) <<"Invalid string dst: " << dstid << " in line : " << line<< endl;

       if (edges_in_28.find(srcid + " " + dstid) != edges_in_28.end()){
         //fout.get_sp() << hash2nodeid[buf1] << " " << hash2nodeid[buf2] << " " << buf1 << " " << buf2 << " "
         //<< duration << " " << timeret << " " << dateret << endl;
         edges_in_28_count[srcid + " " + dstid + " " + hour]++;
         total_selected++;
      }

      line++;
      total_lines++;
      if (total_lines % 1000000 == 0)
        logstream(LOG_INFO) << mytime.current_time() << ") " << vdata.filename << " edges: " << total_lines << " seleted: " << total_selected << endl;

      if (lines && line>=lines)
	 break;

    } 

   logstream(LOG_INFO) <<mytime.current_time() << ") Finished parsing total of " << line << " lines in file " << vdata.filename << " total selected: " << total_selected << endl;

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
  std::string dir = "/mnt/bigbrofs/usr10/haijieg/edge_process/output/"; //"/usr2/bickson/daily.sorted/";
  std::string outdir = "/usr2/bickson/yahoo.graph/"; //"/usr2/bickson/bin.graphs/";
  int unittest = 0;
  int lines = 0;
  int numnodes =121408373; // 1420949264;//121408373;
  std::string filter = ""; //"day";

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
  clopts.attach_option("num_nodes", &numnodes, numnodes, "Number of nodes");
  clopts.attach_option("filter", & filter, filter, "Filter - parse files starting with prefix");
  clopts.attach_option("outdir", &outdir, outdir, "output directory");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  // Create a core
  graphlab::core<graph_type, stringzipparser_update> core;
  core.set_options(clopts); // Set the engine options
  core.set_scope_type("vertex");
  mytime.start();

  std::vector<std::string> in_files;
  if (datafile.size() > 0)
     in_files.push_back(datafile);
  else in_files = list_all_files_in_dir(dir, filter);
  assert(in_files.size() >= 1);
  for (int i=0; i< (int)in_files.size(); i++){
      vertex_data data(in_files[i]);
      core.graph().add_vertex(data);
  }

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(stringzipparser_update());
 
  //accumulator acum;
  //core.add_sync("sync", acum, sync_interval);
  core.add_global("LINES", lines); 
  core.add_global("PATH", dir);
  core.add_global("OUTPATH", outdir);
  core.add_global("NUM_NODES", numnodes);
  core.add_global("OUT_DIR", outdir);

  load_map_from_file(edges_in_28, outdir + ".28.edges");
  logstream(LOG_INFO) << "Loaded a total of " << edges_in_28.size() << " edges" << endl;
  load_map_from_file(hash2nodeid, "/mnt/bigbrofs/usr3/bickson/out_phone_calls/.map");
  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;
  gzip_out_file cnter(outdir + filter + ".28.edges.hour.gz");
  std::map<std::string, uint>::const_iterator it;
  for (it = edges_in_28_count.begin(); it != edges_in_28_count.end(); it++){
    cnter.get_sp() << it->first << " " << it->second << endl;
  }
  logstream(LOG_INFO) << "Found total unique entries: " << edges_in_28_count.size() << endl;
   return EXIT_SUCCESS;
}


#include <graphlab/macros_undef.hpp>

