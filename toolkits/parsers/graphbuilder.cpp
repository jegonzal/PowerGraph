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
boost::unordered_map<string,uint> hash2nodeid;
std::string datafile;
//atomic<unsigned int> conseq_id;
uint conseq_id;
graphlab::mutex mymutex;
graphlab::timer mytime;
unsigned long long total_lines = 0;
unsigned long long self_edges = 0;

struct vertex_data {
  string filename;
  vertex_data(std::string _filename) : filename(_filename) { }
}; // end of vertex_data

struct edge_data {
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
unsigned long int datestr2uint64(const std::string & data, int & dateret, int & timeret);

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

typedef graphlab::graph2<vertex_data2, edge_data2> graph_type2;



void assign_id(uint & outval, const string &name){

  boost::unordered_map<string,uint>::iterator it = hash2nodeid.find(name);
  if (it != hash2nodeid.end()){
     outval = it->second;
     return;
  }
  mymutex.lock();
  outval = hash2nodeid[name];
  if (outval == 0){
      hash2nodeid[name] = ++conseq_id;
      outval = conseq_id;
  }
  mymutex.unlock();
}


void assign_id_quick(const string& name){
  mymutex.lock();
  hash2nodeid[name] = 1;
  mymutex.unlock();
}

void find_ids(uint & from, uint & to, const string &buf1, const string& buf2){

   assign_id(from, buf1);
   assign_id(to, buf2);
   //if (from == to)
   //   logstream(LOG_WARNING)<< " from equals to: " << from << " "  << buf1 << " " <<buf2 << std::endl;
   if (from == to)
     self_edges++;
   assert(from > 0 && to > 0);
}
/*
void save_to_bin(std::string filename, graph_type2 & _graph){
  
typedef graph2<vertex_data2,edge_data2>::edge_list_type edge_list_type2;   
typedef graph2<vertex_data2,edge_data2>::edge_type edge_type2;   
//typedef edge_type edge_type2;

   uint * nodes = new uint[_graph.num_vertices()+1];
   uint * innodes = new uint[_graph.num_vertices()+1];
   uint * edges = new uint[_graph.num_edges()];
   uint * inedges = new uint[_graph.num_edges()];
   nodes[0] = 0;
   innodes[0] = 0;
   int cnt = 0;
   int incnt = 0;
   for (int i=0; i< (int)_graph.num_vertices(); i++){
     nodes[i+1] = nodes[i]+ _graph.out_edges(i).size(); 
     assert(nodes[i+1] <= _graph.num_edges());
     innodes[i+1] = innodes[i] + _graph.in_edges(i).size();
     assert(innodes[i+1] <= _graph.num_edges());
     const edge_list_type2 out_edges = _graph.out_edges(i);
     const edge_list_type2 in_edges = _graph.in_edges(i);
     foreach(const edge_type2 & edge, out_edges){
       edges[cnt++] = (uint)edge.target();
       assert(edge.target() != i);
     } 
     foreach(const edge_type2 & edge, in_edges){
       inedges[incnt++] = edge.source();
       assert(edge.source() != i);
     }
   };
   assert(cnt == (int)_graph.num_edges());
   assert(incnt == cnt);
   write_output_vector_binary(filename + ".bin.nodes", nodes, _graph.num_vertices()+1); 
   write_output_vector_binary(filename + "-r.bin.nodes", innodes, _graph.num_vertices()+1); 
   write_output_vector_binary(filename + ".bin.edges", edges, _graph.num_edges());
   write_output_vector_binary(filename + "-r.bin.edges", inedges, _graph.num_edges());
};*/

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
   int nodes = context.get_global<int>("NUM_NODES");
   std::string outdir = context.get_global<std::string>("OUTPATH");

    //open file
    vertex_data& vdata = context.vertex_data();
    std::ifstream in_file((dir + vdata.filename).c_str(), std::ios::binary);
    logstream(LOG_INFO)<<"Opening input file: " << dir << vdata.filename << std::endl;
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;
    fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);  

    edge_data2 edge; 
    graph_type2 _graph;
    _graph.resize(nodes);

    char linebuf[256], buf1[256], buf2[256];
    char saveptr[1024];
    int line = 1;
    int lines = context.get_global<int>("LINES");
   
    while(true){
      fin.getline(linebuf, 128);
      if (fin.eof())
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
       uint from = boost::lexical_cast<uint>(buf1);
       uint to = boost::lexical_cast<uint>(buf2);
       if (from == DUMMY && to == DUMMY) //placeholder for matrix market size, to be done later
           continue;
       _graph.add_edge(from- 1, to - 1, edge);  //traslate edge offset to start from zero
       //_graph.add_edge(to, from, edge); 

      line++;
      total_lines++;
      if (total_lines % 1000000 == 0)
        logstream(LOG_INFO) << mytime.current_time() << ") " << vdata.filename << " edges: " << total_lines << endl;

      if (lines && line>=lines)
	 break;

    } 

   logstream(LOG_INFO) <<mytime.current_time() << ") Finished parsing total of " << line << " lines in file " << vdata.filename << endl;

    // close file
    fin.pop(); fin.pop();
    _graph.finalize();
    logstream(LOG_INFO) << mytime.current_time() << ") " << outdir + vdata.filename << " Going to save Graph to file" << endl;
    save_to_bin(outdir + vdata.filename, _graph);
    logstream(LOG_INFO) << mytime.current_time() << ") " << outdir + vdata.filename << " Finished saving Graph to file" << endl;
    in_file.close();
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
  int numnodes = 121408373;
  std::string filter = "x"; //"day";

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


  double runtime= core.start();
 
  std::cout << "Finished in " << runtime << std::endl;
  std::cout << "Total number of edges: " << self_edges << std::endl;
   return EXIT_SUCCESS;
}


#include <graphlab/macros_undef.hpp>

