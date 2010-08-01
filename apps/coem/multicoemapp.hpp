

#ifndef COEMAPP_HPP
#define COEMAPP_HPP

#include <cstdio>
#include <string>
#include <vector>
#include <map>

#include <graphlab.hpp>
#include <graphlab/util/generics/blob.hpp>


using namespace graphlab;

typedef graphlab::graph<blob,blob> blob_graph;
std::vector<double> l1_residual(blob_graph &g);


class multicoemapp : public iapp< blob_graph > {
	
public:
  multicoemapp(std::string npsfile, 
               std::string contextfile, 
               std::string matrixfile, 
               std::string seedsdir, 
               std::string negseedsdir);
  ~multicoemapp();
  
  void start();
  
  bool prune_edges;

private:

  void prepare(thread_shared_data<blob_graph>& sdm, iengine<blob_graph> * graphlab);
  void load();
  void loadAndConvertToBin(std::string binfile);
  std::map<std::string, vertex_id_t> loadVertices(short typemask, FILE * f);
  void loadCategories();
  void setseeds(std::map<std::string, vertex_id_t>& nps_map);
  void loadEdges(FILE * fcont_to_nps, 
                 std::map<std::string, vertex_id_t>& nps_map,
                 std::map<std::string, vertex_id_t>& ctx_map);
					
  void output_results();

  
  virtual global_dumper dump_function();  
  virtual std::vector<std::string> dump_headers();  
  virtual int dump_frequency();
  void retrospective_dump();
  
  
  std::string npsfile;
  std::string contextfile; 
  std::string matrixfile; 
  std::string seedsdir;
  std::string negseedsdir;
  
  blob_graph * g;
	
  std::vector<std::string> categories;	
  

};

#endif

