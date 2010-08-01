

#ifndef distmulticoemapp_pull_HPP
#define distmulticoemapp_pull_HPP

#include <cstdio>
#include <string>
#include <vector>
#include <map>

#include <graphlab.hpp>
#include <graphlab/distributed/graph/distributed_graph.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
#include <graphlab/distributed/distributed_engine.hpp>
#include <graphlab/distributed/distributed_scheduler_wrapper.hpp>
#include <graphlab/distributed/distributed_round_robin_scheduler.hpp>

using namespace graphlab;

typedef distributed_graph<blob, unsigned short> coem_distributed_graph;
typedef graph<blob, unsigned  short> coem_graph;

typedef types<coem_distributed_graph> gl_types;

#define SCHEDULER distributed_round_robin_scheduler

std::vector<double> l1_residual(blob_graph &g);



class distmulticoemapp_pull : public iapp<coem_graph> {
	
public:
  distmulticoemapp_pull(distributed_control * _dc, std::string npsfile, 
               std::string contextfile, 
               std::string matrixfile, 
               std::string seedsdir, 
               std::string negseedsdir);
  ~distmulticoemapp_pull();
  
  void start();
  
  bool prune_edges;

private:

 void prepare(distributed_shared_data<coem_distributed_graph>& sdm, 
                           distributed_engine<coem_distributed_graph,
      distributed_scheduler_wrapper<coem_distributed_graph, SCHEDULER<coem_distributed_graph> > > *  graphlab) ;
 
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
  
  bool on_the_fly_partition;

  coem_graph * g; 
  coem_distributed_graph distgraph;
  distributed_control*  dc; 
  std::vector<std::string> categories;	
  

};

#endif

