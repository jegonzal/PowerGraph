/*
 *  \author akyrola
 *  Multi ctageory CoEM implementation, based on material by Justin Betteridge
 *
 */
 
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <dirent.h>

#include <boost/random.hpp>


#include <graphlab.hpp>

 
#include "distmulticoemapp.hpp"

#include <graphlab/distributed/pushy_distributed_engine.hpp>
#include <graphlab/distributed/distributed_scheduler_wrapper.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
#include <graphlab/distributed/distributed_round_robin_scheduler.hpp>
 

#include <graphlab/macros_def.hpp>

using namespace std;
using namespace graphlab;

// Last bit of flag determines type
#define TYPEMASK  (1<<2)
#define NP_MASK  (0) 
#define CTX_MASK  (1<<2)
#define SEEDMASK  0x1

// Hack to define seeds
float NEGSEEDPROB = -1e-30f ;
float POSSEEDPROB = 1.000001f;

float INITIAL_NP_VALUE = 0.0f;

#define TEXT_LENGTH 64

const size_t PARAM_M  = 1;
const size_t NUM_CATS = 2;
const size_t NUM_NPS  = 3;
const size_t NUM_CTX  = 4;

bool ROUNDROBIN = false;

// Hack for speeding up
std::vector<vertex_id_t> ** vertex_lookups;

std::vector<blob> * XSTAR;

double subsamplingRatio = 0.25;

int myprocid, numofprocs;

// lfidc is a modified weight formula for Co-EM (see Justin
// Betteridge's "CoEM results" page)
#define TFIDF(coocc, num_neighbors, vtype_total) (log(1+coocc)*log(vtype_total*1.0/num_neighbors))   

struct vertex_data {
  unsigned short flags;
  int nbcount;
  char text[TEXT_LENGTH];
  float normalizer;
  float p[0];
};


float TARGET_PRECISION = 1e-5;

  
/* Used with analyzer_listener */
void writeDump(blob_graph *g, 
               bool is_xstar, long int budget, 
               double * l1res, double * linf);
  
// Temporary hack 
//distmulticoemapp * coapp;


/// ====== UPDATE FUNCTION ======= ///
void coem_update_function(gl_types::iscope& scope, 
                          gl_types::icallback& scheduler,
                          gl_types::ishared_data* shared_data) {
  // scheduler.disable_buffering();
 
 
  /* Hacky optimization */
  std::vector<vertex_id_t> vlookup = *vertex_lookups[thread::thread_id()];
  vlookup.clear();
   
  /* Get shared vars */
  float param_m = shared_data->get_constant(PARAM_M).as<float>();
  size_t num_nps = shared_data->get_constant(NUM_NPS).as<size_t>();
  size_t num_ctxs = shared_data->get_constant(NUM_CTX).as<size_t>();
  unsigned int num_cats = shared_data->get_constant(NUM_CATS).as<int>();
   
  /* Get vertex data */
  vertex_data * vdata =
    scope.vertex_data().as_ptr<vertex_data>();
   
  /* We use directed edges as indirected, so either this vertex has only incoming  
     or outgoing vertices */
  bool use_outgoing = (scope.in_edge_ids().size()==0);
  std::vector<edge_id_t> edge_ids = 
    (use_outgoing ? scope.out_edge_ids() : scope.in_edge_ids());
   
  if (edge_ids.size() == 0) {
    //printf("Warning : dangling node %d %s\n", scope.vertex(), vdata->text);
    return;
  }	
   
  bool is_np = (vdata->flags & TYPEMASK) == NP_MASK;
  bool was_first_run = false;
   
  unsigned int vtype_total = (is_np ? num_nps : num_ctxs);
  unsigned int vtype_other_total = (is_np ? num_ctxs : num_nps); 
    
  /* First check is my normalizer precomputed. If not, do it. */
  if (vdata->normalizer == 0.0) {
     float norm = param_m * vtype_other_total;
     foreach(edge_id_t eid, edge_ids) {
	  unsigned short cooccurence_count = scope.edge_data(eid);
      const vertex_data * nb_vdata = 
        scope.const_neighbor_vertex_data(use_outgoing ? 
                                   scope.target(eid) :
                                   scope.source(eid)).as_ptr<vertex_data>();
      ASSERT_GT(nb_vdata->nbcount, 0);
      norm += TFIDF(cooccurence_count, nb_vdata->nbcount, vtype_total);
      
      if (scope.vertex() == 100000) {
      	 printf("100000:  edge %d  coo %d nbc %d total %d param_m %lf other_type %d \n", eid, cooccurence_count, nb_vdata->nbcount, vtype_total,
      	 				param_m, vtype_other_total);
      }
    }
    vdata->normalizer = norm;
    was_first_run = true; // Used to cause update to all neighbors
                          // regardless of residual
    if (scope.vertex()%20000 == 0)	
      printf("%d Computed normalizer: %lf \n", scope.vertex(), norm);
   }
	
  /***** COMPUTE NEW VALUE *****/
  float tmp[256];  // assume max 256 cats
  assert(num_cats<256);
  for(int i=0; i<256; i++) tmp[i] = param_m;
	
	
  foreach(edge_id_t eid, edge_ids) {
    unsigned short cooccurence_count = scope.edge_data(eid);
    vertex_id_t nbvid =use_outgoing ? scope.target(eid) : scope.source(eid);
    const vertex_data * nb_vdata = 
      scope.const_neighbor_vertex_data(nbvid).as_ptr<vertex_data>();
    for(unsigned int cat_id=0; cat_id<num_cats; cat_id++) {
      tmp[cat_id] += nb_vdata->p[cat_id] * 
        TFIDF(cooccurence_count, nb_vdata->nbcount, vtype_total);
    }
    // Ensure data is pushed properly
   
   
    vlookup.push_back(nbvid); 
  }
  // Normalize and write data
  float residual = 0.0;
  for(unsigned int cat_id=0; cat_id<num_cats; cat_id++) {
    tmp[cat_id] /= vdata->normalizer;
		 
    /* I am seed for this category? */
    if (!((vdata->p[cat_id] == POSSEEDPROB) || 
          vdata->p[cat_id] == NEGSEEDPROB)) {
      residual += std::fabs(tmp[cat_id] - vdata->p[cat_id]);
      vdata->p[cat_id] = tmp[cat_id];
			 
			 
    }  else {
      // I am seed - do not change
      if (vdata->p[cat_id] == POSSEEDPROB) {
        assert((vdata->flags & SEEDMASK) != 0); 
      }
    }
  }
	
  if (scope.vertex()%20000 == 0) 
    printf("%d %d Entering foreach: %lf \n", myprocid, scope.vertex(), residual);
    
  if (!ROUNDROBIN) {
   // Write data to edges and schedule if threshold reached
   double randomNum = graphlab::random::rand01();
   assert(randomNum >= 0.0);
   assert(randomNum <= 1.0);
 
    int sz = edge_ids.size();
    for(int l = 0; l<sz; l++) {
      vertex_id_t nbvid = vlookup[l];
      gl_types::update_task task(nbvid, coem_update_function);
   
         edge_id_t eid = edge_ids[l];
         unsigned short co_occurencecount = scope.edge_data(eid);
         const vertex_data * nb_vdata = (vertex_data *) 
           scope.neighbor_vertex_data(nbvid).as_ptr<vertex_data>();
         float neighbor_residual = 
           (nb_vdata->normalizer == 0.0f ? 1.0 : 
            residual * 
            TFIDF(co_occurencecount, 
                  vdata->nbcount, vtype_other_total) / 
            nb_vdata->normalizer);	
         if ((neighbor_residual/TARGET_PRECISION >= randomNum) || 
             was_first_run) {
           scheduler.add_task(task, neighbor_residual);
           // if (scope.vertex()%20000 == 0) printf("Added task %d =>
           // %d %f\n", scope.vertex(), nbvid, neighbor_residual);
         }
      
    }
  }
}

distmulticoemapp::distmulticoemapp(distributed_control * _dc, std::string npsfile, std::string contextfile, 
                           std::string matrixfile, std::string seedsdir, 
                           std::string negseedsdir) : distgraph(*_dc){
  dc = _dc;
  this->npsfile = npsfile;
  this->contextfile = contextfile;
  this->matrixfile = matrixfile;
  this->seedsdir = seedsdir;
  this->negseedsdir = negseedsdir;
  prune_edges = false;
  on_the_fly_partition = false;
  
  myprocid = dc->procid();
  numofprocs = dc->numprocs();
  
  std::cout << "Distributed multicoem starting. Procid = " << myprocid << "/" << numofprocs << std::endl;
}

distmulticoemapp::~distmulticoemapp() {
  delete(g);
}


int ncats = 0;
int dumpcount = 0;


void distmulticoemapp::start() { 
  //ncats = 100;
  //dumpcount=3;
  // retrospective_dump();
  // return;

  /* Load XSTAR */
  if (fopen("MC_XSTAR.dat", "r") != NULL) {
    XSTAR = new vector<blob>();
    std::cout << "Loading x_star..." << std::endl;
    std::ifstream fin("MC_XSTAR.dat");
    iarchive iarc(fin);
    iarc >> *XSTAR;
    fin.close();
   }
 /* subsampling ratio */
  if (opts.extra.length() > 0) {
      printf("SUBSAMPLING RATIO: %lf\n", subsamplingRatio);
     sscanf(opts.extra.c_str(), "%lf", &subsamplingRatio);
     printf("SUBSAMPLING RATIO: %lf\n", subsamplingRatio);
  }	

  /**** GRAPH LOADING ****/


  // Load categories
  loadCategories();
  
  dc->barrier();
  
  char ss[8]; 
  sprintf(ss, "%1.2f", subsamplingRatio);
  std::string basefile = npsfile + std::string(ss);
  
  char ss2[10];
  sprintf(ss2, "%d", myprocid);
  std::string procidstr = std::string(ss2);
  
  distgraph.set_constant_edges(true);
  distgraph.set_local_edges(true);

  
  if (opts.engine == "prepartition") {
  	   // Create partitions
	   load();

  	   ASSERT_EQ(dc->numprocs(), 1);
  	   printf("Number of edges: %ld\n", g->num_edges());
  	   printf("Partition method: %s\n", opts.scheduler.c_str());
  	   
  	   if (opts.scheduler == "metis") {
      	   coem_distributed_graph::partition_graph_tofile(*g, opts.ncpus, partition_method::PARTITION_METIS, basefile);
  	   } else if (opts.scheduler == "bfs") {
  	       coem_distributed_graph::partition_graph_tofile(*g, opts.ncpus, partition_method::PARTITION_BFS, basefile);
  	   } else if (opts.scheduler == "random") {
  	       coem_distributed_graph::partition_graph_tofile(*g, opts.ncpus, partition_method::PARTITION_RANDOM, basefile);
  	   } else ASSERT_MSG(false, "You need to set partition method as --scheduler option.");
  	   
  	    exit(0);
  } else if (opts.engine == "onthefly") {
     std::cout << "========== ON THE FLY PARTITION ========== " << std::endl;
     on_the_fly_partition = true;
     load();
  } else {
  	   // Load partitions
  	   distgraph.load(basefile, *dc);
  }
  
  dc->barrier();
  
  /**** GRAPHLAB INITIALIZATION *****/

  // Create engine
  
  
  pushy_distributed_engine<coem_distributed_graph,
      distributed_scheduler_wrapper<coem_distributed_graph, SCHEDULER<coem_distributed_graph> >,
      		  general_scope_factory<coem_distributed_graph> >
      		graphlab(*dc, distgraph, opts.ncpus);
  graphlab.set_default_scope(scope_range::VERTEX_CONSISTENCY);
  
  printf(" =================== \n ");
  
  dc->barrier();
  
  ROUNDROBIN = true;
  // Run 5 iterations
  

  /** Hack **/
  vertex_lookups = (std::vector<vertex_id_t> **) 
    malloc(sizeof(std::vector<vertex_id_t> *) * opts.ncpus);
  for(unsigned int i=0; i<opts.ncpus; i++) {
    vertex_lookups[i] = new std::vector<vertex_id_t>();
    vertex_lookups[i]->reserve(1e4);
  }
  

// Set register shared data
  distributed_shared_data<coem_distributed_graph> sdm(*dc);
  dc->barrier();
  if (myprocid == 0) {
	  float m = 0.01;
  	  sdm.set_constant(PARAM_M, any(m));
  }
  
  dc->barrier();
  std::cout << "going to prepare()" << std::endl;
  prepare(sdm, &graphlab);
  graphlab.set_shared_data_manager(&sdm);
  dc->barrier();
  timer t2;
  t2.start();
  
  /* Round robin barriers etc */
  if (ROUNDROBIN) {
      vertex_id_t zero = 0;
      graphlab.get_scheduler().add_task_to_all(coem_update_function, 1.0);
      graphlab.get_scheduler().set_option(scheduler_options::MAX_ITERATIONS, (void*)3);
      graphlab.get_scheduler().set_option(scheduler_options::DISTRIBUTED_CONTROL, (void*)dc);
      graphlab.get_scheduler().set_option(scheduler_options::BARRIER, &zero);  
      size_t lastNP =  sdm.get_constant(NUM_NPS).as<size_t>();
      graphlab.get_scheduler().set_option(scheduler_options::BARRIER, &lastNP);  
  }
  
  
  
  TARGET_PRECISION = opts.threshold; 
  graphlab.start();
  
  distgraph.send_vertices_to_proczero();
  dc->barrier();
  
  g = &distgraph.mgraph;

  std::cout << "Finished in " << t2.current_time() << " seconds." << std::endl;

  if (opts.visualizer == "analyzer") retrospective_dump();

  std::cout << "Going to output results... " << std::endl;
  
  if (myprocid == 0) output_results();
  /*
  // Write analysis data
  char analyzer_filename[255];
  sprintf(analyzer_filename, "coem_%d_%d.tsv", myprocid, numofprocs);
  FILE * af = fopen(analyzer_filename, "w");
  
  fprintf(af, "vertex,updates,neighbors,remote_neighbors\n");
  
  foreach(vertex_id_t vid, distgraph.my_vertices()) {
     int num_neighbors = 0;
     int num_remote = 0;
     foreach(edge_id_t ev, distgraph.in_edge_ids(vid)) {
        num_neighbors++;
        num_remote += (distgraph.owner(distgraph.source(ev)) != myprocid);
     }
     
     foreach(edge_id_t ev, distgraph.out_edge_ids(vid)) {
        num_neighbors++;
        num_remote += (distgraph.owner(distgraph.target(ev)) != myprocid);
     }
    
     fprintf(af, "%ld,%ld,%ld,%ld\n", (long int) vid, (long int) updatecounts[vid],(long int) num_neighbors, (long int) num_remote);
  }
  
  fclose(af);*/
 }


// Bit dirty... who cares (maybe Joey :) ). Should use bind()...
int sort_cat_id = 0;
bool cmp( vertex_data * a, vertex_data * b ) {
  return a->p[sort_cat_id] > b->p[sort_cat_id]; // Descending order
}

void distmulticoemapp::output_results() {
  int cat_id = 0;
  int n = g->num_vertices();

  foreach(std::string catname, categories) {
    std::string filename = catname + "_results.txt";
    FILE * f = fopen(filename.c_str(), "w");
		
    std::vector<vertex_data *> v = std::vector<vertex_data *>(n);
    int k = 0;
    for(int vid=0; vid<n; vid++) {
      vertex_data * vdata = g->vertex_data(vid).as_ptr<vertex_data>();
      if ((vdata->flags & TYPEMASK) == NP_MASK &&
          (vdata->p[cat_id] != POSSEEDPROB)) {
        v[k++] = g->vertex_data(vid).as_ptr<vertex_data>();
      }
    }
    v.resize(k);
    sort_cat_id = cat_id;
    std::cout << "Sorting for category:: " << catname << std::endl;
    sort(v.begin(), v.end(), cmp);
    // Output top
    for(int i=0; i<50; i++) {
      vertex_data * vdata = v[i];
      fprintf(f, "%d: %s  %1.5f\n", i+1, vdata->text, (float) vdata->p[cat_id]); 
    }
		
    fclose(f);
		 
    cat_id++;
  }
}

void distmulticoemapp::load() {
  
	
  // Check if already converted to binary
  std::string binfile = npsfile + "_graphlabbin.dat";
	
  if (subsamplingRatio != 1.0) {
    char s[255];
    sprintf(s, "%s_%1.2f", binfile.c_str(), subsamplingRatio);
    binfile = std::string(s);
  }
  std::cout << "Binfile: " << binfile << std::endl;
	
  FILE * f = fopen(binfile.c_str(), "r");
  if (f == NULL || on_the_fly_partition) {
    std::cout << "Going to convert to => " << binfile << std::endl;
    loadAndConvertToBin(binfile);
  } else {
    std::cout << "Loading from binary " << binfile << std::endl;
    g = new coem_graph();
    g->load(binfile);
  }       
}

void distmulticoemapp::prepare(distributed_shared_data<coem_distributed_graph>& sdm, 
                            pushy_distributed_engine<coem_distributed_graph,
      distributed_scheduler_wrapper<coem_distributed_graph, SCHEDULER<coem_distributed_graph> >,
      		  general_scope_factory<coem_distributed_graph> >   *  graphlab) {
  // Set register shared data
  /* M is used for "smoothing" */
  
  /*** Count NPs and CONTEXTs, add SEEDS as initial tasks */
  size_t n = distgraph.num_vertices();
  
  std::cout << "Seed mask " << SEEDMASK << std::endl;
  
  size_t num_contexts = 0;
  size_t num_nps = 0;
  size_t negseeds = 0, posseeds = 0;
  int ntasks = 0;
  for(size_t i = 0; i<n; i++) {
    vertex_data * vdata = distgraph.vertex_data(i).as_ptr<vertex_data>();
    if ((vdata->flags & TYPEMASK) == NP_MASK) {
          num_nps++;
    } else {
          num_contexts++;
    }
  
    if (!(distgraph.owner(i) == myprocid)) continue;

    /* Is seed? */
    if ((vdata->flags & SEEDMASK) != 0) {
      for(unsigned int cat_id=0; cat_id<categories.size(); cat_id++) {
        if (vdata->p[cat_id] == POSSEEDPROB) {
          // 	graphlab->get_scheduler().add_task(update_task(i, coem_update_function), 1.0); 
          //	ntasks++;
          posseeds++;
        }
      }
  		 
    }
  	  
    /* Count classes */
    if ((vdata->flags & TYPEMASK) == NP_MASK) {
      if (i%10000 == 0) std::cout << "NP: " << std::string(vdata->text) << std::endl;
  	  	
      /* Set initial values */
      for(unsigned int cat_id=0; cat_id<categories.size(); cat_id++) {
        if (vdata->p[cat_id] == 0.0f) {
          vdata->p[cat_id] = INITIAL_NP_VALUE;
        } 
        else if (vdata->p[cat_id] == NEGSEEDPROB) {
          //	graphlab->get_scheduler().add_task(update_task(i, coem_update_function), 1.0); 
          //	ntasks++;
          negseeds++;
        }
      }
  	  	
    } else  if (((vdata->flags & TYPEMASK) == CTX_MASK)) {
      if (i%10000 == 0) std::cout << "CTX: " << std::string(vdata->text) << std::endl;
      // Start with contexts
      if (!ROUNDROBIN)  graphlab->get_scheduler().
        add_task(update_task<coem_distributed_graph>(i, coem_update_function), 1.0); 
      ntasks++;
    } else assert(false);
  	  
    /* Count neighbors */
    vdata->nbcount = std::max(distgraph.in_edge_ids(i).size(), distgraph.out_edge_ids(i).size());
    if (i%10000 == 0 || i == 141400) std::cout << "Num neighbors: " << i << " : " << vdata->nbcount << std::endl;
	
	// Broadcast update to others
	distgraph.update_vertex(i);
  }
  std::cout << "Positive seeds: " << posseeds << std::endl;
  std::cout << "Negative seeds: " << negseeds << std::endl;
  std::cout << "Num of initial tasks: " << ntasks << std::endl;
  std::cout << "Contexts: " << num_contexts << " nps: " << num_nps << std::endl;
  
  if (myprocid == 0) {
 	 sdm.set_constant(NUM_CTX, any(num_contexts));
  	sdm.set_constant(NUM_NPS, any(num_nps));
  	int numcats = categories.size();
 	 sdm.set_constant(NUM_CATS, any(numcats));
   }
  // if (opts.scheduler == "round_robin") {
  //   ((round_robin_scheduler&) (graphlab->get_scheduler())).set_start_vertex(num_nps);
  // }
  
}


// Removes \n from the end of line
void FIXLINE(char * s) {
  int len = strlen(s)-1; 	  
  if(s[len] == '\n') s[len] = 0;
}


std::map<std::string, vertex_id_t> distmulticoemapp::loadVertices(short typemask, FILE * f) {
  /* Allocation size depends on number of categories */
  int vsize = sizeof(vertex_data) + categories.size()*sizeof(float);
  map<std::string, vertex_id_t> map;
  char s[255];
  char delims[] = "\t";	
  char *t = NULL;
     
  std::cout << "Vsize: " << vsize << " typemask:" << typemask <<  std::endl;
     
  // Set same seed
  graphlab::random::seed(999);
  
  while(fgets(s, 255, f) != NULL) {
    // Remove new line
    FIXLINE(s);
    t = strtok(s,delims);

    // Create vertex
    vertex_data * vdata = (vertex_data*) calloc(1, vsize);
    vdata->flags = typemask;
    int len = strlen(t);
    // rather ugly.....
    memcpy(vdata->text, t, (TEXT_LENGTH < len+1 ? TEXT_LENGTH : len+1)); 
        
    // Read count (Not used for anything?)
    //     t = strtok(NULL,delims);
    //  sscanf(t, "%d", &vdata->count);
    double randomNum = graphlab::random::rand01();
		
    if (randomNum <= subsamplingRatio) {
      vertex_id_t vid = (on_the_fly_partition ? distgraph.add_vertex(graphlab::random::rand_int(numofprocs-1), blob(vsize,vdata)) :
                                                    g->add_vertex(blob(vsize, vdata)));
      if (vid%50000 == 0) printf("Vertex: %d  %s\n", vid, vdata->text); 
      map.insert(make_pair(std::string(vdata->text), vid));
    } 
  }
  return map;
}

void distmulticoemapp::loadCategories() {
  /* Iterate files in seeds dir */
  DIR *dp;
  struct dirent *ep;   
  dp = opendir (seedsdir.c_str());
  while(  (ep = readdir (dp))  ) {
    std::string s(ep->d_name);
    if (!((s == ".")|| (s == ".."))) {
      categories.push_back(s);
      std::cout << "Category [" << s << "]" << std::endl;
    }
  }
  closedir(dp);
}

void distmulticoemapp::setseeds(std::map<std::string, vertex_id_t>& nps_map) {
  short catid = 0;
  size_t num_of_seeds=0;
  foreach(std::string cat, categories) {
    std::cout << "Loading seeds for category: " << cat << std::endl;
	 	
    std::string seedfilename = seedsdir + "/" + cat;
    FILE * seedf = fopen(seedfilename.c_str(), "r");
    char cs[255];
    while(fgets(cs, 255, seedf) != NULL) {
      FIXLINE(cs);
      std::string s(cs);
      map<string,vertex_id_t>::iterator iter = nps_map.find(s);
      if (iter == nps_map.end()) {
        //	std::cout << "Seed not found: [" << s << "]" << std::endl;
        continue;
      }
      vertex_data * vdata = (on_the_fly_partition ? distgraph.vertex_data(iter->second).as_ptr<vertex_data>() :
        g->vertex_data(iter->second).as_ptr<vertex_data>());
      //if ((vdata->flags & SEEDMASK)) {
      // 	std::cout << "Warning : Multiple seed: " << s << " was:" << vdata->flags << std::endl;
      //}
      vdata->flags = vdata->flags | SEEDMASK;
      vdata->p[catid] = POSSEEDPROB;
      //	 std::cout << "Pos seed: [" << s << "]" << std::endl;
      num_of_seeds++;
    }	 
    fclose(seedf);
	 	
    // Negative seeds - copypastecode, bad
    std::string negseedfilename = negseedsdir + "/" + cat;
    FILE * negseedf = fopen(negseedfilename.c_str(), "r");
    while(fgets(cs, 255, negseedf) != NULL) {
      FIXLINE(cs);
      std::string s(cs);
      map<string,vertex_id_t>::iterator iter = nps_map.find(s);
      if (iter == nps_map.end()) {
        //std::cout << "Seed not found: [" << s << "]" << std::endl;
        continue;
      }			
      vertex_data * vdata =(on_the_fly_partition ? distgraph.vertex_data(iter->second).as_ptr<vertex_data>() :
        g->vertex_data(iter->second).as_ptr<vertex_data>());
      vdata->p[catid] = NEGSEEDPROB;
      // std::cout << "Neg seed: [" << s << "]" << std::endl;
      //num_of_seeds++;
    }	 
    fclose(negseedf);
    catid++;
  }
  std::cout << "Num of seeds (pos): " << num_of_seeds << std::endl;
}

void distmulticoemapp::loadEdges(FILE * fcont_to_nps, 
                             std::map<std::string, vertex_id_t>& nps_map,
                             std::map<std::string, vertex_id_t>& ctx_map) {
  size_t MAXBUF = 5*1000000;
  char  * s = (char*) malloc(MAXBUF); 
  char delims[] = "\t";	
  char *t = NULL, *t2;
  int vsize = sizeof(vertex_data) + categories.size()*sizeof(float);

  std::cout << "Start to load edges... "<< std::endl;

  /* Get file size*/
  fseek(fcont_to_nps, 0L, SEEK_END);
  long sz = ftell(fcont_to_nps);
  fseek(fcont_to_nps, 0L, SEEK_SET);

  unsigned int i = 0, randj=0;
  unsigned int prunecount = 0;
  std::vector<std::string> * nplist = new std::vector<std::string>();

  while(fgets(s, MAXBUF, fcont_to_nps) != NULL) {
    if (strlen(s) >= MAXBUF-2) {
      // Too long!
      printf("WARNING Skipping too long line!! %s \n", strtok(s, delims));
      // Read while end of line reached 
      while(fgets(s, MAXBUF, fcont_to_nps) != NULL) {
        // Ugly hack	
        if(s[strlen(s)-1] == '\n') break;
      }
      continue;
    }
    FIXLINE(s);
    t = strtok(s, delims);
		
    // Truncate
    if (strlen(t) > TEXT_LENGTH) {
      t[TEXT_LENGTH+1] = 0;
    }
    if (t == NULL) {
      std::cout << "ERROR t was null!!! " << std::endl;
      continue;
    }
    std::string ctxname(t);
    map<string,vertex_id_t>::iterator iter = ctx_map.find(ctxname);
    if (iter == ctx_map.end()) {
      if (subsamplingRatio < 1.0) {
        // Not found (subsampling in effect)
        continue;
      } else {
        // Create vertex
        vertex_data * vdata = (vertex_data*) calloc(1, vsize);
        vdata->flags = CTX_MASK;
        int len = strlen(ctxname.c_str());
        // rather ugly.....
        memcpy(vdata->text, ctxname.c_str(), 
               (TEXT_LENGTH < len+1 ? TEXT_LENGTH : len+1)); 
        vertex_id_t vid = (on_the_fly_partition ?  distgraph.add_vertex(len%numofprocs, blob(vsize, vdata))  : g->add_vertex(blob(vsize, vdata)));
        if (vid%50000 == 0) printf("NEW Vertex: %d  %s\n", vid, vdata->text); 
        ctx_map.insert(make_pair(ctxname, vid));
        iter = ctx_map.find(ctxname);
      }	
    }
    vertex_id_t ctx_vid = iter->second;
		
    //  Each line contains the NP pairs for one pattern, with
    //  format <pattern>[\t<NP> || <NP> -#- <count>[...]]\n
    nplist->clear();
    while ((t = strtok(NULL, delims)) != NULL) {
      nplist->push_back(std::string(t));
    }
 		
    // Parse each np entry on the line
    foreach(std::string npentry, *nplist) {
      t2 = strtok((char*) npentry.c_str(), " ");
      if (t2 == NULL) {
        std::cout << "ERROR t2 was null!!! " << std::endl;
        break;
      }
      std::string npname(t2);
      unsigned short cooccurence_count;
      while(true) {
        t2 = strtok(NULL, " ");
        std::string word(t2);
        if (word == "-#-") {
          t2 = strtok(NULL, " ");
          // This is coocc count
          sscanf(t2, "%hu", &cooccurence_count);
          break;
        } else {
          npname += " " + word; 
        }
      }
 		
 	  ASSERT_GT(cooccurence_count, 0);	
 		
      if (prune_edges && cooccurence_count == 1) {
        // HACK
        if ((ctx_vid+(randj++)) % 3 == 1) {
          prunecount++;
          continue;
        }
      }
	
       map<string,vertex_id_t>::iterator iter = nps_map.find(npname);
		    
      if (iter == nps_map.end()) {
        if (subsamplingRatio < 1.0) {
          continue;
        } else {
          // Create vertex
          vertex_data * vdata = (vertex_data*) calloc(1, vsize);
          vdata->flags = NP_MASK;
          int len = strlen(npname.c_str());
          // rather ugly.....
          memcpy(vdata->text, npname.c_str(), 
                 (TEXT_LENGTH < len+1 ? TEXT_LENGTH : len+1)); 
          vertex_id_t vid = (on_the_fly_partition ? distgraph.add_vertex(len%numofprocs,blob(vsize,vdata)) 
                        : g->add_vertex(blob(vsize,vdata)));
          if (vid%50000 == 0) printf("NEW Vertex: %d  %s\n", vid, vdata->text); 
          nps_map.insert(make_pair(npname, vid));
          iter = nps_map.find(npname);
        }
      }
      vertex_id_t np_vid = iter->second;
	
	   if (!on_the_fly_partition || (distgraph.owner(np_vid) == myprocid || distgraph.owner(ctx_vid) == myprocid)) {				
		  // Create edge data
		  
				
		  if (i++ % 100000 == 0) {
			std::cout << "For ctx [" << ctxname 
					  << "<" << ctx_vid << ">] linked to np: [" << npname << 
			  "<" << np_vid << ">] with weight " << cooccurence_count << std::endl;
			long filepos = ftell(fcont_to_nps);
			printf("Edge #%d (%2.1f perc) pruned:%u\n", i, 
				   filepos*1.0/sz*100, prunecount);
		  }
				
		   if (on_the_fly_partition) {
			  distgraph.add_edge(ctx_vid, np_vid, cooccurence_count);
		   } else {
			  g->add_edge(ctx_vid, np_vid, cooccurence_count);
		   }
       }
    }
  }
  delete(nplist);
  free(s);
}

void distmulticoemapp::loadAndConvertToBin(std::string binfile) {
  /* Open files */
  FILE * fnps = fopen(npsfile.c_str(), "r");
  assert(fnps != NULL);
  FILE * fctx = fopen(contextfile.c_str(), "r");
  assert(fctx != NULL);
  FILE * fmatrix = fopen(matrixfile.c_str(), "r");
  assert(fmatrix != NULL);
	
	
  g = new coem_graph(0); // Todo: read from data num of vertices
    	
  /* Load NPS */
  std::map<std::string, vertex_id_t> nps_map = loadVertices(NP_MASK, fnps);
  /* Load CONTEXTS */
  std::map<std::string, vertex_id_t> ctx_map = loadVertices(CTX_MASK, fctx);

  std::cout << "NPS map " << nps_map.size() << " entries" << std::endl;
  std::cout << "CTX map " << ctx_map.size() << " entries" << std::endl;

  if (on_the_fly_partition)   dc->barrier();
	
  /* Set seeds */
  setseeds(nps_map);
	
  if (on_the_fly_partition)   dc->barrier();
   
	
  /* Load edges */
  loadEdges(fmatrix, nps_map, ctx_map);
  if (on_the_fly_partition)   dc->barrier();

  std::cout << "(After) NPS map " << nps_map.size() << " entries" << std::endl;
  std::cout << "(AFter) CTX map " << ctx_map.size() << " entries" << std::endl;

	
	
  /* Save */
  if (!on_the_fly_partition) {
      std::cout << "Saving graph to disk.... file:"<< binfile << std::endl;
      g->finalize();
      g->save(binfile);
  } else {
     distgraph.finalize();
  }
  fclose(fnps);
  fclose(fctx);
  fclose(fmatrix);
  
}



/**** ANALYZER FUNCTIONS ****/
void calcdiff(std::vector<blob>& x, 
              std::vector<blob>& x_star, 
              double * l1, double * linf) {
  double l1res = 0, l1inf = 0;
  for(unsigned int i=0; i<x.size(); i++) {
    vertex_data *v1 = x_star[i].as_ptr<vertex_data>();
    vertex_data *v2 = x[i].as_ptr<vertex_data>();
		
    for(int catid=0; catid<ncats; catid++) {
      l1res += std::fabs(v1->p[catid] - v2->p[catid]);
      l1inf = std::max((double)std::fabs(v1->p[catid] - v2->p[catid]), l1inf);

    }
    // Blob destructor does not destroy data
    x[i].clear();
  }  
	
  *l1 = l1res;
  *linf = l1inf;
	
  std::cout << "Diff: "<< l1res << " " << l1inf << std::endl;
}

void writeDump(blob_graph * g, 
               bool is_xstar, long int budget=0, 
               double *l1res=NULL, double * linf=NULL) {
  std::vector<blob> dump(0);
  dump.reserve(g->num_vertices());
    
  for(unsigned int i=0; i<g->num_vertices(); i++) {
    dump.push_back(g->vertex_data(i).copy());
  }   
    
  char fname[255];
  if (!is_xstar) sprintf(fname, "mc_dump.%d", dumpcount++);
  else sprintf(fname, "mc_dump.star.%ld", budget);
    
  printf("==========> Dumping to %s\n", fname);
  std::ofstream fout(fname);
  oarchive oarc(fout);
  oarc << dump;
  fout.close();
    
  /* Calculate diff to real x-star */
  if (l1res != NULL && linf != NULL && XSTAR != NULL) {
    calcdiff(dump, *XSTAR, l1res, linf);
  }
    
  foreach(blob b, dump) {
    b.clear();
  }
}



]