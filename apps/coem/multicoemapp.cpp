/*
 *  \author akyrola
 *  Multi category CoEM implementation, based on material by Justin Betteridge
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
 
// lfidc is a modified weight formula for Co-EM (see Justin
// Betteridge's "CoEM results" page)
#define TFIDF(coocc, num_neighbors, vtype_total) (log(1+coocc)*log(vtype_total*1.0/num_neighbors))   

struct vertex_data {
  unsigned short flags;
  int nbcount;
  char text[TEXT_LENGTH];
  float normalizer;
  float p[200];
 
  vertex_data() {
  	memset(p, 0, sizeof(float)*200);
  	normalizer = 0.0f;
  	nbcount = 0;
  	flags = 0;
  }
  
   // TODO: save p!
  void save(graphlab::oarchive& archive) const { 
    archive << flags << nbcount  << normalizer;
    serialize(archive, text, sizeof(char)*TEXT_LENGTH);
   	serialize(archive, p, sizeof(float)*200);
  } 
  
  void load(graphlab::iarchive& archive) { 
    archive >> flags >> nbcount  >> normalizer; 
    deserialize(archive, text, sizeof(char)*TEXT_LENGTH);
    deserialize(archive, p, sizeof(float)*200);
  }
};

struct edge_data {
  unsigned short cooccurence_count;
  
  void save(graphlab::oarchive& archive) const { 
    archive << cooccurence_count; 
  } 
  
  void load(graphlab::iarchive& archive) { 
    archive >> cooccurence_count; 
  }
};

// Graph type defintiions
typedef graphlab::graph<vertex_data, edge_data> coem_graph;
typedef coem_graph graph_type;
typedef graphlab::types<graph_type> gl_types;

// Accuracy 
float TARGET_PRECISION = 1e-5;

// Files
std::string npsfile;
std::string contextfile; 
std::string matrixfile; 
std::string seedsdir;
std::string negseedsdir;

// Category names (entity categories)
std::vector<std::string> categories;	


// Subsampling (to get smaller problems)
double subsamplingRatio = 1.0;
bool prune_edges;
  
// Hack for speeding up
std::vector<vertex_id_t> ** vertex_lookups;

// For convergence analysis
std::vector<vertex_data> * XSTAR;


// Temporary hack 
//multicoemapp * coapp;

// Forward declarations
void loadCategories();
std::map<std::string, vertex_id_t>  loadVertices(gl_types::core& core, short typemask, FILE * f);
void prepare(gl_types::core& core);
void load(gl_types::core& core);
void loadEdges(gl_types::core& core, FILE * fcont_to_nps, 
                             std::map<std::string, vertex_id_t>& nps_map,
                             std::map<std::string, vertex_id_t>& ctx_map);
void setseeds(gl_types::core& core, std::map<std::string, vertex_id_t>& nps_map);
void loadAndConvertToBin(gl_types::core& core, std::string binfile);
void output_results(gl_types::core& core);
void writeDump(coem_graph& g, 
               bool is_xstar, long int budget, 
               double *l1res, double * linf);
std::vector<std::string> dump_headers();
std::vector<double> l1_residual(coem_graph &g);

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
  size_t num_cats = shared_data->get_constant(NUM_CATS).as<size_t>();
   
  /* Get vertex data */
  vertex_data& vdata =   scope.vertex_data();
   
  /* We use directed edges as indirected, so either this vertex has only incoming  
     or outgoing vertices */
  bool use_outgoing = (scope.in_edge_ids().size()==0);
  edge_list edge_ids = 
    (use_outgoing ? scope.out_edge_ids() : scope.in_edge_ids());
   
  if (edge_ids.size() == 0) {
    //printf("Warning : dangling node %d %s\n", scope.vertex(), vdata.text);
    return;
  }	
   
  bool is_np = (vdata.flags & TYPEMASK) == NP_MASK;
  bool was_first_run = false;
   
  unsigned int vtype_total = (is_np ? num_nps : num_ctxs);
  unsigned int vtype_other_total = (is_np ? num_ctxs : num_nps); 
    
  /* First check is my normalizer precomputed. If not, do it. */
  if (vdata.normalizer == 0.0) {
    float norm = param_m * vtype_other_total;
    foreach(edge_id_t eid, edge_ids) {
      const edge_data& edata = scope.const_edge_data(eid);
      const vertex_data& nb_vdata = 
        scope.const_neighbor_vertex_data(use_outgoing ? 
                                   scope.target(eid) :
                                   scope.source(eid));
       norm += TFIDF(edata.cooccurence_count, nb_vdata.nbcount, vtype_total);
    }
    vdata.normalizer = norm;
    was_first_run = true; // Used to cause update to all neighbors
                          // regardless of residual
    if (scope.vertex()%20000 == 0)	
      printf("%d Computed normalizer: %lf \n", scope.vertex(), norm);
  }
	
  /***** COMPUTE NEW VALUE *****/
  float tmp[200];  // assume max 200 cats
  assert(num_cats<200);
  for(int i=0; i<200; i++) tmp[i] = param_m;
	
	
  foreach(edge_id_t eid, edge_ids) {
    const edge_data& edata = scope.const_edge_data(eid);
    vertex_id_t nbvid = use_outgoing ? scope.target(eid) : scope.source(eid);
    const vertex_data& nb_vdata = scope.const_neighbor_vertex_data(nbvid);
    for(unsigned int cat_id=0; cat_id<num_cats; cat_id++) {
      tmp[cat_id] += nb_vdata.p[cat_id] * 
        TFIDF(edata.cooccurence_count, nb_vdata.nbcount, vtype_total);
    }
    vlookup.push_back(nbvid); 
  }
  // Normalize and write data
  float residual = 0.0;
  for(unsigned int cat_id=0; cat_id<num_cats; cat_id++) {
    tmp[cat_id] /= vdata.normalizer;
		 
    /* I am seed for this category? */
    if (!((vdata.p[cat_id] == POSSEEDPROB) || 
          vdata.p[cat_id] == NEGSEEDPROB)) {
      residual += std::fabs(tmp[cat_id] - vdata.p[cat_id]);
      vdata.p[cat_id] = tmp[cat_id];
   }  else {
      // I am seed - do not change
      if (vdata.p[cat_id] == POSSEEDPROB) {
        assert((vdata.flags & SEEDMASK) != 0); 
      }
    }
  }
	
  // Write data to edges and schedule if threshold reached
  double randomNum = graphlab::random::rand01();
  assert(randomNum >= 0.0);
  assert(randomNum <= 1.0);
  if (scope.vertex()%20000 == 0) 
    printf("%d Entering foreach: %lf \n", scope.vertex(), residual);
    
 if (!ROUNDROBIN) {
    int sz = edge_ids.size();
    for(int l = 0; l<sz; l++) {
        vertex_id_t nbvid = vlookup[l];
        gl_types::update_task task(nbvid, coem_update_function);
   		edge_id_t eid = edge_ids[l];
        const edge_data& edata = scope.const_edge_data(eid);
           
        const vertex_data nb_vdata = scope.const_neighbor_vertex_data(nbvid);
        float neighbor_residual = 
           (nb_vdata.normalizer == 0.0f ? 1.0 : 
            residual * 
            TFIDF(edata.cooccurence_count, 
                  vdata.nbcount, vtype_other_total) / 
            nb_vdata.normalizer);	
         if ((neighbor_residual/TARGET_PRECISION >= randomNum) || 
             was_first_run) {
           scheduler.add_task(task, neighbor_residual);
         }
    }
  }
}

/**
  * Used for analysis. Ignore.
  */
void analysis_load_ground_truth() {
	/* Load XSTAR - i.e the "ground truth". Used for convergence anlaysis. */
  if (fopen("MC_XSTAR.dat", "r") != NULL) {
    XSTAR = new vector<vertex_data>();
    std::cout << "Loading x_star..." << std::endl;
    std::ifstream fin("MC_XSTAR.dat");
    iarchive iarc(fin);
    iarc >> *XSTAR;
    fin.close();
  }
}

int ncats = 0;
int dumpcount = 0;
int task_budget;


/**
  *  MAIN FUNCTION 
  */
int main(int argc,  char ** argv) {
  
  /**** GRAPHLAB INITIALIZATION *****/
  std::string root = "/mnt/bigbrofs/usr5/graphlab/testdata/coem/justin/";
  
  // Setup the parser
  graphlab::command_line_options   clopts("Run the CoEM algorithm.");
  
  clopts.attach_option("data_root", &root, root,
                       "Root for data.");
  clopts.attach_option("subsampling_ratio", &subsamplingRatio, subsamplingRatio,
  				"Subsampling ratio.");
  clopts.attach_option("task_budget", &task_budget, 0,
  										"Task budget (needed for round robin).");
  clopts.attach_option("target_precision", &TARGET_PRECISION, 0.0f,
  										"Termination threshold.");
  										
  // TODO: do not hard-code
  npsfile     = root + "/cat_nps.txt";
  contextfile = root + "/cat_contexts.txt";
  matrixfile  = root + "/cat_pairs_cont-idx.txt";
  seedsdir    = root + "/seeds/";
  negseedsdir = root + "/seeds-neg/";
 
  // Create a graphlab core
  gl_types::core core;


  
  if(!clopts.parse(argc, argv)) {
     std::cout << "Error in parsing input." << std::endl;
     return EXIT_FAILURE;
  }
  
  // Set the engine options
  core.set_engine_options(clopts);
  
  
  /**** GRAPH LOADING ****/
  timer t;
  t.start();
  load(core);
  printf("Loading data took %lf secs\n", t.current_time());
  analysis_load_ground_truth();

  
  /** Hack **/
  vertex_lookups = (std::vector<vertex_id_t> **) 
    malloc(sizeof(std::vector<vertex_id_t> *) * clopts.ncpus);
  for(unsigned int i=0; i<clopts.ncpus; i++) {
    vertex_lookups[i] = new std::vector<vertex_id_t>();
    vertex_lookups[i]->reserve(1e6);
  }
  
  prepare(core);
 
  /* Special handling for round_robin */
  if (clopts.scheduler_type == "round_robin") {
    core.add_task_to_all(coem_update_function, 1.0);
    assert(task_budget > 0);  // Need to use budget with round robin
    ROUNDROBIN = true;
  }
  
  if (task_budget > 0) {
 	  core.engine().set_task_budget(task_budget);
  }
  
  
  
  /*** TODO:
  if (clopts.monitor == "analyzer") {
    ncats = categories.size();
    l1_residual(core.graph());
  }
  */
  
  double runtime = core.start();
  
  
  std::cout << "Finished in " << runtime << " seconds." << std::endl;

  //if (opts.visualizer == "analyzer") retrospective_dump();

  std::cout << "Going to output results... " << std::endl;
  
  output_results(core);
  
  if (clopts.scheduler_type == "round_robin" && task_budget >= 20000000) {
    writeDump(core.graph(), true, task_budget, NULL, NULL);
  }  

//  if (clopts.extra == "dump") {
//    writeDump(core.graph(), true, 0, NULL, NULL);
//  }
  
  /*** Write output ***/
  /*
  std::cout << "==== DUMPING FOR JUSTIN ===\n" << std::endl;
  FILE * fnp = fopen("justin_np.txt", "w");
  FILE * fct = fopen("justin_ct.txt", "w");
    
  int n = core.graph().num_vertices();
  loadCategories();
  int ncats = categories.size();
  cout << "Categories: " << ncats << std::endl;
  for(int i=0; i<n; i++) {
    vertex_data& v1 = core.graph().vertex_data(i);
    FILE * f = (((v1.flags & TYPEMASK) == NP_MASK) ? fnp : fct);
    fprintf(f, "%s", v1.text);
    for(int c=0; c<ncats; c++) {
      fprintf(f, "\t%s^^%f", categories[c].c_str(), v1.p[c]);
    }
    fprintf(f, "\n");
    if (i%1000 == 0) fflush(f);
    if (i%500 == 0) printf("%d\n", i);
  }
    
  fclose(fnp);
  fclose(fct);*/
}


// Bit dirty... who cares (maybe Joey :) ). Should use bind()...
int sort_cat_id = 0;
bool cmp(vertex_data * a, vertex_data * b) {
  return a->p[sort_cat_id] > b->p[sort_cat_id]; // Descending order
}

void output_results(gl_types::core& core) {
  int cat_id = 0;
  int n = core.graph().num_vertices();

  foreach(std::string catname, categories) {
    std::string filename = catname + "_results.txt";
    FILE * f = fopen(filename.c_str(), "w");
		
    std::vector<vertex_data *> v = std::vector<vertex_data *>(n);
    int k = 0;
    for(int vid=0; vid<n; vid++) {
      vertex_data& vdata = core.graph().vertex_data(vid);
      if ((vdata.flags & TYPEMASK) == NP_MASK &&
          (vdata.p[cat_id] != POSSEEDPROB)) {
        v[k++] = &core.graph().vertex_data(vid);
      }
    }
    v.resize(k);
    sort_cat_id = cat_id;
    std::cout << "Sorting for category:: " << catname << std::endl;
    sort(v.begin(), v.end(), cmp);
    // Output top
    for(int i=0; i<50; i++) {
      vertex_data vdata = *v[i];
      fprintf(f, "%d: %s  %1.5f\n", i+1, vdata.text, (float) vdata.p[cat_id]); 
    }
		
    fclose(f);
    cat_id++;
  }
}

void load(gl_types::core& core) {
  // Load categories
  loadCategories();
	
  // Check if already converted to binary
  std::string binfile = npsfile + "_graphlabbin.dat";
	
  if (subsamplingRatio != 1.0) {
    char s[255];
    sprintf(s, "%s_%1.2f", binfile.c_str(), subsamplingRatio);
    binfile = std::string(s);
  }
  std::cout << "Binfile: " << binfile << std::endl;
	
  FILE * f = fopen(binfile.c_str(), "r");
  if (f == NULL) {
    std::cout << "Going to convert to => " << binfile << std::endl;
    loadAndConvertToBin(core, binfile);
  } else {
    std::cout << "Loading from binary " << binfile << std::endl;
    core.graph().load(binfile);
  }       
}

void prepare(gl_types::core& core) {
  // Set register shared data
  /* M is used for "smoothing" */
  float m = 0.01;
  core.shared_data().set_constant(PARAM_M, float(m));
  
  /*** Count NPs and CONTEXTs, add SEEDS as initial tasks */
  size_t n = core.graph().num_vertices();
  
  std::cout << "Seed mask " << SEEDMASK << std::endl;
  
  size_t num_contexts = 0;
  size_t num_nps = 0;
  size_t negseeds = 0, posseeds = 0;
  int ntasks = 0;
  for(size_t i = 0; i<n; i++) {
    vertex_data& vdata = core.graph().vertex_data(i);
  	  
    /* Is seed? */
    if ((vdata.flags & SEEDMASK) != 0) {
      for(unsigned int cat_id=0; cat_id<categories.size(); cat_id++) {
        if (vdata.p[cat_id] == POSSEEDPROB) {
          // 	graphlab->get_scheduler().add_task(update_task(i, coem_update_function), 1.0); 
          //	ntasks++;
          posseeds++;
        }
      }
  		 
    }
  	  
    /* Count classes */
    if ((vdata.flags & TYPEMASK) == NP_MASK) {
      if (i%10000 == 0) std::cout << "NP: " << std::string(vdata.text) << std::endl;
      num_nps++;
  	  	
      /* Set initial values */
      for(unsigned int cat_id=0; cat_id<categories.size(); cat_id++) {
        if (vdata.p[cat_id] == 0.0f) {
          vdata.p[cat_id] = INITIAL_NP_VALUE;
        } 
        else if (vdata.p[cat_id] == NEGSEEDPROB) {
          //	graphlab->get_scheduler().add_task(update_task(i, coem_update_function), 1.0); 
          //	ntasks++;
          negseeds++;
        }
      }
  	  	
    } else  if ((vdata.flags & TYPEMASK) == CTX_MASK) {
      if (i%10000 == 0) std::cout << "CTX: " << std::string(vdata.text) << std::endl;
      num_contexts++;
      // Start with contexts
      
      if (!ROUNDROBIN) {
      	core.add_task(gl_types::update_task(i, coem_update_function), 1.0); 
      }
      ntasks++;
    } else assert(false);
  	  
    /* Count neighbors */
    vdata.nbcount = std::max(core.graph().in_edge_ids(i).size(), core.graph().out_edge_ids(i).size());
    if (i%10000 == 0) std::cout << "Num neighbors: " << vdata.nbcount << std::endl;

  }
  std::cout << "Positive seeds: " << posseeds << std::endl;
  std::cout << "Negative seeds: " << negseeds << std::endl;
  std::cout << "Num of initial tasks: " << ntasks << std::endl;
  std::cout << "Contexts: " << num_contexts << " nps: " << num_nps << std::endl;
  
  core.shared_data().set_constant(NUM_CTX, size_t(num_contexts));
  core.shared_data().set_constant(NUM_NPS, size_t(num_nps));
  int numcats = categories.size();
  core.shared_data().set_constant(NUM_CATS, size_t(numcats));
}


// Removes \n from the end of line
void FIXLINE(char * s) {
  int len = strlen(s)-1; 	  
  if(s[len] == '\n') s[len] = 0;
}


std::map<std::string, vertex_id_t>  loadVertices(gl_types::core& core, short typemask, FILE * f) {
  /* Allocation size depends on number of categories */
  int vsize = sizeof(vertex_data) + categories.size()*sizeof(float);
  map<std::string, vertex_id_t> map;
  char s[255];
  char delims[] = "\t";	
  char *t = NULL;
     
  std::cout << "Vsize: " << vsize << " typemask:" << typemask <<  std::endl;
     
  while(fgets(s, 255, f) != NULL) {
    // Remove new line
    FIXLINE(s);
    t = strtok(s,delims);

    // Create vertex
    vertex_data vdata = vertex_data();
    vdata.flags = typemask;
    int len = strlen(t);
    // rather ugly.....
    memcpy(vdata.text, t, (TEXT_LENGTH < len+1 ? TEXT_LENGTH : len+1)); 
        
    // Read count (Not used for anything?)
    //     t = strtok(NULL,delims);
    //  sscanf(t, "%d", &vdata.count);
    double randomNum = graphlab::random::rand01();
		
    if (randomNum <= subsamplingRatio) {
      vertex_id_t vid = core.graph().add_vertex(vdata);
      if (vid%50000 == 0) printf("Vertex: %d  %s\n", vid, vdata.text); 
      map.insert(make_pair(std::string(vdata.text), vid));
    } 
  }
  return map;
}

void loadCategories() {
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

void setseeds(gl_types::core& core, std::map<std::string, vertex_id_t>& nps_map) {
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
        continue;
      }
      vertex_data& vdata = core.graph().vertex_data(iter->second);
      vdata.flags = vdata.flags | SEEDMASK;
      vdata.p[catid] = POSSEEDPROB;
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
      vertex_data& vdata = core.graph().vertex_data(iter->second);
      vdata.p[catid] = NEGSEEDPROB;
    }	 
    fclose(negseedf);
    catid++;
  }
  std::cout << "Num of seeds (pos): " << num_of_seeds << std::endl;
}

void loadEdges(gl_types::core& core, FILE * fcont_to_nps, 
                             std::map<std::string, vertex_id_t>& nps_map,
                             std::map<std::string, vertex_id_t>& ctx_map) {
  size_t MAXBUF = 5*1000000;
  char  * s = (char*) malloc(MAXBUF); 	// 10 meg buffer :)
  char delims[] = "\t";	
  char *t = NULL, *t2;

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
        vertex_data vdata = vertex_data();
        vdata.flags = CTX_MASK;
        int len = strlen(ctxname.c_str());
        // rather ugly.....
        memcpy(vdata.text, ctxname.c_str(), 
               (TEXT_LENGTH < len+1 ? TEXT_LENGTH : len+1)); 
        vertex_id_t vid = core.graph().add_vertex(vdata);
        if (vid%50000 == 0) printf("NEW Vertex: %d  %s\n", vid, vdata.text); 
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
      unsigned int cooccurence_count;
      while(true) {
        t2 = strtok(NULL, " ");
        std::string word(t2);
        if (word == "-#-") {
          t2 = strtok(NULL, " ");
          // This is coocc count
          sscanf(t2, "%u", &cooccurence_count);
          break;
        } else {
          npname += " " + word; 
        }
      }
 		
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
          vertex_data vdata = vertex_data();
          vdata.flags = NP_MASK;
          int len = strlen(npname.c_str());
          // rather ugly.....
          memcpy(vdata.text, npname.c_str(), 
                 (TEXT_LENGTH < len+1 ? TEXT_LENGTH : len+1)); 
          vertex_id_t vid = core.graph().add_vertex(vdata);
          if (vid%50000 == 0) printf("NEW Vertex: %d  %s\n", vid, vdata.text); 
          nps_map.insert(make_pair(npname, vid));
          iter = nps_map.find(npname);
        }
      }
      vertex_id_t np_vid = iter->second;

 			
      // Create edge data
      edge_data edata;
      edata.cooccurence_count = cooccurence_count;
			
      if (i++ % 100000 == 0) {
        std::cout << "For ctx [" << ctxname 
                  << "<" << ctx_vid << ">] linked to np: [" << npname << 
          "<" << np_vid << ">] with weight " << cooccurence_count << std::endl;
        long filepos = ftell(fcont_to_nps);
        printf("Edge #%d (%2.1f perc) pruned:%u\n", i, 
               filepos*1.0/sz*100, prunecount);
      }
			
			
			
      core.graph().add_edge(ctx_vid, np_vid, edata);
    }
  }
  delete(nplist);
  free(s);
}

void loadAndConvertToBin(gl_types::core& core, std::string binfile) {
  /* Open files */
  FILE * fnps = fopen(npsfile.c_str(), "r");
  assert(fnps != NULL);
  FILE * fctx = fopen(contextfile.c_str(), "r");
  assert(fctx != NULL);
  FILE * fmatrix = fopen(matrixfile.c_str(), "r");
  assert(fmatrix != NULL);
    	
  /* Load NPS */
  std::map<std::string, vertex_id_t> nps_map = loadVertices(core, NP_MASK, fnps);
  /* Load CONTEXTS */
  std::map<std::string, vertex_id_t> ctx_map = loadVertices(core, CTX_MASK, fctx);

  std::cout << "NPS map " << nps_map.size() << " entries" << std::endl;
  std::cout << "CTX map " << ctx_map.size() << " entries" << std::endl;

	
  /* Set seeds */
  setseeds(core, nps_map);
	
  /* Load edges */
  loadEdges(core, fmatrix, nps_map, ctx_map);
	
  std::cout << "(After) NPS map " << nps_map.size() << " entries" << std::endl;
  std::cout << "(AFter) CTX map " << ctx_map.size() << " entries" << std::endl;

	
  std::cout << "Saving graph to disk.... file:"<< binfile << std::endl;
	
  /* Save */
  core.graph().finalize();
  core.graph().save(binfile);
	
  fclose(fnps);
  fclose(fctx);
  fclose(fmatrix);
	
  exit(0);
}



/**** ANALYZER FUNCTIONS ****/
 void calcdiff(std::vector<vertex_data>& x, 
              std::vector<vertex_data>& x_star, 
              double * l1, double * linf) {
  double l1res = 0, l1inf = 0;
  for(unsigned int i=0; i<x.size(); i++) {
    vertex_data& v1 = x_star[i];
    vertex_data& v2 = x[i];
		
    for(int catid=0; catid<ncats; catid++) {
      l1res += std::fabs(v1.p[catid] - v2.p[catid]);
      l1inf = std::max((double)std::fabs(v1.p[catid] - v2.p[catid]), l1inf);

    }
  }  
	
  *l1 = l1res;
  *linf = l1inf;
	
  std::cout << "Diff: "<< l1res << " " << l1inf << std::endl;
}

void writeDump(coem_graph& g, 
               bool is_xstar, long int budget=0, 
               double *l1res=NULL, double * linf=NULL) {
  std::vector<vertex_data> dump(0);
  dump.reserve(g.num_vertices());
    
  for(unsigned int i=0; i<g.num_vertices(); i++) {
    dump.push_back(g.vertex_data(i));
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
}



std::vector<double> l1_residual(coem_graph &g) {
  double l1res=0, linf=0;

  // Dumps a snapshot to disk
  writeDump(g, false, 0, &l1res, &linf);
   
  // Just fake numbers - actual residuals computed after convergence
  std::vector<double> v;
  v.push_back(log(l1res/g.num_vertices()));
  v.push_back(log(linf));
  return v;
}

/* Used with analyzer_listener */
//global_dumper dump_function() {
//  return multicoemdumper;
//}

std::vector<std::string> dump_headers() {
  std::vector<std::string> h;
  h.push_back("l1_residual");
  h.push_back("inf_residual");
  return h;
}

int dump_frequency() {
  return 50000; 
} 

void retrospective_dump() {
  std::string filename = "multicoem_dump2.dat";
  FILE * dumpfile = fopen(filename.c_str(), "w");
  /* Write headers */
  fprintf(dumpfile, "updates\t");
  for(unsigned int i=0; i<dump_headers().size(); i++) {
    fprintf(dumpfile, "\t");
    fprintf(dumpfile, "%s", dump_headers()[i].c_str());
  }
  fprintf(dumpfile,"\n");
  std::vector<vertex_data> x_star;
    
  char fname[255];

  std::cout << "Reading x*" << std::endl;
  int lastdumpid = dumpcount-1;
  {
    sprintf(fname, "mc_dump.%d", lastdumpid);
    std::ifstream fin(fname);
    iarchive iarc(fin);
    iarc >> x_star;
    fin.close();
  }
    
    
  for(int dumpid=0; dumpid<lastdumpid; dumpid++) {
    std::vector<vertex_data> x;
    sprintf(fname, "mc_dump.%d", dumpid);
    std::cout << "Reading " << fname << std::endl;
    std::ifstream fin(fname);
    iarchive iarc(fin);
    iarc >> x;
    fin.close();
        
    std::cout << "x size " << x.size() << std::endl;
        
    /* Calculate l1 and linf */
    double l1res = 0, l1inf = 0;
   		
    calcdiff(x, x_star, &l1res, &l1inf);
   		
    printf("l1res=%lf\n", l1res);
    l1res = l1res/x_star.size();
    fprintf(dumpfile, "%d\t%lf\t%lf\n", dumpid*dump_frequency(), log(l1res), log(l1inf));
        
  }
  fclose(dumpfile);
}


/* Writes current graph L1 norm and max value */
std::vector<double> multicoemdumper(coem_graph &g) {
  return l1_residual(g);
}


