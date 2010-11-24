/*
 *  \author akyrola
 *  CoEM implementation, based on ReadTheWeb course source material by Tom Mitchell.
 *
 */

#include <cmath>
#include <cstdio>


#include <graphlab.hpp>



typedef graphlab::types<graphlab::blob_graph> gl_types;


#include "coemapp.hpp"


#include <graphlab/macros_def.hpp>

using namespace std;
using namespace graphlab;


// Allow easy testing with float/double.
typedef double probability_t;

bool useSumResiduals = false;


size_t memory_writes = 0;
size_t memory_reads = 0;


struct coem_vertex_data {
  uint8_t flags;
  probability_t prob;
  int iterations;
  char text[64];
};

struct coem_edge_data {
  probability_t prob;
  probability_t target_norm;
  probability_t lastusedval;
  unsigned int cooccurence_count;
};


double TARGET_PRECISION = 1e-5;
int NUM_ITERATIONS = 10000;

uint8_t NOUN_PHRASE_FLAG = 1;
uint8_t CONTEXT_FLAG = 1 << 1;
uint8_t POSITIVESEED_FLAG = 1 << 2;
uint8_t NEGATIVESEED_FLAG = 1 << 3;

// Temporary hack 
coemapp * coapp;


/// ====== UPDATE FUNCTION ======= ///
void coem_update_function(gl_types::iscope& scope, 
                          gl_types::icallback& scheduler,
                          gl_types::ishared_data* shared_data) {
  //scheduler.disable_buffering();
  probability_t p = 0.0;
  coem_vertex_data * vdata =
    scope.vertex_data().as_ptr<coem_vertex_data>();
    
  // Iteration check
  if (vdata->iterations >= NUM_ITERATIONS-1) return;
    
  size_t reads = 8;
  size_t writes = 8;
    
  // Check seeds
  bool isseed = false;
  if ((vdata->flags & POSITIVESEED_FLAG) != 0) {
    p = 1.0;
    isseed = true;
  } else if ((vdata->flags & NEGATIVESEED_FLAG) != 0) {
    p = 0.0;
    isseed = true;
  } else {
    // Calculate EM step
    unsigned int norm = 0;
    foreach(edge_id_t eid, scope.in_edge_ids()) {        
      coem_edge_data * edata =
        scope.edge_data(eid).as_ptr<coem_edge_data>();
      p += edata->cooccurence_count * edata->prob;
      norm += edata->cooccurence_count;
        
      // Update residual and last used value
      edata->lastusedval = edata->prob;
    }
    reads += scope.in_edge_ids().size() *
      (sizeof(edge_id_t) + sizeof(probability_t)); // prob and edge_id
    writes +=  scope.in_edge_ids().size() * sizeof(probability_t); // lastusedval
   
    p = p/norm;
    if (norm == 0) {
      logger(LOG_INFO, "What??, %d", scope.in_edge_ids().size());
    }
  }
    
  // Schedule outbound edges
  vdata->prob = p;
  writes += sizeof(probability_t);
    
  /* EXPERIMENTAL: Create random number to use for stochastic scheduling.
     Randomize only once for vertex for performance reasons. Have to think if it
     it is ok. */
  double randomNum = 0.0;
  if (!useSumResiduals) randomNum = graphlab::random::rand01();
  assert(randomNum >= 0.0);
  assert(randomNum <= 1.0);
    
  foreach(edge_id_t eid, scope.out_edge_ids()) {        
    coem_edge_data * edata =
      scope.edge_data(eid).as_ptr<coem_edge_data>();
    edata->prob = p;
    double residual = std::fabs(edata->prob-edata->lastusedval) *
      edata->cooccurence_count/edata->target_norm;
    if (useSumResiduals || residual/TARGET_PRECISION >= randomNum ||
        (isseed && vdata->iterations == 0)) {
      scheduler.add_task(gl_types::update_task(scope.target(eid), 
                                                 coem_update_function),
                         residual );
      writes += sizeof(gl_types::update_task) + 
        sizeof(double) + 8; // 8 is for the pointer

    }
  }
  writes += scope.out_edge_ids().size() *
    (sizeof(edge_id_t) + sizeof(probability_t));
  reads += scope.out_edge_ids().size() *
    sizeof(double); // last_used_val

  vdata->iterations++;
    
  // Debug
  if (scope.vertex() % 20000 == 0) {
    printf( "vertex %u iter %d, prob = %lf\n",
            (unsigned int) scope.vertex(),
            vdata->iterations,  (double) vdata->prob);
  }
  __sync_add_and_fetch(&memory_reads, reads);
  __sync_add_and_fetch(&memory_writes, writes);
}



// Comparator
bool cmp( coem_vertex_data * a, coem_vertex_data * b ) {
  return a->prob > b->prob; // Descending order
}


// Outputs top vertices of given type
void coemapp::outputTopVertices(blob_graph &g, uint8_t type) {
  int n = g.num_vertices();
  int k = 0;
  std::vector<coem_vertex_data *> v = std::vector<coem_vertex_data *>(n);
  for(int vid=0; vid<n; vid++) {
    coem_vertex_data * vdata = 
      g.vertex_data(vid).as_ptr<coem_vertex_data>();
    if ((vdata->flags & POSITIVESEED_FLAG) == 0 &&
        (vdata->flags & type) != 0) {
      v[k++] = g.vertex_data(vid).as_ptr<coem_vertex_data>();
    }
  }
  v.resize(k);
  sort(v.begin(), v.end(), cmp);
  
  // Output top
  for(int i=0; i<40; i++) {
    coem_vertex_data * vdata = v[i];
    printf("%d: %s  %1.5f\n", i+1, vdata->text, (float) vdata->prob); 
  }
}



// Helper to load vertices
char s[255];
int coemapp::loadVertices(blob_graph * g, const  char * filename, uint8_t vtype) {
  printf("Loading.. %s\n", filename);
  set_status("Loading: %s", filename); 
  int lastid = -1;
  FILE * f = fopen(filename, "r");
  if (f == NULL) {
    set_status("File %s not found - exiting", filename);
    exit(1);
  }
  
  while(fgets(s, 255, f) != NULL) {
    // Remove new line
    int len = strlen(s)-1;
    if(s[len] == '\n') s[len] = 0;
    
    // Set vertex fields
    coem_vertex_data vdata;
    vdata.flags = vtype;
    memcpy(vdata.text, s, (64 < len+1 ? 64 : len+1)); // rather ugly.....
    vdata.iterations = 0;
    vdata.prob = 0.0;
    
    lastid = g->add_vertex(blob(sizeof(coem_vertex_data), &vdata));
    if (lastid % 10000 == 0) {
      set_status("Load vertex id: %d", lastid);
    }
  }
  fclose(f);
  return lastid;
}

// Loads a file of positive or negative seeds and marks
// appropiate vertices.
void  coemapp::loadSeeds(const char * filename, bool positive_seed, 
                         blob_graph &g) {
  FILE * f = fopen(filename, "r");
  if (f == NULL) {
    set_status( "File %s not found - exiting", filename);
    exit(1);
  }
  

  while(fgets(s, 255, f) != NULL) {
    int len = strlen(s)-1;
    
    // Remove regex stuff from the strings, we just match
    if(s[len] == '\n') s[len--] = 0; // Remove trailing \n
    if(s[len] == '$') s[len] = 0;
    
    int n = g.num_vertices();
    // look for correct vertex and set appropiate seed flag
    char * t = &s[1]; // Hack to remove leading regex '^'
    logger(LOG_INFO, "Looking for seed: %s", t);

    for(int vid=0; vid<n; vid++) {
      coem_vertex_data * vd = 
        g.vertex_data(vid).as_ptr<coem_vertex_data>();
      if (strcmp(t, vd->text) == 0) {
        vd->flags |= (positive_seed ? POSITIVESEED_FLAG : NEGATIVESEED_FLAG);
        logger(LOG_INFO, "Found seed: %s, flag -> %d", t, vd->flags);
        break;
      } 
    } 
  }
  fclose(f);
}

void  coemapp::loadCoEMDataset(const char * npsfile, 
                               const  char * contextfile, 
                               const  char * matrixfile, 
                               const  char * seedsfile, 
                               const  char * negseedsfile) {
  
  printf("Creating graph\n");
  
  /* Register needed fields */
  g = new blob_graph();

  /* Load noun phrases and contexts - and make them as vertices */
  int last_nps_id = loadVertices(g, npsfile, NOUN_PHRASE_FLAG);
  loadVertices(g, contextfile, CONTEXT_FLAG);
  
  /* Load seeds */
  loadSeeds(seedsfile, true, *g);
  loadSeeds(negseedsfile, false, *g);

  /* Create edges (co-occurences) */
  set_status("Loading matrix file");
  FILE * fmatrix = fopen(matrixfile, "r");
  if (fmatrix == NULL) {
    set_status("File %s not found - exiting", matrixfile);
    exit(1);
  }
  
  char delims[] = "\t";
  int j=0;
  char *t = NULL;

  /* Get file size - for debug - can be removed later */
  fseek(fmatrix, 0L, SEEK_END);
  long sz = ftell(fmatrix);
  fseek(fmatrix, 0L, SEEK_SET);
  
  while(fgets(s, 255, fmatrix) != NULL) {
    int v1, v2;
    int w;
    t = strtok(s,delims);
    sscanf(t, "%d", &v1);
    t = strtok(NULL, delims);
    sscanf(t, "%d", &v2);
    t = strtok(NULL, delims);
    sscanf(t, "%d", &w);
      
    // Debug
    if (++j%100000==0) {
      long filepos = ftell(fmatrix);
      printf("Edge #%d (%2.1f perc)\n", j, filepos*1.0/sz*100); 
    }
   
    /* Create edge between context and nps - two ways  */
    coem_edge_data edata;
    edata.cooccurence_count = w;
    edata.prob = 0;
    edata.lastusedval = 10000;
    g->add_edge(v1-1, last_nps_id+v2, blob(sizeof(coem_edge_data), &edata)); 

    edata.cooccurence_count = w;
    edata.prob = 0;
    edata.lastusedval = 10000;
    g->add_edge(last_nps_id+v2, v1-1, blob(sizeof(coem_edge_data), &edata)); 

  }

  fclose(fmatrix);
  
  /* Compute normalizers */
  printf("Computing normalizers\n");
  int n = g->num_vertices();
  for(int vid=0; vid<n; vid++) {
    probability_t norm = 0.0;
    /* Sum up all cooccurence_counts for incoming edges */
    foreach(edge_id_t eid, g->in_edge_ids(vid)) {
      coem_edge_data * edata =
        g->edge_data(eid).as_ptr<coem_edge_data>();
      norm += edata->cooccurence_count;
    }
	  
    /* Write norm to all outgoing edges */
    foreach(edge_id_t eid, g->out_edge_ids(vid)) {
      coem_edge_data * edata =
        g->edge_data(eid).as_ptr<coem_edge_data>();
      edata->target_norm = norm;     	
    }
  }	
  printf("Going to finalize graph...\n");
  g->finalize();
}

coemapp::coemapp(const char * npsfile, const char * contextfile, 
                 const char * matrixfile, const char * seedsfile, 
                 const char * negseedsfile) {
  this->npsfile = npsfile;
  this->contextfile = contextfile;
  this->matrixfile = matrixfile;
  this->seedsfile = seedsfile;
  this->negseedsfile = negseedsfile;
  coapp = this;
}

coemapp::~coemapp() {
  delete(g);
}


/**** GLOBAL CONVERGENCE CALCULATOR ****/
// Calculates |Ax-x|_1
std::vector<double> l1_residual(blob_graph &g) {
  unsigned int n = g.num_vertices();
  double l1res = 0.0;
  double l1inf = 0.0;
  double l1norm = 0.0;
  for(vertex_id_t vid=0; vid<n; vid++) {
    coem_vertex_data * vdata =
      g.vertex_data(vid).as_ptr<coem_vertex_data>();

    // Check seeds
    double p = 0.0;
    bool isseed = false;
    if ((vdata->flags & POSITIVESEED_FLAG) != 0) {
      p = 1.0;
      isseed = true;
    } else if ((vdata->flags & NEGATIVESEED_FLAG) != 0) {
      p = 0.0;
      isseed = true;
    } else {
      // Calculate EM step
      unsigned int norm = 0;
      foreach(edge_id_t eid, g.in_edge_ids(vid)) {        
        coem_edge_data * edata =
          g.edge_data(eid).as_ptr<coem_edge_data>();  
        p += edata->cooccurence_count * edata->prob;
        norm += edata->cooccurence_count;
      }
      p = p/norm;
      l1norm += vdata->prob;
      l1res += std::fabs(p-vdata->prob);
      l1inf = std::max(l1inf, std::fabs(p-vdata->prob));
    }
  }
    
  printf("L1 norm = %lf\n", l1norm);
  printf("|A-x|_inf = %lf\n", log(l1inf/l1norm));
    
  std::vector<double> v;
  v.push_back(log(l1res/l1norm));
  v.push_back(log(l1inf/l1norm));
    
  return v;
}



void coemapp::start() {  
  /**** GRAPH LOADING ****/
  loadCoEMDataset(npsfile, contextfile, matrixfile, seedsfile, negseedsfile);

  /**** GRAPHLAB INITIALIZATION *****/
  iengine<blob_graph> * graphlab = create_engine(g, NULL);
  
  int n = g->num_vertices();
  timer t2;
  t2.start();
  
  TARGET_PRECISION = opts.threshold;
  typedef gl_types::ischeduler ischeduler_type;
  graphlab->set_option("update_function",coem_update_function);
  if (typeid(*graphlab) == typeid(synchronous_engine<blob_graph>)) {
    ((synchronous_engine<blob_graph>*)(graphlab))->
      set_update_function(coem_update_function);
  } else {
    if (opts.scheduler == "round_robin") {
      graphlab->add_task_to_all(coem_update_function, 1.0);
    } else {
      /* Initial tasks - only positive seeds */
      for(int vid=0; vid<n; vid++) {
        coem_vertex_data * vdata =
          g->vertex_data(vid).as_ptr<coem_vertex_data>();
        if ((vdata->flags & POSITIVESEED_FLAG) != 0) {
          graphlab->add_task(gl_types::update_task(vid, coem_update_function), 1.0);
          set_status("Add initial task: %d\n", vid);
        }
        //   printf("Flag: %d %d", vid, vdata->flags); 
      }
    }
    printf("Add tasks: %lf\n", t2.current_time());
  }
  timer t;
  t.start();
  
  
  // if (opts.scheduler == "linear_threshold") {
  //   std::cout << "Use linear linear_threshold scheduler" << std::endl;
  //   useSumResiduals = true;

  //   // Set TARGET_PRECISION
  //   ((linear_threshold_scheduler&) graphlab->get_scheduler()).set_threshold(1e-5);
  // }	
  
  
  /* Start */
  graphlab->start();
  
  printf("Finished: %lf secs\n", t.current_time());
  
  /* Output l1 residual for benchmark */
  write_benchmark_value("residual", l1_residual(*g)[0]); 
  write_benchmark_value("memory_writes_mb", memory_writes * 1.0 / 1024.0/1024.0); 
  write_benchmark_value("memory_reads_mb", memory_reads * 1.0 / 1024.0/1024.0); 

  // Outout
  outputTopVertices(*g, NOUN_PHRASE_FLAG);
  delete(graphlab);
}


/**** ANALYZER FUNCTIONS ****/

/* Writes current graph L1 norm and max value */
std::vector<double> coemdumper(blob_graph &g) {
  return l1_residual(g);
}


/* Used with analyzer_listener */
global_dumper coemapp::dump_function() {
  return coemdumper;
}

std::vector<std::string> coemapp::dump_headers() {
  std::vector<std::string> h;
  h.push_back("l1_residual");
  h.push_back("linf_residual");
  return h;
}

int coemapp::dump_frequency() {
  return 100000; // Every 100000 updates
}

