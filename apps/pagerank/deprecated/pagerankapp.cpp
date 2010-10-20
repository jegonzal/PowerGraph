/*
 * \author akyrola
 */
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <graphlab.hpp>

#include "pagerankapp.hpp"

#include <graphlab/macros_def.hpp>

using namespace graphlab;


size_t memory_writes = 0;
size_t memory_reads = 0;



#define TARGET_PRECISION 1e-30

#define RANDOM_SUMRESIDUALS 1
#define SUMRESIDUALS 2
#define DEGREERESIDUALSWITHMULTIPLE 3
#define DEGREERESIDUALSNOMULTIPLE 4

bool DISABLE_ADD_TASKS = false;
int TASK_SCHEDULING_TYPE = RANDOM_SUMRESIDUALS;


float jumpweight = 0.1;
pagerank_graph * __g;


const size_t NUM_VERTICES = 1;
const size_t RANDOM_JUMP_WEIGHT = 2;

gl_types::imonitor * _listener(NULL);

//void threshold_limit_scheduler(gl_types::set_scheduler &sched);





/**** GLOBAL CONVERGENCE CALCULATOR ****/
// Calculates |Ax-x|_1
double l1_residual(gl_types::graph &g) {
  unsigned int n = g.num_vertices();
  double l1res = 0.0;
  double l1inf = 0.0;
  double l1norm = 0.0;
  for(vertex_id_t vid=0; vid<n; vid++) {
    vertex_data vdata = g.vertex_data(vid);

    /* Calculate Ax_i */
    double x = vdata.value * vdata.selfweight;

    foreach(edge_id_t eid, g.in_edge_ids(vid)) {
      edge_data edata = g.edge_data(eid);
      x += edata.weight * edata.srcvalue;
    }
    // compute the jumpweight
    x=jumpweight/n + (1- jumpweight)*x;
    l1norm+=x;
    l1res += std::fabs(x-vdata.value);
    l1inf = std::max(std::fabs(x-vdata.value), l1inf);
    if (std::fabs(x-vdata.value) > 1e-4) {
      //printf("%lf, %lf => %lf\n", std::fabs(x-vdata->value),vdata->value,x);
    }
  }
        
  printf("L1 norm = %lf\n", l1norm);
  printf("|A-x|_inf = %lf\n", log(l1inf/l1norm));
  return log(l1res/l1norm);
}




void pagerankapp::loadGraph(std::string filename) {
  // if this is a text file
  if (filename.substr(filename.length()-3,3)=="txt") {
    std::ifstream fin(filename.c_str());
    size_t edgecount = 0;
    while(fin.good()) {
      // read a line. if it begins with '#' it is a comment
      std::string line;
      std::getline(fin, line);
      if (line[0] == '#') {
        std::cout << line << std::endl;
        continue;
      }

      // read the vertices
      std::stringstream s(line);
      size_t srcv, destv;
      s >> srcv >> destv;
      // make sure we have all the vertices we need
      while (g.num_vertices() <= srcv || g.num_vertices() <= destv) {
        vertex_data v;
        v.value = 1.0;
        v.selfweight = 0.0;
        g.add_vertex(v);
      }
      // add the edge. BUT
      // The graph has multiedges, AND selfedges.. this is unbelievably annoying
      if (srcv == destv) {
        // selfedge
        vertex_data& v = g.vertex_data(srcv);
        // we use selfweight temporarily to store the number of self edges
        v.selfweight += 1.0;
      }
      else {
        // check if the edge already exists
        std::pair<bool, edge_id_t> edgefind = g.find(srcv, destv);
        if (edgefind.first) {
          edge_data& edata = g.edge_data(edgefind.first);
          // if it does, increment the weight
          edata.weight += 1.0;
        }
        else {
          edge_data e;
          e.srcvalue = 1.0;
          // we use weight temporarily as a counter
          e.weight = 1.0; 
          e.residual = 0.0;
          g.add_edge(srcv, destv, e);
        }
      }
      ++edgecount;
      if (edgecount % 1000000 == 0) {
        std::cout << edgecount << " Edges inserted" << std::endl;
      }
    }
    // now we have to go through the edges again and figure out the
    // weights
    for (vertex_id_t i = 0;i < g.num_vertices(); ++i) {
      vertex_data& v = g.vertex_data(i);
      // count the effective number of out edges
      double weight = 0;
      foreach(edge_id_t e, g.out_edge_ids(i)){
        edge_data edata = g.edge_data(e);
        weight += edata.weight;
      }
      // remember to count selfweights
      weight += v.selfweight;
      if (weight == 0) continue;
      // now the weights should all be scaled to 1
      foreach(edge_id_t e, g.out_edge_ids(i)) {
        edge_data& edata = g.edge_data(e);
        edata.weight /= weight;
      }
      v.selfweight /= weight;
    }
    g.finalize();
  } else {
    g.load(filename);
  }
}


pagerankapp::pagerankapp(std::string ifile, std::string bfile, bool optimize) {
  inputfile = ifile;
  binoutfile = bfile;
  graphoptimize=optimize;
}

pagerankapp::~pagerankapp() {
  delete graphlab;
}


void pagerankapp::start() {
  loadGraph(inputfile);
  if (binoutfile.length() > 0) {
    g.save(binoutfile);
  }
  std::cout << "Graph has " << g.num_vertices() << " vertices and "
            << g.num_edges() << " edges." <<std::endl;
  const size_t NUM_VERTICES = 1;
  const size_t RANDOM_JUMP_WEIGHT = 2;
  
  size_t numv = g.num_vertices();

  // Set register shared data
  gl_types::thread_shared_data sdm;
  sdm.set_constant(NUM_VERTICES, numv);
  sdm.set_constant(RANDOM_JUMP_WEIGHT, jumpweight);
  graphlab = create_engine(&g, NULL);
  
  // Hack
  __g = &g;
  
  if (opts.extra == "indegree_multiple") {
    printf(" ==== Using indegree multiple in residual computation ====\n");
    TASK_SCHEDULING_TYPE = DEGREERESIDUALSWITHMULTIPLE;
  } else if (opts.extra == "simple_residual") {
    printf(" ==== Using simple degree residuals ====\n");
    TASK_SCHEDULING_TYPE = DEGREERESIDUALSNOMULTIPLE;
  }
  

  // register shared data
  graphlab->set_shared_data_manager(&sdm);
  
  _listener = this->listener;
  

  for(vertex_id_t vid=0; vid<g.num_vertices(); vid+=100000) {
    vertex_data vdata = g.vertex_data(vid);
    printf("%d, log:%lf\n ", vid, log(vdata.value));
  }
  
  timer t2;
  t2.start();



  if (typeid(*graphlab) == typeid(synchronous_engine<gl_types::graph>)) {
    ((synchronous_engine<gl_types::graph>*)(graphlab))->
      set_update_function(pagerank_function);
  } else {
    graphlab->get_scheduler().set_option(scheduler_options::UPDATE_FUNCTION,
                                       (void*)pagerank_function);

      //    graphlab->get_scheduler().set_update_function(pagerank_function);
    // Removed by Aapo: this will not use the sorting. TODO: use sorting
    /* for (size_t i = 0;i < g.num_vertices(); ++i) {
       if (g.in_edge_ids(i).size() > 0 || g.out_edge_ids(i).size() > 0 ) {
       graphlab->get_scheduler().add_task(update_task(i, pagerank_function), 1000.0);
       ++numvertices;
       }
       }
       std::cout << "Effective num vertices = " << numvertices << std::endl;
    */
    // Special handling for set scheduler. 
//     if (opts.scheduler == "set") {
//       DISABLE_ADD_TASKS = true;
//                 
//       /* Set residuals to quick start set scheduler */
//       for (size_t i = 0;i < g.num_vertices(); ++i) {
//         if (g.in_edge_ids(i).size() > 0) {
//           edge_data * edata = 
//             g.edge_data(g.in_edge_ids(i)[0]).as_ptr<edge_data>();
//           edata->residual = 1000.0;
//         }
//       }   
//       ((gl_types::set_scheduler&)graphlab->get_scheduler()).
//         begin_schedule(threshold_limit_scheduler);
//     } 
      graphlab->get_scheduler().add_task_to_all(pagerank_function, 1000.0);   
  }
  
  // if (opts.scheduler == "linear_twophase") {
  //   std::cout << "Use twophase threshold scheduler, main threshold=" << TARGET_PRECISION << std::endl;
  //   TASK_SCHEDULING_TYPE = SUMRESIDUALS;
       
  //   /* Set higher thresholds that are 1e-4*outdegree. I.e, lower outdegree vertices will be scheduled       
  //      earlier. */
  //   for (size_t i = 0; i < g.num_vertices(); ++i) {
  //     if (g.out_edge_ids(i).size() > 0) 
  //       ((linear_twophase_scheduler&)graphlab->get_scheduler()).
  //         set_vertex_higher_threshold(i, TARGET_PRECISION*g.out_edge_ids(i).size());
  //   }
  //   ((linear_twophase_scheduler&) graphlab->get_scheduler()).set_threshold(TARGET_PRECISION);

  //   std::cout << "Succesfully set vertex-specific lower thresholds" << std::endl;
  // }  

   
  // if (opts.scheduler == "linear_threshold") {
  //   std::cout << "Use linear threshold scheduler = sum residuals" << std::endl;
  //   assert(TASK_SCHEDULING_TYPE == RANDOM_SUMRESIDUALS);
  //   TASK_SCHEDULING_TYPE = SUMRESIDUALS;
  //   ((linear_threshold_scheduler&) graphlab->get_scheduler()).set_threshold(TARGET_PRECISION);

  // }   
  
  printf("Add tasks: %lf\n", t2.current_time());

  timer t;
  t.start();
  
  /* Start */
  exec_status status = graphlab->start();
  float runtime = t.current_time();
  
  printf("Finished in %lf (status: %d)\n", runtime, (int)status);
        

  /* Output in Matlab copypasteable format */
  for(vertex_id_t vid=0; vid<g.num_vertices(); vid+=100000) {
    vertex_data  vdata = g.vertex_data(vid);
    printf("%d, log:%lf\n ", vid, log(vdata.value));
  }

  printf("Runtime %f\n", runtime);
  
  /* Output l1 residual for benchmark */
  write_benchmark_value("residual", l1_residual(g)); 
  write_benchmark_value("memory_writes_mb", memory_writes * 1.0 / 1024.0/1024.0); 
  write_benchmark_value("memory_reads_mb", memory_reads * 1.0 / 1024.0/1024.0); 

}



/**** SET SCHEDULER STUFF ****/
bool compute_vertex_priority(graphlab::vertex_id_t v,
                             gl_types::iscope& scope,
                             bool& reschedule) {
  // Get the vertex data
  float sumresidual = 0.0;
  foreach(graphlab::edge_id_t edgeid, scope.in_edge_ids()) {
    //Sum residuals of each edge
    const edge_data edata =
      scope.edge_data(edgeid);
    sumresidual =+ edata.residual;
  }

  return (sumresidual >= TARGET_PRECISION);
}

// void threshold_limit_scheduler(gl_types::set_scheduler &sched) {
//   // All sets must be created before scheduling calls
//   printf("==== Set Scheduler =====\n");
//   gl_types::ivertex_set &vs1 = sched.attach(gl_types::rvset(compute_vertex_priority),sched.root_set());
//   sched.init();
//   sched.execute_rep(vs1, pagerank_function);
// }

/**** ANALYZER FUNCTIONS ****/

/* Writes current graph L1 norm and max value */
std::vector<double> pagerankdumper(gl_types::graph &g) {
  double l1res = l1_residual(g);
  std::vector<double> v;
  v.push_back(l1res);
  return v;
}


/* Used with analyzer_listener */
global_dumper pagerankapp::dump_function() {
 // return pagerankdumper;
 return NULL;
}

std::vector<std::string> pagerankapp::dump_headers() {
  std::vector<std::string> h;
  h.push_back("l1_residual");
  return h;
}

int pagerankapp::dump_frequency() {
  return 200000; // Every 200000 updates
}


