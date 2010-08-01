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

#include "distpagerankapp.hpp"
#include <graphlab/distributed/pushy_distributed_engine.hpp>
#include <graphlab/distributed/distributed_collaborative_scheduler_wrapper.hpp>

#include <graphlab/distributed/distributed_shared_data.hpp>

#include <graphlab/macros_def.hpp>

using namespace graphlab;

struct edge_data {
  float weight;
  float lastusedval;
  float srcvalue;
};

struct vertex_data {
  float value;
  float selfweight;
};

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
blob_distributed_graph * __g;


const size_t NUM_VERTICES = 1;
const size_t RANDOM_JUMP_WEIGHT = 2;

int myprocid, numofprocs;

//gl_types::imonitor * _listener(NULL);

//void threshold_limit_scheduler(gl_types::set_scheduler &sched);


void pagerank_function(gl_types::iscope &scope,
                       gl_types::icallback &scheduler,
                       gl_types::ishared_data* shared_data) {
  assert(shared_data != NULL);
  float jumpweight =
    shared_data->get_constant(RANDOM_JUMP_WEIGHT).as<float>();
  size_t numv =
    shared_data->get_constant(NUM_VERTICES).as<size_t>();
  // Get the data associated with the vertex
  vertex_data* vdata =
    scope.vertex_data().as_ptr<vertex_data>();
  // sum incoming weights
  
  size_t reads = 8;
  size_t writes = 0;
    
  
  float x = vdata->value * vdata->selfweight;

  foreach(edge_id_t eid, scope.in_edge_ids()) {
    edge_data* edata =
      scope.edge_data(eid).as_ptr<edge_data>();
    x += edata->weight * edata->srcvalue;
    
      if (scope.vertex() == 0) {
	  	printf("%ld %lf \n", eid, edata->weight);
	  }
  }
  writes += scope.in_edge_ids().size() * sizeof(float) * 2;
  reads += scope.in_edge_ids().size() * (sizeof(edge_id_t) + sizeof(edge_data));

  // compute the jumpweight
  x=jumpweight/numv + (1- jumpweight)*x;
  //std::cout << scope.vertex() << ": " << x << "\n";
  double residual = std::abs(x- vdata->value);
  vdata->value = x;
  
  if (scope.vertex() == 5) {
    std::cout << myprocid << "/"  << thread::thread_id() << " vertex: " <<  scope.vertex() << " ==> " << x << std::endl;
  }
  
  
  /* EXPERIMENTAL: Create random number to use for stochastic scheduling.
     Randomize only once for vertex for performance reasons. Have to think if it
     it is ok. */

  double randomNum = 0.0;
    
  if (TASK_SCHEDULING_TYPE == RANDOM_SUMRESIDUALS) {
    randomNum = graphlab::random::rand01();
    assert(randomNum >= 0.0);
    assert(randomNum <= 1.0);
  }

  double outedge_sum = 0.0;

  /* Send new value and schedule neighbors */
  foreach(edge_id_t eid, scope.out_edge_ids()) {
    edge_data* outedgedata = scope.edge_data(eid).as_ptr<edge_data>();
    outedgedata->srcvalue = x;
    outedge_sum += outedgedata->weight;
            if (scope.vertex() == 0) {
	  	printf("out %ld %lf \n", eid, outedgedata->weight);
	  }   
    switch(TASK_SCHEDULING_TYPE) {   
    case RANDOM_SUMRESIDUALS:    
      /* Schedule always if residual > target precision, and schedule
         randomly in proportion to the ratio to the precision. For
         example if residual is 1e-5, there is 1e-5/1e-4=1e-1 chance
         of scheduling the target. */
       
      if (residual / TARGET_PRECISION >= randomNum) {
        scheduler.
          add_task(gl_types::update_task(scope.target(eid), pagerank_function),
                   outedgedata->weight * residual);
        writes += sizeof(gl_types::update_task) + sizeof(double);
      }         
      break;
    case SUMRESIDUALS:
      scheduler.
        add_task(gl_types::update_task(scope.target(eid), pagerank_function),
                 residual);
      break;
    case DEGREERESIDUALSWITHMULTIPLE:
      {
        int indegree = __g->in_edge_ids(__g->target(eid)).size();
        float target_selfweight = 
          __g->vertex_data(__g->target(eid)).as<vertex_data>().selfweight;
        if ((1.0-jumpweight) * residual * indegree /
            (1.0 - (1-jumpweight) * target_selfweight) >= TARGET_PRECISION) {
          scheduler.add_task(gl_types::update_task(scope.target(eid),
                                                   pagerank_function),
                             outedgedata->weight *  residual);
          writes += sizeof(gl_types::update_task) + sizeof(double);
        }
      }
      break;
    case DEGREERESIDUALSNOMULTIPLE:
      if (residual > TARGET_PRECISION) {
        scheduler.add_task(gl_types::update_task(scope.target(eid),
                                                 pagerank_function),
                           residual);
        writes += sizeof(gl_types::update_task) + sizeof(double);
      }
      break;
    }
  }
  writes += scope.out_edge_ids().size() * sizeof(float) * 2;
  reads += scope.out_edge_ids().size() * sizeof(float) * 2;

  if (outedge_sum>1.00001) {
    printf("Edge weights too large: %lf (%d, %d)\n", outedge_sum, scope.vertex(), myprocid);
	 assert(false);
  }
  
//  if (_listener != NULL) 
 //   _listener->app_set_vertex_value(scope.vertex(), x);
    
  __sync_add_and_fetch(&memory_reads, reads);
  __sync_add_and_fetch(&memory_writes, writes);    
} 




/**** GLOBAL CONVERGENCE CALCULATOR ****/
// Calculates |Ax-x|_1
double l1_residual(blob_distributed_graph &distgraph) {
    return 0;
}




void pagerankapp::loadGraph(std::string filename, blob_graph &mgraph) {  
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
      while (mgraph.num_vertices() <= srcv || mgraph.num_vertices() <= destv) {
        vertex_data v;
        v.value = 1.0;
        v.selfweight = 0.0;
        mgraph.add_vertex( blob(sizeof(v), &v));
      }
      // add the edge. BUT
      // The graph has multiedges, AND selfedges.. this is unbelievably annoying
      if (srcv == destv) {
        // selfedge
        vertex_data *v = mgraph.vertex_data(srcv).as_ptr<vertex_data>();
        // we use selfweight temporarily to store the number of self edges
        v->selfweight += 1.0;
      }
      else {
        // check if the edge already exists
        std::pair<bool, edge_id_t> edgefind = mgraph.find(srcv, destv);

        if (edgefind.first) {
          edge_data* edata = mgraph.edge_data(edgefind.first).as_ptr<edge_data>();
 	      edata->weight += 1.0;
        }
        else {
          edge_data e;
          e.srcvalue = 1.0;
          // we use weight temporarily as a counter
          e.weight = 1.0; 
          mgraph.add_edge(srcv, destv, blob(sizeof(e), &e));
        }
      }
      ++edgecount;
      if (edgecount % 1000000 == 0) {
        std::cout << edgecount << " Edges inserted" << std::endl;
      }
    }
    

    
    // now we have to go through the edges again and figure out the
    // weights
    for (vertex_id_t i = 0;i < mgraph.num_vertices(); ++i) {
      // Check only my vertices
      
      vertex_data *v = mgraph.vertex_data(i).as_ptr<vertex_data>();
      // count the effective number of out edges
      double weight = 0;
      foreach(edge_id_t e, mgraph.out_edge_ids(i)){
        edge_data* edata = mgraph.edge_data(e).as_ptr<edge_data>();
        weight += edata->weight;
      }
      // remember to count selfweights
      weight += v->selfweight;
      if (weight == 0) continue;
      // now the weights should all be scaled to 1
      foreach(edge_id_t e, mgraph.out_edge_ids(i)) {
        edge_data* edata = mgraph.edge_data(e).as_ptr<edge_data>();
        edata->weight /= weight;
      }
      v->selfweight /= weight;
 //     printf("%ld => %lf\n", i, weight); 
      
    }
    mgraph.finalize();
  } else {
    mgraph.load(filename);
  }
}



pagerankapp::pagerankapp(graphlab::distributed_control * _dc, std::string ifile, std::string bfile, bool optimize) :
		distgraph(*_dc){
  inputfile = ifile;
  binoutfile = bfile;
  dc = _dc;
  graphoptimize=optimize;
  
  myprocid = dc->procid();
  numofprocs = dc->numprocs();
  
  std::cout << "Distributed pagerank starting. Procid = " << myprocid << "/" << numofprocs << std::endl;
}

pagerankapp::~pagerankapp() {
}


void pagerankapp::start() {
  dc->barrier();
  
  if (opts.extra == "prepartition") {
  	   // Create partitions
  	   ASSERT_EQ(dc->numprocs(), 1);
	   blob_graph mgraph;
	   loadGraph(inputfile, mgraph);
  	   blob_distributed_graph::partition_graph_tofile(mgraph, opts.ncpus, partition_method::PARTITION_METIS, inputfile);
  	   exit(0);
  } else {
  	   // Load partitions
  	   distgraph.load(inputfile, *dc);
  }
  dc->barrier();

  
  std::cout << "Graph has " << distgraph.num_vertices() << " vertices and "
            << " (? edges)" <<std::endl;
  const size_t NUM_VERTICES = 1;
  const size_t RANDOM_JUMP_WEIGHT = 2;
  
  size_t numv = distgraph.num_vertices();

  // Set register shared data
  distributed_shared_data<blob_distributed_graph> sdm(*dc);
  dc->barrier();

	if (myprocid == 0) {
 		sdm.set_constant(NUM_VERTICES, numv);
 		sdm.set_constant(RANDOM_JUMP_WEIGHT, jumpweight);
	}
	dc->barrier();

  
  // Create engine
  pushy_distributed_engine<blob_distributed_graph,
      distributed_collaborative_scheduler_wrapper<blob_distributed_graph, fifo_scheduler<blob_distributed_graph> >,
      		  general_scope_factory<blob_distributed_graph> >
      		graphlab(*dc, distgraph, opts.ncpus);
  graphlab.set_default_scope(scope_range::VERTEX_CONSISTENCY);

   dc->barrier();
   		
   printf("Created distributed engine...\n");
      		
  
  // Hack
  __g = &distgraph;
  
  if (opts.extra == "indegree_multiple") {
    printf(" ==== Using indegree multiple in residual computation ====\n");
    TASK_SCHEDULING_TYPE = DEGREERESIDUALSWITHMULTIPLE;
  } else if (opts.extra == "simple_residual") {
    printf(" ==== Using simple degree residuals ====\n");
    TASK_SCHEDULING_TYPE = DEGREERESIDUALSNOMULTIPLE;
  }
  
   printf("Creating shared data...\n");

  // register shared data
  graphlab.set_shared_data_manager(&sdm);
  
  //_listener = this->listener;
  
  
  timer t2;
  t2.start();

  dc->barrier();


   graphlab.get_scheduler().set_option(scheduler_options::UPDATE_FUNCTION,
                                       (void*)pagerank_function);
   if (myprocid == 0) {
     graphlab.get_scheduler().add_task_to_all( pagerank_function, 1000.0);   
      printf("Add tasks: %lf\n", t2.current_time());

   }


  timer t;
  t.start();
  
  
  for(vertex_id_t vid=0; vid<distgraph.num_vertices(); vid+=1) {
  	  double w = 0;
  	  if (distgraph.owner(vid) == myprocid) {
	  foreach(edge_id_t e, distgraph.out_edge_ids(vid)) {
		 edge_data * ed = distgraph.edge_data(e).as_ptr<edge_data>();
		 w += ed->weight;	
		 
		 if (vid == 0) printf("%ld %lf\n", (long int) e, ed->weight);
	  }
	  }
  }
  printf("Starting graphlab...\n");
  
  dc->barrier();
  
  /* Start */
  exec_status status = graphlab.start();
  float runtime = t.current_time();
  
  dc->barrier();
  
  printf("Finished in %lf (status: %d)\n", runtime, (int)status);
        

  for(vertex_id_t vid=0; vid<distgraph.num_vertices(); vid+=1000) {
    vertex_data * vdata = distgraph.vertex_data(vid).as_ptr<vertex_data>();
    printf("%d, log:%lf\n ", vid, log(vdata->value));
  }

  printf("Runtime %f\n", runtime);
  
  /* Output l1 residual for benchmark */
  //write_benchmark_value("residual", l1_residual(distgraph)); 
  write_benchmark_value("memory_writes_mb", memory_writes * 1.0 / 1024.0/1024.0); 
  write_benchmark_value("memory_reads_mb", memory_reads * 1.0 / 1024.0/1024.0); 

}

 

// void threshold_limit_scheduler(gl_types::set_scheduler &sched) {
//   // All sets must be created before scheduling calls
//   printf("==== Set Scheduler =====\n");
//   gl_types::ivertex_set &vs1 = sched.attach(gl_types::rvset(compute_vertex_priority),sched.root_set());
//   sched.init();
//   sched.execute_rep(vs1, pagerank_function);
// }

/**** ANALYZER FUNCTIONS ****/



