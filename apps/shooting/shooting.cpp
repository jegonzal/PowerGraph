/**
 * This file contains a parallel implementation of the shooting
 * algorithm using graphlab
 *
 *  \author Joseph Gonzalez, Yucheng Low
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>


#include <graphlab.hpp>


// Include the macro for the for each operation
#include <graphlab/macros_def.hpp>

enum std_sync_keys { WCHANGE_KEY };

enum sdt_const_keys {LAMBDA_KEY,
                     WBOUND_KEY,
                     NUMDIM_KEY,
                     NUMDAT_KEY};

enum vertex_type {DATA_VERTEX, FEATURE_VERTEX};

struct vertex_data {
  vertex_type type;
  // Feature members
  float weight;
  float wchange;
  float inv_xtx;

  // Data members
  float yvalue;
  float yresidual;
  
  vertex_data() :
    type(DATA_VERTEX),
    weight(0), wchange(0), inv_xtx(0),
    yvalue(0), yresidual(0) { }
};

struct edge_data {
  float x;
  edge_data() : x(0) { }
};


// Define the graphlab types
typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;


// Termination function
bool nsamples_terminator(const gl_types::ishared_data& shared_data) {
  // Get the termination bound
  const double& bound =
    shared_data.get_constant(WBOUND_KEY).as<double>();
  // Get the ammount of change
  double change = shared_data.get(WCHANGE_KEY).as<double>();
  bool terminate = (change != 0) && (change < bound);
  if(terminate) {
    std::cout << "Termination condition reached" << std::endl;
  }
  return terminate;
}




// Primary update function
void weight_update_function(gl_types::iscope& scope,
                            gl_types::icallback& scheduler,
                            gl_types::ishared_data* shared_data) {

  // get some constants
  assert(shared_data != NULL);
  const double& lambda = shared_data->get_constant(LAMBDA_KEY).as<double>();
  const double& wbound = shared_data->get_constant(WBOUND_KEY).as<double>();

  // Get the vertex data
  vertex_data& vdata = scope.vertex_data();
  assert(vdata.type == FEATURE_VERTEX);
  
  // if we do not have inverse x'x we need to compute it
  if (vdata.inv_xtx == 0)  {
    foreach(gl_types::edge_id_t eid, scope.out_edge_ids()) {        
      const edge_data& edata = scope.edge_data(eid);
      vdata.inv_xtx += (edata.x * edata.x);
    }
    vdata.inv_xtx = 1.0 / vdata.inv_xtx;
  }
  
  assert(vdata.inv_xtx != 0);
  
  // compute my new weight
  // remember to add back my old contributions
  // localcontribution = w(i) * x(:, i)
  float least_squares_estimator = 0;
  foreach(gl_types::edge_id_t eid, scope.out_edge_ids()) { 
    const edge_data& edata = scope.const_edge_data(eid);
    const vertex_data& nvdata =
      scope.neighbor_vertex_data(scope.target(eid));
    assert(nvdata.type == DATA_VERTEX);
    // here we "undo" the previous contribution restoring some of the
    // original "yvalue" 
    float local_contribution = edata.x * vdata.weight;
    least_squares_estimator +=
      (nvdata.yresidual + local_contribution) * edata.x;
  }
  least_squares_estimator = vdata.inv_xtx * least_squares_estimator;
  
  // reproject
  if (least_squares_estimator > lambda)
    least_squares_estimator -= lambda;
  else if (least_squares_estimator < -lambda)
    least_squares_estimator += lambda;
  else least_squares_estimator = 0;

  
  // Update the weight and record the weight change
  float oldw = vdata.weight;
  vdata.weight = least_squares_estimator;
  vdata.wchange = fabs(oldw - vdata.weight);


  // now subtract the local contribution from the remote vertices
  foreach(gl_types::edge_id_t eid, scope.out_edge_ids()) {
    const edge_data& edata = scope.edge_data(eid);
    gl_types::vertex_id_t targetv = scope.target(eid);
    vertex_data& nvdata = scope.neighbor_vertex_data(targetv);
    assert(nvdata.type == DATA_VERTEX);
    // Add back the old contribution
    nvdata.yresidual += edata.x * oldw;
    // Subtract out the new contribution
    nvdata.yresidual -= edata.x * vdata.weight;
  }


  // Rescheduling phase =============================================>
  if (vdata.wchange < wbound) {
    //scheduler.add_task(update_task(scope.vertex(), w_update_function), wchange);
    return;
  }

  // Compute the priorities of all the weights that are neighbors of
  // any of the data that depend on this weight
  std::map<gl_types::vertex_id_t, double> reschedule;
  // Determine the neighbor priorities
  foreach(gl_types::edge_id_t eid, scope.out_edge_ids()) {
    const edge_data& edata = scope.edge_data(eid);
    gl_types::vertex_id_t targetv = scope.target(eid);
    if (vdata.wchange * edata.x > wbound)  {
      foreach(gl_types::edge_id_t eid2, scope.in_edge_ids(targetv) ) {
        gl_types::vertex_id_t targetv2 = scope.source(eid2);
        if(targetv2 != scope.vertex())
          reschedule[targetv2] += vdata.wchange * edata.x;
      }
    }
  }

  // Actually schedule
  typedef std::pair<gl_types::vertex_id_t, double> reschedpair;  
  foreach(reschedpair pair, reschedule) {
    gl_types::update_task task(pair.first, weight_update_function);
    scheduler.add_task(task, pair.second);
  }
}



bool build_graph(graph_type& graph,
                 size_t& numdat,
                 size_t& numdim,
                 const std::string& data_file,
                 const std::string& yfile) {
  std::ifstream fin(data_file.c_str());
  if(!fin.good()) return false;
  // first pass. Count the number of data points number of features
  // and number of edges
  size_t numedges = 0;
  numdat = 0;
  numdim = 0;
  while (fin.good()) {
    size_t dat, feat;
    double val;
    fin >> dat >> feat >> val;
    numdat = std::max(numdat, dat);
    numdim = std::max(numdim, feat);
    ++numedges;
  }
  fin.close();
  // Increment the number of data points to accomodate the 0 data
  // point
  numdat++;
  // Do not need to increment the number of dimensions since 0
  // represents the y value
  // numdim++;
  
  
  std::cout << numdim << " dimensions" << std::endl
            << numdat << " data points" << std::endl
            << numedges << " edges" << std::endl;

  
  // Add the data vertices
  for (size_t i = 0;i < numdat; ++i) {
    vertex_data vdata;
    vdata.type = DATA_VERTEX;
    graph.add_vertex(vdata);
  }
  // Adding feature vertices
  for (size_t i = 0;i < numdim; ++i) {
    vertex_data vdata;
    vdata.type = FEATURE_VERTEX;
    graph.add_vertex(vdata);
  }
  
  // Read in the actual data adding edges and initializing vertex
  // values
  size_t numread = 0;
  fin.open(data_file.c_str()); 
  while (fin.good()) {
    size_t dat, feat;
    double val;
    fin >> dat >> feat >> val;
    if (fin.fail()) break;    
    numread++;
    if (numread % 100000==0) {
      std::cout << numread << "/" << numedges << "\n";
    }
    // this is class label
    if (feat == 0) {
      vertex_data& vdata = graph.vertex_data(dat);
      assert(vdata.type == DATA_VERTEX);
      vdata.yvalue = val;
      vdata.yresidual = val;
    } else {
      // add an edge from feature to the data
      edge_data edata;
      edata.x = val;
      graph.add_edge(feat + numdat - 1, dat, edata);
    }    
  }
  fin.close();

  // A separte yfile is optionally provided
  if (yfile.length() > 0) {
    size_t datid = 0;
    fin.open(yfile.c_str()); 
    while (fin.good()) {
      double yval;
      fin >> yval;
     // std::cout << yval << std::endl;
      if (fin.fail()) break;    
      vertex_data& vdata = graph.vertex_data(datid);
      assert(vdata.type == DATA_VERTEX);
      vdata.yvalue = yval;
      vdata.yresidual = yval;
      ++datid;
    }
    std::cout << "Read " << datid << " y-values" << std::endl; 
  }
  graph.finalize();
  std::cout << graph.num_vertices() << " vertices" << std::endl
            << graph.num_edges() << " edges" << std::endl;
  return true;
}



double get_wchange(const vertex_data& vdata) {
  assert(vdata.type == FEATURE_VERTEX);
  return vdata.wchange;
}




double l2error(graph_type& graph, gl_types::vertex_id_t vertexid){ 
  // compute alpha
  double sum = 0;
  foreach(gl_types::edge_id_t eid, graph.in_edge_ids(vertexid)) {        
    const edge_data& edata = graph.edge_data(eid);
    const vertex_data& vdata = graph.vertex_data(graph.source(eid));
    assert(vdata.type == FEATURE_VERTEX);
    sum += edata.x * vdata.weight;
  }
  const vertex_data& vdata = graph.vertex_data(vertexid);
  assert(vdata.type == DATA_VERTEX);
  return (sum - vdata.yvalue) * (sum - vdata.yvalue);  
}





int main(int argc, char** argv) {
  
  // sets the logging level of graphlab
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(false);

  // Parse command line arguments
  graphlab::command_line_options clopts("Shooting algorithm input arguments.");
  // Get the input files
  std::string data_file;
  std::string yfile;
  double lambda = 0.5;
  double wbound = 1E-4;
  clopts.attach_option("datafile", &data_file, "The data file in csv format.");
  clopts.attach_option("yfile", &yfile, "An optional file containing y values.");
  clopts.attach_option("lambda", &lambda, lambda, "The regularization term");
  clopts.attach_option("wbound", &wbound, wbound, "The termination bound.");
  // Set data file as the first positional (mandatory) op
  clopts.add_positional("datafile");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Failed parsing command line args" << std::endl;
    return EXIT_FAILURE;
  }

  if(!clopts.is_set("datafile")) {
    std::cout << "No datafile provided!" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }

  std::cout << "Lambda: " << lambda << "\n"
            << "Bound:  " << wbound << std::endl;
  
  // create a graphlab core
  gl_types::core core;
  // Set the engine options for the core
  core.set_engine_options(clopts);


  std::cout << "Building Graph:" << std::endl;
  size_t numdat = 0;
  size_t numdim = 0;
  bool success = 
    build_graph(core.graph(), numdat, numdim, data_file, yfile);
  if(!success) {
    std::cout << "Failed to parse file: " << data_file << std::endl;
    if(yfile.size() > 0) 
      std::cout << "\t Optional y-file: " << yfile << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Finished!" << std::endl;


  std::cout << "Configuring Shared Data:" << std::endl;
  core.shared_data().set_constant(LAMBDA_KEY, graphlab::any(lambda));
  core.shared_data().set_constant(WBOUND_KEY, graphlab::any(wbound));
  // Adding sync
  core.shared_data().
    set_sync(WCHANGE_KEY,
             gl_types::sync_ops::max< double, get_wchange >,
             gl_types::apply_ops::identity_print< double >,
             double(0),
             100,
             numdat, numdat+numdim);
  std::cout << "Finished!" << std::endl;

  
  std::cout << "Scheduling all tasks." << std::endl;
  for(size_t i = numdat; i < numdat+numdim; ++i) 
    core.add_task(i, weight_update_function, 1000.0);
  std::cout << "Finished!" << std::endl;


  std::cout << "Running the engine" << std::endl;
  double runtime = core.start();
  std::cout << "Finished! Runtime:" << runtime << std::endl;


  std::cout << "Saving results:" << std::endl;
  std::ofstream fout("weights.txt");
  double score = 0;
  for (size_t i = 0; i < numdim; ++i) {
    const vertex_data& vdata = core.graph().vertex_data(numdat+i);
    fout << vdata.weight ;
    std::cout << vdata.weight ;
    if( i + 1 < numdim ) {
      fout << ", ";
      std::cout << ", ";
    } 
    score += std::fabs(lambda * vdata.weight);
  }
  fout << std::endl;
  std::cout << std::endl;
  fout.close();
  std::cout << "regularization loss " << score << std::endl;


  double l2loss = 0;
  for (size_t i = 0;i < numdat; ++i) {
    l2loss += l2error(core.graph(), i);
  }
  std::cout << "l2 loss   " << l2loss << "\n"
            << "Taskcount " << core.last_update_count() << "\n"
            << "Runtime   " << runtime << std::endl;

 
  return EXIT_SUCCESS;
}















#include <graphlab/macros_undef.hpp>

