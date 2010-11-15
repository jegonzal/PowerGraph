/**
 *
 * Parallel blocked gibbs using graphlab
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>

// Including Standard Libraries
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <set> 
#include <algorithm>
#include <limits>
#include <cmath>



#include <boost/program_options.hpp>
#include <boost/bind.hpp>



#include <graphlab.hpp>



// Image reading/writing code
#include "image.hpp"
#include "util.hpp"
#include "data_structures.hpp"
#include "update_functions.hpp"
#include "image_denoise.hpp"
#include "drawing.hpp"

// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>

// #define DRAW_IMAGE

bool draw = false;




// Command Line Parsing =======================================================>


std::string tree_results_fn = "tree_results.tsv";
std::string async_results_fn = "async_results.tsv";
std::string colored_results_fn = "colored_results.tsv";


size_t get_next_experiment_id(const std::string& experiment_file) {
  std::ifstream fin(experiment_file.c_str());
  size_t lines = 0;
  std::string line;
  while(getline(fin, line)) lines++;
  return lines;
}




void run_colored_samples(const factorized_model& model,
                         gl::core& core, 
                         const std::vector<size_t>& nsamples) {
  // set some of the core options
  core.set_scheduler_type("colored");
  core.set_scope_type("null");
  core.engine().enable_sched_yield(false);

  // setup the update function 
  const bool use_callback = false;
  gl::update_function update_function = single_sample_update<use_callback>;
  core.scheduler().set_option(gl::scheduler_options::UPDATE_FUNCTION, 
                              (void*) update_function);



  // Precolor the graph
  std::cout << "Computing coloring " << std::endl;
  size_t colors = core.graph().compute_coloring();
  assert(core.graph().valid_coloring());
  std::cout << "Colors: " << colors << std::endl;

  std::cout << "Adding all factors to shared data manager" << std::endl;
  add_factors_to_sdm(core.shared_data(), model.factors());
  
  core.shared_data().set_sync(NSAMPLES_ID,
                              gl::sync_ops::sum< size_t, get_nsamples >,
                              gl::apply_ops::identity_print< size_t >,
                              size_t(0),
                              1000);
  core.engine().add_terminator(nsamples_terminator);

  graphlab::timer time;
  double runtime = 0;
  // Run the experiments
  for(size_t i = 0; i < nsamples.size(); ++i) {
    size_t samples = nsamples[i] * core.graph().num_vertices();
    core.shared_data().set_constant(MAX_NSAMPLES_ID,  
                                    graphlab::any( samples ) );

    size_t experiment_id = get_next_experiment_id(colored_results_fn);
    std::cout << "Running Colored sampler experiment: " 
              << experiment_id << std::endl;
    std::cout << "With samples: " << samples << std::endl;

    // Run the engine
    time.start();
    core.start();
    runtime += time.current_time();

    // Get the number of samples
    core.shared_data().sync(core.graph(), NSAMPLES_ID);
    size_t actual_samples = core.shared_data().get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = core.engine().last_update_count();

    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Samples requested: " << samples << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;
    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(core.graph(),
                 make_filename("colored_blfs_", ".tsv", experiment_id));


    std::cout << "Saving Assignment." << std::endl;
    save_asg(core.graph(),
             make_filename("colored_asg_", ".asg", experiment_id));



    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(core.graph(), model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;
    
    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(core.graph(), min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    
    std::ofstream fout(colored_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << actual_samples << '\t'
         << core.get_engine_options().ncpus << '\t'
	 << colors << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();




    if(draw) {
    // Plot the final answer
    size_t rows = std::sqrt(core.graph().num_vertices());
    std::cout << "Rows: " << rows << std::endl;
    image img(rows, rows);
    std::vector<double> values(1);

    factor_t belief;
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
      const vertex_data& vdata = core.graph().vertex_data(vid);
      belief = vdata.belief;
      belief.normalize();
      belief.expectation(values);
      img.pixel(vid) = values[0];
    }
    img.pixel(0) = 0;
    img.pixel(1) = core.graph().vertex_data(0).variable.arity-1;
    img.save(make_filename("pred", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates;
    }
    img.save(make_filename("updates", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates == 0;
    }
    img.save(make_filename("unsampled", ".pgm", experiment_id).c_str());
  
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).asg;
    }
    img.save(make_filename("final_sample", ".pgm", experiment_id).c_str());
    }




  } // end of for loop 
  std::cout << "Finished." << std::endl;
}









void run_colored_times(const factorized_model& model,
                       gl::core& core,
                       const std::vector<size_t>& times) {


  // set some of the core options
  core.set_scheduler_type("colored");
  core.set_scope_type("null");
  core.engine().enable_sched_yield(false);

  // setup the update function 
  const bool use_callback = false;
  gl::update_function update_function = single_sample_update<use_callback>;
  core.scheduler().set_option(gl::scheduler_options::UPDATE_FUNCTION, 
                              (void*) update_function);




  // Precolor the graph
  std::cout << "Computing coloring " << std::endl;
  size_t colors = core.graph().compute_coloring();
  std::cout << "Colors: " << colors << std::endl;


  std::cout << "Adding all factors to shared data manager" << std::endl;
  add_factors_to_sdm(core.shared_data(), model.factors());
  
  core.shared_data().set_sync(NSAMPLES_ID,
                              gl::sync_ops::sum< size_t, get_nsamples >,
                              gl::apply_ops::identity_print< size_t >,
                              size_t(0),
                              -1);

  graphlab::timer timer;
  double runtime = 0;
  double last_time = 0;
  // Run the experiments
  for(size_t i = 0; i < times.size(); ++i) {
    size_t timeout = times[i] - last_time;
    last_time = times[i];
    if(timeout == 0) continue;
    assert(timeout <= times[i]);
    std::cout << "Timeout: " << timeout << std::endl;
    core.engine().set_timeout(timeout);

    size_t experiment_id = get_next_experiment_id(colored_results_fn);
    std::cout << "Running Colored sampler Time experiment:" 
              << experiment_id << std::endl;
    std::cout << "With runtime: " << times[i] << std::endl;

    // Run the engine
    timer.start();
    core.start();
    runtime += timer.current_time();

    // Get the number of samples
    core.shared_data().sync(core.graph(), NSAMPLES_ID);
    size_t actual_samples = 
      core.shared_data().get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = core.engine().last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(core.graph(),
                 make_filename("colored_blfs_", ".tsv", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(core.graph(),
             make_filename("colored_asg_", ".asg", experiment_id));


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(core.graph(), model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(core.graph(), min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;


    std::ofstream fout(colored_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << actual_samples << '\t'
         << core.get_engine_options().ncpus << '\t'
	 << colors << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();




    if(draw) {
    // Plot the final answer
    size_t rows = std::sqrt(core.graph().num_vertices());
    std::cout << "Rows: " << rows << std::endl;
    image img(rows, rows);
    std::vector<double> values(1);
    factor_t belief;
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
      const vertex_data& vdata = core.graph().vertex_data(vid);
      belief = vdata.belief;
      belief.normalize();
      belief.expectation(values);
      img.pixel(vid) = values[0];
    }
    img.pixel(0) = 0;
    img.pixel(1) = core.graph().vertex_data(0).variable.arity-1;
    img.save(make_filename("pred", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates;
    }
    img.save(make_filename("updates", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates == 0;
    }
    img.save(make_filename("unsampled", ".pgm", experiment_id).c_str());
  
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).asg;
    }
    img.pixel(0) = 0;
    img.pixel(1) = core.graph().vertex_data(0).variable.arity-1;
    img.save(make_filename("final_sample", ".pgm", experiment_id).c_str());
    }



  } // end of for loop 
  std::cout << "Finished." << std::endl;
}













void run_async_samples(const factorized_model& model,
                       gl::core& core, 
                       const std::vector<size_t>& nsamples) {

  add_factors_to_sdm(core.shared_data(), model.factors());

  
  core.shared_data().set_sync(NSAMPLES_ID,
                              gl::sync_ops::sum< size_t, get_nsamples >,
                              gl::apply_ops::identity_print< size_t >,
                              size_t(0),
                              1000);
  core.engine().add_terminator(nsamples_terminator);


  // Set the initial schedule
  std::vector<vertex_id_t> vertices(core.graph().num_vertices(), 0);
  for(vertex_id_t i = 0; i < vertices.size(); ++i) vertices[i] = i;
  std::random_shuffle(vertices.begin(), vertices.end());
  const double residual = 1.0;
  const bool UseCallback = true;
  for(size_t i = 0; i < vertices.size(); ++i) {
    gl::update_task task(vertices[i], single_sample_update<UseCallback>);
    core.add_task(task, residual);
  }

  graphlab::timer time;
  double runtime = 0;
  // Run the experiments
  for(size_t i = 0; i < nsamples.size(); ++i) {
    size_t samples = nsamples[i] * core.graph().num_vertices();
    core.shared_data().set_constant(MAX_NSAMPLES_ID,  graphlab::any( samples ) );

    size_t experiment_id = get_next_experiment_id(async_results_fn);
    std::cout << "Running Asynchronous sampler experiment: " 
              << experiment_id << std::endl;
    std::cout << "With samples: " << samples << std::endl;

    // Run the engine
    time.start();
    core.start();
    runtime += time.current_time();

    // Get the number of samples
    core.shared_data().sync(core.graph(), NSAMPLES_ID);
    size_t actual_samples = core.shared_data().get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = core.engine().last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Samples requested: " << samples << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(core.graph(),
                 make_filename("async_blfs_", ".tsv", experiment_id));


    std::cout << "Saving Assignment." << std::endl;
    save_asg(core.graph(),
             make_filename("async_asg_", ".asg", experiment_id));


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(core.graph(), model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(core.graph(), min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;
    
    

    std::ofstream fout(async_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << core.get_engine_options().scheduler_type << '\t'
         << actual_samples << '\t'
         << core.get_engine_options().ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();







    if(draw) {
    // Plot the final answer
    size_t rows = std::sqrt(core.graph().num_vertices());
    std::cout << "Rows: " << rows << std::endl;
    image img(rows, rows);
    std::vector<double> values(1);
    factor_t belief;
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
      const vertex_data& vdata = core.graph().vertex_data(vid);
      belief = vdata.belief;
      belief.normalize();
      belief.expectation(values);
      img.pixel(vid) = values[0];
    }
    img.pixel(0) = 0;
    img.pixel(1) = core.graph().vertex_data(0).variable.arity-1;
    img.save(make_filename("pred", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates;
    }
    img.save(make_filename("updates", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates == 0;
    }
    img.save(make_filename("unsampled", ".pgm", experiment_id).c_str());
  
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).asg;
    }
    img.save(make_filename("final_sample", ".pgm", experiment_id).c_str());
    }




  } // end of for loop 
  std::cout << "Finished." << std::endl;
}













void run_async_times(const factorized_model& model,
                     gl::core& core, 
                     const std::vector<size_t>& times) {
  add_factors_to_sdm(core.shared_data(), model.factors());  
  core.shared_data().set_sync(NSAMPLES_ID,
                              gl::sync_ops::sum< size_t, get_nsamples >,
                              gl::apply_ops::identity_print< size_t >,
                              size_t(0),
                              -1);
 

 
  // Set the initial schedule
  std::vector<vertex_id_t> vertices(core.graph().num_vertices(), 0);
  for(vertex_id_t i = 0; i < vertices.size(); ++i) vertices[i] = i;
  std::random_shuffle(vertices.begin(), vertices.end());
  const double residual = 1.0;
  const bool UseCallback = true;
  for(size_t i = 0; i < vertices.size(); ++i) {
    gl::update_task task(vertices[i], single_sample_update<UseCallback>);
    core.add_task(task, residual);
  }

  graphlab::timer timer;
  double runtime = 0;
  double last_time = 0;
  // Run the experiments
  for(size_t i = 0; i < times.size(); ++i) {
    size_t timeout = times[i] - last_time;
    last_time = times[i];
    if(timeout == 0) continue;
    assert(timeout <= times[i]);
    std::cout << "Timeout: " << timeout << std::endl;
    core.engine().set_timeout(timeout);

    size_t experiment_id = get_next_experiment_id(async_results_fn);
    std::cout << "Running Asynchronous sampler Time experiment:" 
              << experiment_id << std::endl;
    std::cout << "With runtime: " << times[i] << std::endl;

    // Run the engine
    timer.start();
    core.start();
    runtime += timer.current_time();

    // Get the number of samples
    core.shared_data().sync(core.graph(), NSAMPLES_ID);
    size_t actual_samples = core.shared_data().get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = core.engine().last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(core.graph(),
                 make_filename("async_blfs_", ".tsv", experiment_id));

    
    std::cout << "Saving Assignment." << std::endl;
    save_asg(core.graph(),
             make_filename("async_asg_", ".asg", experiment_id));


    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(core.graph(), model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;
    
    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(core.graph(), min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    
    std::ofstream fout(async_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << core.get_engine_options().scheduler_type << '\t'
         << actual_samples << '\t'
         << core.get_engine_options().ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << std::endl;
    fout.close();






    if(draw) {
    // Plot the final answer
    size_t rows = std::sqrt(core.graph().num_vertices());
    std::cout << "Rows: " << rows << std::endl;
    image img(rows, rows);
    std::vector<double> values(1);
    factor_t belief;
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
      const vertex_data& vdata = core.graph().vertex_data(vid);
      belief = vdata.belief;
      belief.normalize();
      belief.expectation(values);
      img.pixel(vid) = values[0];
    }
    img.pixel(0) = 0;
    img.pixel(1) = core.graph().vertex_data(0).variable.arity-1;
    img.save(make_filename("pred", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates;
    }
    img.save(make_filename("updates", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates == 0;
    }
    img.save(make_filename("unsampled", ".pgm", experiment_id).c_str());
  
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).asg;
    }
    img.save(make_filename("final_sample", ".pgm", experiment_id).c_str());
    }






  } // end of for loop 
  std::cout << "Finished." << std::endl;
}






















void run_tree_samples(const factorized_model& model,
                      gl::core& core,
                      std::vector<size_t> nsamples,
                      size_t nroots,
                      size_t tree_height,
                      bool use_weights,
                      double pruning) {

 
  // Setup the engine
  set_tree_sampler_constants(core.shared_data(),
                             nroots,
                             tree_height,
                             use_weights,
                             pruning);
  global_graph = &core.graph();

  add_factors_to_sdm(core.shared_data(), model.factors());
  core.shared_data().set_sync(NSAMPLES_ID,
                              gl::sync_ops::sum< size_t, get_nsamples >,
                              gl::apply_ops::identity_print< size_t >,
                              size_t(0),
                              1000);
  core.engine().add_terminator(nsamples_terminator);


  // Create the initial set of tasks
  for(size_t j = 0; j < nroots; ++j) {
    // Add the root
    vertex_id_t root = get_next_root(j, nroots, core.graph().num_vertices());
    // Create a task to grow a tree from the root
    gl::update_task task(root, grow_root_update);
    core.add_task(task, grow_root_residual);
  }

  graphlab::timer time;
  double runtime = 0;
  // Run the experiments
  for(size_t i = 0; i < nsamples.size(); ++i) {
    size_t samples = nsamples[i] * core.graph().num_vertices();
    core.shared_data().set_constant(MAX_NSAMPLES_ID,  graphlab::any( samples ) );

    size_t experiment_id = get_next_experiment_id(tree_results_fn);
    std::cout << "Running Tree sampler Samples experiment: " 
              << experiment_id << std::endl;
    std::cout << "With samples: " << samples << std::endl;

    // Run the engine
    time.start();
    core.start();
    runtime += time.current_time();

    // Get the number of samples
    core.shared_data().sync(core.graph(), NSAMPLES_ID);
    size_t actual_samples = core.shared_data().get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = core.engine().last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Samples requested: " << samples << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(core.graph(),
                 make_filename("tree_blfs_", ".tsv", experiment_id));

    std::cout << "Saving Assignment." << std::endl;
    save_asg(core.graph(),
             make_filename("tree_asg_", ".asg", experiment_id));    

    
    std::cout << "Saving tree state." << std::endl;
    save_tree_state(core.graph(),
                    make_filename("tree_state_", ".state", experiment_id));    

    

    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(core.graph(), model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(core.graph(), min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    
    std::ofstream fout(tree_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << core.get_engine_options().scheduler_type << '\t'
         << actual_samples << '\t'
         << core.get_engine_options().ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << nroots << '\t'       
         << use_weights << '\t'
         << pruning << '\t'
         << tree_height << std::endl;
    fout.close();






    if(draw) {
    // Plot the final answer
    size_t rows = std::sqrt(core.graph().num_vertices());
    std::cout << "Rows: " << rows << std::endl;
    image img(rows, rows);
    std::vector<double> values(1);
    factor_t belief;
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
      const vertex_data& vdata = core.graph().vertex_data(vid);
      belief = vdata.belief;
      belief.normalize();
      belief.expectation(values);
      img.pixel(vid) = values[0];
    }
    img.pixel(0) = 0;
    img.pixel(1) = core.graph().vertex_data(0).variable.arity-1;
    img.save(make_filename("pred", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates;
    }
    img.save(make_filename("updates", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates == 0;
    }
    img.save(make_filename("unsampled", ".pgm", experiment_id).c_str());
  
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).asg;
    }
    img.save(make_filename("final_sample", ".pgm", experiment_id).c_str());
    }





  } // end of for loop 
  std::cout << "Finished." << std::endl;
}










void run_tree_times(const factorized_model& model,
                    gl::core& core,
                    std::vector<size_t> times,
                    size_t nroots,
                    size_t tree_height,
                    bool use_weights,
                    double pruning) {
  
  set_tree_sampler_constants(core.shared_data(),
                             nroots,
                             tree_height,
                             use_weights,
                             pruning);

  add_factors_to_sdm(core.shared_data(), model.factors());

  global_graph = & core.graph();
  
  core.shared_data().set_sync(NSAMPLES_ID,
                              gl::sync_ops::sum< size_t, get_nsamples >,
                              gl::apply_ops::identity_print< size_t >,
                              size_t(0),
                              -1);

  // Create the initial set of tasks
  for(size_t j = 0; j < nroots; ++j) {
    // Add the root
    vertex_id_t root = get_next_root(j, nroots, core.graph().num_vertices());
    // Create a task to grow a tree from the root
    gl::update_task task(root, grow_root_update);
    core.add_task(task, grow_root_residual);
  }

  graphlab::timer timer;
  double runtime = 0;
  double last_time = 0;
  // Run the experiments
  for(size_t i = 0; i < times.size(); ++i) {
    size_t timeout = times[i] - last_time;
    last_time = times[i];
    if(timeout == 0) continue;
    assert(timeout <= times[i]);
    std::cout << "Timeout: " << timeout << std::endl;
    core.engine().set_timeout(timeout);

    size_t experiment_id = get_next_experiment_id(tree_results_fn);
    std::cout << "Running Tree sampler Time experiment:" 
              << experiment_id << std::endl;
    std::cout << "With runtime: " << times[i] << std::endl;

    // Run the engine
    timer.start();
    core.start();
    runtime += timer.current_time();

    // Get the number of samples
    core.shared_data().sync(core.graph(), NSAMPLES_ID);
    size_t actual_samples = core.shared_data().get(NSAMPLES_ID).as<size_t>();    
    
    // Get the number of updates
    size_t update_counts = core.engine().last_update_count();


    std::cout << "Runtime:           " << runtime << std::endl;
    std::cout << "Actual samples:    " << actual_samples << std::endl;    
    std::cout << "Update Counts:     " << update_counts  << std::endl;

    
    std::cout << "Saving belief." << std::endl;
    save_beliefs(core.graph(),
                 make_filename("tree_blfs_", ".tsv", experiment_id));
    
    std::cout << "Saving Assignment." << std::endl;
    save_asg(core.graph(),
             make_filename("tree_asg_", ".tsv", experiment_id));

    std::cout << "Saving tree state." << std::endl;
    save_tree_state(core.graph(),
                    make_filename("tree_state_", ".state", experiment_id));    

    
    
    std::cout << "Computing loglikelihood of final assignment." << std::endl;
    double loglik = unnormalized_likelihood(core.graph(), model.factors());
    std::cout << "Loglikelihood:    " << loglik << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(core.graph(), min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;

    
    std::ofstream fout(tree_results_fn.c_str(), std::ios::app);
    fout << experiment_id << '\t'
         << core.get_engine_options().scheduler_type << '\t'
         << actual_samples << '\t'
         << core.get_engine_options().ncpus << '\t'
         << runtime << '\t'
         << update_counts << '\t'
         << loglik << '\t'
         << nroots << '\t'       
         << use_weights << '\t'
         << pruning << '\t'
         << tree_height << std::endl;
    fout.close();


    if(draw) {
    // Plot the final answer
    size_t rows = std::sqrt(core.graph().num_vertices());
    std::cout << "Rows: " << rows << std::endl;
    image img(rows, rows);
    std::vector<double> values(1);
    factor_t belief;
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
      const vertex_data& vdata = core.graph().vertex_data(vid);
      belief = vdata.belief;
      belief.normalize();
      belief.expectation(values);
      img.pixel(vid) = values[0];
    }
    img.pixel(0) = 0;
    img.pixel(1) = core.graph().vertex_data(0).variable.arity-1;
    img.save(make_filename("pred", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates;
    }
    img.save(make_filename("updates", ".pgm", experiment_id).c_str());
    
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).updates == 0;
    }
    img.save(make_filename("unsampled", ".pgm", experiment_id).c_str());
  
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {   
      img.pixel(vid) = core.graph().vertex_data(vid).asg;
    }
    img.save(make_filename("final_sample", ".pgm", experiment_id).c_str());
    }





  } // end of for loop 
  std::cout << "Finished." << std::endl;
}






// MAIN =======================================================================>
int main(int argc, char** argv) { 

  // set the global logger
  global_logger().set_log_level(LOG_WARNING);
  global_logger().set_log_to_console(true);


  std::srand ( graphlab::timer::usec_of_day() );
  graphlab::random::seed();



  
  // Parse command line arguments --------------------------------------------->
  std::string model_filename;  
  size_t treeheight = 1000;
  size_t nroots = 3;
  bool priorities = false;
  double pruning  = 0;
  std::vector<size_t> nsamples;
  std::vector<size_t> runtimes(1,10);  
  std::string experiment_type = "colored";

  graphlab::command_line_options clopts;

  clopts.attach_option("model",
                       &model_filename, model_filename,
                       "model file name");
  clopts.add_positional("model");

  clopts.attach_option("treeheight", 
                       &treeheight, treeheight,
                       "The height of the tree.");

  clopts.attach_option("nroot",
                       &nroots, nroots,
                       "number of simultaneous tree roots");

  clopts.attach_option("priorities",
                       &priorities, priorities,
                       "Use weights when constructing trees");

  clopts.attach_option("pruning",
                       &pruning, pruning,
                       "Pruning threshold (== 0 implies no pruning)");

  clopts.attach_option("nsamples",
                       &nsamples, nsamples,
                       "number of samples to compute");

  clopts.attach_option("runtime", 
                       &runtimes, runtimes,
                       "total runtime in seconds");

  clopts.attach_option("draw", 
                       &draw, draw,
                       "draw pictures");



  clopts.attach_option("experiment", 
                       &experiment_type, experiment_type,
                       "the type of experiment to run "
                       "{tree, async, colored}");

  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Application Options" << std::endl;
  std::cout 
    << "model:          " << model_filename << std::endl
    << "treeheight:     " << treeheight << std::endl
    << "nroots:         " << nroots << std::endl
    << "priorities:     " << priorities << std::endl
    << "pruning:        " << pruning << std::endl
    << "nsamples:       "
    << boost::lexical_cast<std::string>(nsamples) << std::endl
    << "runtime:        "
    << boost::lexical_cast<std::string>(runtimes) << std::endl
    << "experiment:     " << experiment_type << std::endl;
  std::cout << "Graphlab Options" << std::endl;
  clopts.print();


  // create model filename
  std::cout << "Load alchemy file." << std::endl;
  factorized_model factor_graph;
  factor_graph.load_alchemy(model_filename);
  std::cout << "Building graphlab MRF." << std::endl;
  gl::core core;
  construct_clique_graph(factor_graph, core.graph());
  std::cout << "Setting core options." << std::endl;
  core.set_engine_options(clopts);



  // for(size_t i = 0; i < problem.graph.num_vertices(); ++i) {
  //   vertex_data& vdata = problem.graph.vertex_data(i);
  //   vdata.asg.uniform_sample();
  // }
  
  // Create synthetic images -------------------------------------------------->
  std::cout << "Running experiment" << std::endl;


  if(experiment_type == "tree") {
    if(runtimes.empty()) {
      run_tree_samples(factor_graph,
                       core,
                       nsamples,
                       nroots,
                       treeheight,
                       priorities,
                       pruning);
    } else {
      run_tree_times(factor_graph,
                     core,
                     runtimes,
                     nroots,
                     treeheight,
                     priorities,
                     pruning);
    }
  } else if(experiment_type == "async") {
    if(runtimes.empty()) {
      run_async_samples(factor_graph,
                        core,
                        nsamples);
    } else {
      run_async_times(factor_graph,
                      core,
                      runtimes);
                      
    }
  } else if(experiment_type == "colored") {
    if(runtimes.empty()) {
      run_colored_samples(factor_graph, 
                          core,
                          nsamples);
    } else {
      run_colored_times(factor_graph,
                        core,
                        runtimes);
    }
  } else {
    std::cout << "Invalid experiment: " << experiment_type << std::endl;
    return EXIT_FAILURE;
  }  
  



  
  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;  
} // End of main










