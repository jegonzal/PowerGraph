/**
 * Run parallel junction tree gibbs sampling on a factorized model
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>


#include <graphlab.hpp>


#include "image.hpp"
#include "data_structures.hpp"
#include "sequential_jt_gibbs.hpp"



#include "jt_worker.hpp"





#include <graphlab/macros_def.hpp>


#define DRAW_IMAGE


std::string results_fn = "experiment_results.tsv";


size_t get_next_experiment_id(const std::string& experiment_file) {
  std::ifstream fin(experiment_file.c_str());
  size_t lines = 0;
  std::string line;
  while(getline(fin, line)) lines++;
  fin.close();
  return lines;
}

void gibbs_sweep(mrf::graph_type& mrf, 
		 const factorized_model& factor_graph) {
  factor_t conditional;
  factor_t belief;
  for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {
    mrf::vertex_data& vdata = mrf.vertex_data(vid);
    belief.set_args(vdata.variable);
    belief.uniform();
    conditional.set_args(vdata.variable);
    const std::set<vertex_id_t>& factor_ids = 
      factor_graph.factor_ids(vdata.variable);
    foreach(vertex_id_t fid, factor_ids) {
      // get the factor
      const factor_t& factor = factor_graph.factors()[fid];
      // build the conditional
      assignment_t conditional_asg = factor.args() - vdata.variable;
      for(size_t i = 0; i < conditional_asg.num_vars(); ++i) {
	const mrf::vertex_data& other_vdata = 
	  mrf.vertex_data(conditional_asg.args().var(i).id);
	assert(conditional_asg.args().var(i) == other_vdata.variable);
	conditional_asg &= 
	  assignment_t(other_vdata.variable, other_vdata.asg);
      }
      conditional.set_args(vdata.variable);
      conditional.condition(factor, conditional_asg);
      belief *= conditional;
    }
    belief.normalize();
    vdata.asg = belief.sample().asg_at(0);
  }

}






int main(int argc, char** argv) {
  // set the global logger
  // global_logger().set_log_level(LOG_WARNING);
  // global_logger().set_log_to_console(true);


  std::cout << "This program runs junction tree blocked MCMC "
            << "inference on large factorized models."
            << std::endl;

  std::srand ( graphlab::timer::usec_of_day() );
  graphlab::random::seed();
  // std::srand ( 123  );
  // graphlab::random::seed( 123 );



  std::string model_filename = "";

  size_t treesize = 1000;
  size_t treeheight = 0;
  bool priorities = false;
  std::vector<float> runtimes(1,10);
  size_t treewidth = 3;
  size_t factorsize = 0;
  size_t subthreads = 1; 
  bool draw = false;


  // Command line parsing
  graphlab::command_line_options clopts("Parallel Junction Tree MCMC");
  clopts.attach_option("model", 
                       &model_filename, model_filename,
                       "Alchemy formatted model file");
  clopts.add_positional("model");

  clopts.attach_option("runtime", 
                       &runtimes, runtimes,
                       "total runtime in seconds");

  clopts.attach_option("treesize", 
                       &treesize, treesize,
                       "The number of variables in a junction tree");

  clopts.attach_option("treeheight", 
                       &treeheight, treeheight,
                       "The height of the tree.");



  clopts.attach_option("treewidth", 
                       &treewidth, treewidth,
                       "The maximum treewidth");

  clopts.attach_option("factorsize", 
                       &factorsize, factorsize,
                       "The maximum factorsize");

  clopts.attach_option("subthreads", 
                       &subthreads, subthreads,
                       "The number of threads to use inside each tree");


  clopts.attach_option("priorities",
                       &priorities, priorities,
                       "Use priorities?");

  clopts.attach_option("draw",
                       &draw, draw,
                       "Draw a picture?");




  clopts.scheduler_type = "fifo";
  clopts.scope_type = "edge";
  if( !clopts.parse(argc, argv) ) { 
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    return EXIT_FAILURE;
  }


  std::cout << "Load alchemy file." << std::endl;
  factorized_model factor_graph;
  factor_graph.load_alchemy(model_filename);

  std::cout << "Building graphlab MRF." << std::endl;
  mrf::graph_type mrf_graph;
  construct_mrf(factor_graph, mrf_graph);
  std::cout << "Quick init sweep" << std::endl;
  //  gibbs_sweep(mrf_graph, factor_graph);
  

  parallel_sampler sampler(factor_graph,
                           mrf_graph,
                           clopts,
                           treesize,
                           treewidth,
                           factorsize,
                           treeheight,
                           subthreads,
                           priorities);

  float run_so_far = 0;
  foreach(float runtime, runtimes) {
    
    // Get the experiment id
    size_t experiment_id = get_next_experiment_id(results_fn);

    std::cout << "Settings: ======================" << std::endl
              << "Experiment:    " << experiment_id << std::endl
              << "Model:         " << model_filename << std::endl
              << "runtime:       " << runtime << std::endl
              << "treesize:      " << treesize << std::endl
              << "treewidth:     " << treewidth << std::endl
              << "treeheight:    " << treeheight << std::endl
              << "factorsize:    " << factorsize << std::endl
              << "subthreads:    " << subthreads << std::endl
              << "priorities:    " << priorities << std::endl;   
    clopts.print();

    // run the fully parallel sampler
    float remaining_time = runtime - run_so_far;
    if(remaining_time <= 0) remaining_time = 0;

    graphlab::timer timer;
    timer.start();
    sampler.sample_once(remaining_time);
    double actual_runtime = timer.current_time();
    std::cout << "Local Runtime: " << actual_runtime << std::endl;
    
    run_so_far += actual_runtime;
    std::cout << "Total Runtime: " << run_so_far << std::endl;
    

    double loglik = unnormalized_loglikelihood(mrf_graph,
                                               factor_graph.factors());
    
    std::cout << "LogLikelihood: " << loglik << std::endl;
    mrf::save_beliefs(mrf_graph,  
                      make_filename("beliefs",".tsv", experiment_id).c_str());

    std::cout << "Total Samples: " << sampler.total_samples() 
	      << std::endl;
    std::cout << "Total Changes: " << sampler.total_changes() 
	      << std::endl;
    std::cout << "Total Trees: " << sampler.total_trees() 
	      << std::endl;
    std::cout << "Total Collisions: " << sampler.total_collisions() 
	      << std::endl;

    size_t min_samples = 0;
    size_t max_samples = 0;
    min_max_samples(mrf_graph, min_samples, max_samples);
    std::cout << "Min Samples:  " << min_samples << std::endl
              << "Max Samples:  " << max_samples << std::endl;


    // check mrf graph
    for(size_t i = 0; i < mrf_graph.num_vertices(); ++i) {
      assert(!mrf_graph.vertex_data(i).in_tree);
      assert(mrf_graph.vertex_data(i).tree_id == NULL_VID);
    }


    std::ofstream fout(results_fn.c_str(),  std::ios::app);
    fout.precision(16);
    fout << experiment_id << '\t'
         << clopts.ncpus << '\t'
         << run_so_far << '\t'
         << runtime << '\t'
         << treesize << '\t'
         << treewidth << '\t'
         << factorsize << '\t'
         << treeheight << '\t'
         << subthreads << '\t'
         << priorities << '\t'
         << actual_runtime << '\t'
         << sampler.total_samples() << '\t'
         << sampler.total_changes() << '\t'
         << sampler.total_trees() << '\t'
         << loglik << std::endl;
    fout.close();



    if(draw) {
      // Plot the final answer
      size_t rows = std::sqrt(mrf_graph.num_vertices());
      std::cout << "Rows: " << rows << std::endl;
      image img(rows, rows);
      std::vector<double> values(1);
      factor_t belief;
      for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {
        const mrf::vertex_data& vdata = mrf_graph.vertex_data(vid);
        belief = vdata.belief;
        belief.normalize();
        belief.expectation(values);
        img.pixel(vid) = values[0];
      }
      img.pixel(0) = 0;
      img.pixel(1) = mrf_graph.vertex_data(0).variable.arity-1;
      img.save(make_filename("pred", ".pgm", experiment_id).c_str());
    
      for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {   
        img.pixel(vid) = mrf_graph.vertex_data(vid).updates;
      }
      img.save(make_filename("updates", ".pgm", experiment_id).c_str());
    
      for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {   
        img.pixel(vid) = mrf_graph.vertex_data(vid).updates == 0;
      }
      img.save(make_filename("unsampled", ".pgm", experiment_id).c_str());
  
      for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {   
        img.pixel(vid) = mrf_graph.vertex_data(vid).asg;
      }
      img.pixel(0) = 0;
      img.pixel(1) = mrf_graph.vertex_data(0).variable.arity-1;
      img.save(make_filename("final_sample", ".pgm", experiment_id).c_str());

      // for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {   
      //   img.pixel(vid) = mrf_graph.vertex_data(vid).height;
      // }
      // img.save(make_filename("heights", ".pgm", experiment_id).c_str());

      for(vertex_id_t vid = 0; vid < mrf_graph.num_vertices(); ++vid) {   
        img.pixel(vid) = tanh(std::max(0.0, mrf_graph.vertex_data(vid).priority));
      }
      img.save(make_filename("priorities", ".pgm", experiment_id).c_str());

    }

  } // end of for loop over runtimes

  

  




  return EXIT_SUCCESS;
} // end of main













#include <graphlab/macros_undef.hpp>


