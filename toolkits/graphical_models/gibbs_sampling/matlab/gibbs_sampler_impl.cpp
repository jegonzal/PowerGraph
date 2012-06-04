/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


/**
 * This file contains the code used to run the chromatic gibbs sampler
 * from matlab
 */

#include "mex.h"

#include <iostream>

#include <graphlab.hpp>



#include "../factorized_model.hpp"
#include "../mrf.hpp"
#include "../global_variables.hpp"
#include "../junction_tree.hpp"

#include "../chromatic_sampler.hpp"
#include "../jt_splash_sampler.hpp"



#include "matwrap.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////
// Struct fieldnames

const char* vars_field_name = "vars";
const char* logP_field_name = "logP";
int vars_field_id = -1;
int logP_field_id = -1;





///////////////////////////////////////////////////////////////////////////


struct options {
  enum { CHROMATIC, SPLASH } alg_type;
  size_t nsamples;
  //  double tskip;
  size_t nskip;
  size_t ncpus;

  //  size_t ntrees;
  size_t treewidth;
  size_t treeheight;
  size_t treesize;
  bool priorities;
  size_t vanish;
  bool save_alchemy;
  size_t ncpus_per_splash;

  options(matwrap args =  matwrap(NULL)) :
    alg_type(CHROMATIC),  nsamples(10), nskip(10), 
    // tskip(5), 
    ncpus(2),
    // ntrees(ncpus),
    treewidth(3), 
    treeheight(std::numeric_limits<size_t>::max()), 
    treesize(std::numeric_limits<size_t>::max()), 
    priorities(false),  
    vanish(10),
    save_alchemy(false),
    ncpus_per_splash(1) {
    if(args.is_null()) return;
    safe_assert(args.is_struct(), 
                "Additional arguments must be in a struct");
    { // parse the sampler algorithm type
      matwrap arg(args.get_field("alg_type"));
      if(!arg.is_null()) {
        const size_t str_len(255);
        char sampler_type_c_str[str_len];
        arg.as_string(sampler_type_c_str, str_len);
        std::string sampler_type_str(sampler_type_c_str, arg.size());

        if(sampler_type_str == "CHROMATIC") alg_type = CHROMATIC;
        else if(sampler_type_str == "SPLASH") alg_type = SPLASH;
        else {
          char error_str[2*str_len];
          std::sprintf(error_str, "Invalid sampler type: (%s)\n", 
                       sampler_type_c_str);
          mexErrMsgTxt(error_str);
        }        
      }
    } // end of parse field name
    { // parse the number of samples
      matwrap arg(args.get_field("nsamples"));
      if(!arg.is_null()) {
        nsamples = size_t(arg.get_double_array()[0]);    
      }
    } // end of parse field name
    { // parse the number of skipped samples
      matwrap arg(args.get_field("nskip"));
      if(!arg.is_null()) {
        nskip = size_t(arg.get_double_array()[0]);    
      }
    } // end of parse field name
    // { // parse the number of skipped samples
    //   matwrap arg(args.get_field("tskip"));
    //   if(!arg.is_null()) {
    //     tskip = arg.get_double_array()[0];    
    //   }
    // } // end of parse field name                        
    { // parse the number of cpus
      matwrap arg(args.get_field("ncpus"));
      if(!arg.is_null()) {
        ncpus = size_t(arg.get_double_array()[0]);    
      }
    } // end of parse field name
    // { // parse the number of trees
    //   matwrap arg(args.get_field("ntrees"));
    //   if(!arg.is_null()) {
    //     ntrees = size_t(arg.get_double_array()[0]);    
    //   }
    // } // end of parse field name
    { // parse the treewidth
      matwrap arg(args.get_field("treewidth"));
      if(!arg.is_null()) {
        treewidth = size_t(arg.get_double_array()[0]);    
      }
    } // end of parse field name
    { // parse the treeheight
      matwrap arg(args.get_field("treeheight"));
      if(!arg.is_null()) {
        treeheight = size_t(arg.get_double_array()[0]);    
      }
      ASSERT_GT(treeheight, 0);
    } // end of parse field name
    { // parse the treesize
      matwrap arg(args.get_field("treesize"));
      if(!arg.is_null()) {
        treesize = size_t(arg.get_double_array()[0]);    
      }
      ASSERT_GT(treesize, 0);
    } // end of parse field name
    { // parse the priorities
      matwrap arg(args.get_field("priorities"));
      if(!arg.is_null()) {
        priorities = bool(arg.get_double_array()[0]);    
      }
    } // end of parse field name                         
    { // parse the vanish
      matwrap arg(args.get_field("vanish"));
      if(!arg.is_null()) {
        vanish = size_t(arg.get_double_array()[0]);    
      }
    } // end of parse field name
    { // parse the ncpus_per_splash
      matwrap arg(args.get_field("npcus_per_splash"));
      if(!arg.is_null()) {
        ncpus_per_splash = size_t(arg.get_double_array()[0]);    
      }
    } // end of parse field name                         
    { // parse the ncpus_per_splash
      matwrap arg(args.get_field("save_alchemy"));
      if(!arg.is_null()) {
        save_alchemy = bool(arg.get_double_array()[0]);    
      }
    } // end of parse field name                         

  } // end of constructor

  void print() {
    std::cout << "Generating "
              << nsamples << " samples " << std::endl;
    std::cout << "Skipping every " << nskip << " samples." << std::endl;
    switch(alg_type) {
    case CHROMATIC: mexPrintf("Using the chromatic sampler.\n"); break;
    case SPLASH: mexPrintf("Using the splash sampler.\n"); 
      std::cout << "ncpus:        "  << ncpus << std::endl
        //                << "ntrees:       "  << ntrees << std::endl
                << "treewidth:    "  << treewidth << std::endl
                << "treeheight:   " << treeheight << std::endl
                << "prioritizeds:   " << (priorities ? "enabled" : "disabled")
                << std::endl
                << "ncpus/splash: " << ncpus_per_splash << std::endl;

      break;
    default: mexErrMsgTxt("No algorithm selected!\n"); break;
    }
    flush_screen();
  }

}; // end of options struct









///////////////////////////////////////////////////////////////////////////
// Factor Graph generation code

/**
 * Build a graphlab table factor from the matlab table factor
 */
void add_factor(factorized_model& model,
                matwrap matlab_factor) {                
  if(matlab_factor.is_null()) {
    mexErrMsgTxt("Null factor argument to build factor!\n");    
  }
  if(!matlab_factor.is_struct()) {
    mexErrMsgTxt("Invalid factor type");
  }
  // Get the field numbers
  if(vars_field_id < 0) {
    vars_field_id = matlab_factor.get_field_number(vars_field_name);
  }
  if(logP_field_id < 0) {
    logP_field_id = matlab_factor.get_field_number(logP_field_name);
  }
  safe_assert(vars_field_id >= 0, "No field named vars in factor struct!");
  safe_assert(logP_field_id >= 0, "No field named logP in factor struct!");

  // Load the members from the factor
  matwrap variables = matlab_factor.get_field(vars_field_id);
  matwrap logP = matlab_factor.get_field(logP_field_id);

  // Check the fields 
  safe_assert(variables.is_uint32(), "Variables are not uint32t");
  safe_assert(logP.is_double(), "logP must be of type double");

  // Get information about each field
  const mwSize num_vars = variables.size();
  const mwSize* var_dims = logP.get_dimensions();
  const mwSize num_dims = logP.get_num_dimensions();
  safe_assert(num_vars <= MAX_DIM,
              "Too many (>32) variables in factor.");


  // Build up the domain
  domain_t dom;
  const uint32_t* varids = variables.get_data<uint32_t>();
  for(size_t i = 0, j = 0; i < num_vars && j < num_dims; ++i, ++j) {
    // Skip all the empty dimensions
    while(var_dims[j] <= 1 && j < num_dims) j++;
    safe_assert(j < num_dims, "Factor dimensions do not match variables");
    safe_assert(varids[i] > 0, "Variabile ids must start at 1");
    variable_t var(varids[i]-1, var_dims[j]);
    dom += var;
  }

  //  Create the factor
  factor_t& factor(model.add_factor(dom));
  const size_t num_logP(logP.size());
  const double* data(logP.get_double_array());
  safe_assert(num_logP == dom.size(), 
          "Insufficient factor data to match variables");
  for(size_t i = 0; i < num_logP; ++i) {
    factor.logP(i) = data[i];
  }
  factor.normalize();

}


void build_factorized_model(factorized_model& model, matwrap factors) {
  const size_t num_factors(factors.size());
  model.reserve(num_factors);
  // Load all the factors
  for(size_t i = 0; i < num_factors; ++i) {
    add_factor(model, factors.get_cell(i));   
  }
}







///////////////////////////////////////////////////////////////////////////
// Code for sample collection
struct result_collector {
  size_t sample_id;
  matwrap samples;
  matwrap beliefs;
  matwrap nsamples; 
  matwrap nchanges; 
 
  result_collector(matwrap samples = NULL, 
                   matwrap beliefs = NULL,
                   matwrap nsamples = NULL,
                   matwrap nchanges = NULL) :
    sample_id(0), samples(samples), beliefs(beliefs), 
    nsamples(nsamples), nchanges(nchanges)  { }
};

graphlab::glshared<result_collector> glshared_collector;

void collector_sync(mrf_gl::iscope& scope,
                    graphlab::any& accumulator) {
  result_collector collector = glshared_collector.get_val();
  const size_t& sample_id = collector.sample_id;
  const mrf_vertex_data& vdata(scope.const_vertex_data());
  if(!collector.samples.is_null() && 
     sample_id < collector.samples.cols()) {
    collector.samples.mat_index2d(vdata.variable.id(), 
                                  sample_id) = vdata.asg + 1;     
  }
  if(!collector.beliefs.is_null() && 
     sample_id < collector.beliefs.cols()) {
    matwrap blf(collector.beliefs.
                       get_cell_index2d(vdata.variable.id(), sample_id));
    double sum = 0;
    for(size_t i = 0; i < vdata.variable.size(); ++i) {
      sum += (blf.get_double_array()[i] = exp(vdata.belief.logP(i)));
    }
    for(size_t i = 0; i < vdata.variable.size(); ++i) {
      blf.get_double_array()[i] /= sum;
    }
  }
  if(!collector.nsamples.is_null() && 
     sample_id < collector.nsamples.cols()) {
    collector.nsamples.mat_index2d(vdata.variable.id(),
                                   sample_id) = vdata.nsamples;
  }
  if(!collector.nchanges.is_null() && 
     sample_id < collector.nchanges.cols()) {
    collector.nchanges.mat_index2d(vdata.variable.id(),
                                   sample_id) = vdata.nchanges;
  }
}

void collector_apply(graphlab::any& current_data, 
                     const graphlab::any& param) {
  result_collector& collector = current_data.as<result_collector>();
  collector.sample_id++;
}

void collector_merge(graphlab::any& merge_dest, 
                     const graphlab::any& merge_src) {
  // nop
}


size_t global_termination_nsamples = 0;
bool nsamples_terminator() {
  return glshared_collector.get_val().sample_id >= 
    global_termination_nsamples;
}









void run_chromatic_sampler(mrf_gl::core& core, 
                           const options& opts) {
  
  core.sched_options().add_option("update_function", 
                                  single_gibbs_update);
  // core.sched_options().add_option("max_iterations", 
  //                                 opts.nsamples * opts.nskip);
  global_termination_nsamples = opts.nsamples;
  core.engine().add_terminator(nsamples_terminator);
  std::cout << "Running." << std::endl;
  flush_screen();
  double runtime = core.start();    
  std::cout << "Runtime: " << runtime << std::endl;
  flush_screen();
}  // end of run_chromatic_sampler



void run_jt_splash_sampler(mrf_gl::core& core, 
                           const options& opts) {

  std::cout << "Starting Splash Sampler." << std::endl;
  flush_screen();
  splash_settings settings;
  settings.ntrees           = opts.ncpus;
  settings.max_tree_width   = opts.treewidth;
  settings.max_tree_height  = opts.treeheight;
  settings.max_tree_size    = opts.treesize;
  settings.priorities       = opts.priorities;
  settings.subthreads       = opts.ncpus_per_splash;

  jt_splash_sampler jtsplash_sampler(core.graph(),
                                     settings);
  
  const size_t samples_per_iteration = 
    core.graph().num_vertices() * opts.nskip;
  for(size_t i = 0; i < opts.nsamples; ++i) {
    //    std::cout << "Running an iteration" << std::endl;
    flush_screen();
    // run the splash sampler
    jtsplash_sampler.sample_updates(samples_per_iteration);
    // jtsplash_sampler.sample_seconds(opts.tskip);
    // compute the sync
    core.sync_now(glshared_collector);
    // std::cout << "Ntrees: " << jtsplash_sampler.total_trees() 
    //           << std::endl
    //           << "Nsamples: " << jtsplash_sampler.total_samples()
    //           << std::endl
    //           << "NCollisions: " << jtsplash_sampler.total_collisions()
    //           << std::endl;
    
  }
  
} // end of run splash sampler





/**
 * See parallel_gibbs.m for identical arguments
 * 
 *   [samples, nupdates, nchanges, marginals] = ...
 *     gibbs_sampler(factors, options); 
 * 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  safe_assert(nrhs > 0, 
              "Invalid number of arguments. See documentation for details");  

  { // Seed the random number generator using rand from within matlab
    mxArray* plhs[1];
    mxArray* prhs[4];
    prhs[0] = mxCreateNumericMatrix(1,1, mxUINT32_CLASS, mxREAL);
    prhs[1] = mxCreateNumericMatrix(1,1, mxUINT32_CLASS, mxREAL);
    prhs[2] = mxCreateNumericMatrix(1,1, mxUINT32_CLASS, mxREAL);
    matwrap(prhs[0]).get_data<uint32_t>()[0] = 
      std::numeric_limits<uint32_t>::max();
    matwrap(prhs[1]).get_data<uint32_t>()[0] = 1;
    matwrap(prhs[2]).get_data<uint32_t>()[0] = 2;
    prhs[3] = mxCreateString("uint32");
    mexCallMATLAB(1, plhs, 4, prhs, "randi");
    const size_t seed_value = 
      matwrap(plhs[0]).get_data<size_t>()[0];
    //    std::cout << "Seed value: " << seed_value << std::endl;
    mxDestroyArray(plhs[0]);
    mxDestroyArray(prhs[0]); 
    mxDestroyArray(prhs[1]); 
    mxDestroyArray(prhs[2]);
    mxDestroyArray(prhs[3]);
    graphlab::random::seed(seed_value);
  }


  // Get first argument the cell array of factors
  const matwrap matlab_factors(const_cast<mxArray*>(prhs[0]));
  if(matlab_factors.is_null()) { 
    mexErrMsgTxt("No factors provided!\n"); 
  }
  if(!matlab_factors.is_cell()) { 
    mexErrMsgTxt("Factors must be in cell array form!\n"); 
  }
  
  options opts;
  if(nrhs > 1) opts = options(const_cast<mxArray*>(prhs[1]));
  opts.print();
 


  // Load the factorized model
  factorized_model model;
  build_factorized_model(model, matlab_factors);
  std::cout << "Finished Loading Factors" << std::endl;
  flush_screen();
  if(opts.save_alchemy) {
    std::cout << "Saving Alchemy file \"problem.alchemy\"" << std::endl;
    model.save_alchemy("problem.alchemy");
  }

  // mexPrintf("Finished Saving Alchemy File\n");

  // Set the global factors
  SHARED_FACTORS_PTR = &(model.factors());


  // Construct the markov random field
  mrf_gl::core mrf_core;
  mrf_from_factorized_model(model, mrf_core.graph());    

  if(opts.alg_type == options::CHROMATIC) {
    const size_t ncolors(mrf_core.graph().compute_coloring());  
    std::cout << "Finished coloring graph with " << ncolors 
              << " colors." << std::endl;
    flush_screen();
    mrf_core.set_ncpus(opts.ncpus);
    mrf_core.set_scheduler_type("chromatic");
    mrf_core.set_scope_type("null");
    mrf_core.engine().set_sched_yield(false);
  }


  // allocate the return matrix for samples
  matwrap matlab_samples;
  if(nlhs > 0) {
    matlab_samples = 
      matwrap::create_matrix(model.variables().size(),
                                    opts.nsamples);
    plhs[0] = matlab_samples.array;
    safe_assert(!matlab_samples.is_null(), 
                "Error initializing return samples");
  }
  if( !matlab_samples.is_null() ) {
    double* entries(matlab_samples.get_double_array());
     const size_t num_entries(matlab_samples.size());
    for(size_t i = 0; i < num_entries; ++i) {
      entries[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }


  // allocate the return matrix for nsamples
  matwrap matlab_nsamples;
  if(nlhs > 1) {
    matlab_nsamples = 
      matwrap::create_matrix(model.variables().size(),
                                    opts.nsamples);
    plhs[1] = matlab_nsamples.array;
    safe_assert(!matlab_nsamples.is_null(), 
                "Error initializing return nsamples");
  }
  if( !matlab_nsamples.is_null() ) {
    double* entries(matlab_nsamples.get_double_array());
     const size_t num_entries(matlab_nsamples.size());
    for(size_t i = 0; i < num_entries; ++i) {
      entries[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  // allocate the return matrix for nsamples
  matwrap matlab_nchanges;
  if(nlhs > 2) {
    matlab_nchanges = 
      matwrap::create_matrix(model.variables().size(),
                                    opts.nsamples);
    plhs[2] = matlab_nchanges.array;
    safe_assert(!matlab_nchanges.is_null(), 
                "Error initializing return nchanges");
  }
  if( !matlab_nchanges.is_null() ) {
    double* entries(matlab_nchanges.get_double_array());
     const size_t num_entries(matlab_nchanges.size());
    for(size_t i = 0; i < num_entries; ++i) {
      entries[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }


  // allocate the return cell array for beliefs
  matwrap matlab_beliefs;
  if(nlhs > 3) {
    matlab_beliefs = 
      matwrap::create_cell(model.variables().size(),
                                  opts.nsamples);
    plhs[3] = matlab_beliefs.array;
    safe_assert(!matlab_beliefs.is_null(), 
                "Error initializing return beliefs");
  }
  if( !matlab_beliefs.is_null() ) {
    //populate each of the entries
    for(size_t i = 0; i < mrf_core.graph().num_vertices(); ++i) {
      const mrf_vertex_data& vdata = 
        mrf_core.graph().vertex_data(i);      
      for(size_t j = 0; j < opts.nsamples; ++j) {
        matwrap blf = 
          matwrap::create_matrix(vdata.variable.size(),
                                        1);
        safe_assert(!blf.is_null(), "Unable to allocate beliefs");
        matlab_beliefs.set_cell_index2d(i, j, blf);
      }
    }
  }




  //! initialize results collection
  glshared_collector.set(result_collector(matlab_samples, 
                                          matlab_beliefs, 
                                          matlab_nsamples, 
                                          matlab_nchanges));


  // if(!matlab_samples.is_null()) {
  //   std::cout << "Enabling sample collection." << std::endl;
  //   flush_screen();
    //    last_tic = graphlab::lowres_time_seconds();
  const size_t sync_interval = 
    opts.nskip * mrf_core.graph().num_vertices();
  mrf_core.set_sync(glshared_collector,
                    collector_sync,
                    collector_apply,
                    graphlab::any(size_t(0)),
                    sync_interval,
                    collector_merge);
    //}



  if(opts.alg_type == options::CHROMATIC) {
    run_chromatic_sampler(mrf_core, opts);
  } else if( opts.alg_type == options::SPLASH) {
    run_jt_splash_sampler(mrf_core, opts);
  } else {
    mexErrMsgTxt("Invalid algorithm type.\n");
  }


   

} // end of main


