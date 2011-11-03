#ifndef FACTOR_GIBBS_HPP
#define FACTOR_GIBBS_HPP

#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>


// graphlab utitilities


#include <graphlab.hpp>



#include "factor_graph.hpp"



#include <graphlab/macros_def.hpp>


static const double DEFAULT_RESIDUAL(1.0);
static const size_t MAX_COLOR_ID(0);
static const size_t FACTOR_OFFSET(1);
static const size_t BUFFER_SIZE(sizeof(double)*5000);


// These should be set by the main
extern size_t total_colors;
extern size_t nsamples;

// Structs
// ===========================================================================>

struct edge_data {
  /// No edge data
};

/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data {
  graphlab::discrete_variable variable;
  uint16_t asg;
  uint16_t color;
  uint32_t samples;
  factor_graph::factor_type belief;
  std::vector<uint32_t> counts;
  std::vector<uint32_t> factors;
  vertex_data() :
    asg(0), color(0), samples(0) { }
}; // End of vertex data


typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;


// Update Functions
// ===========================================================================>

/**
 * The gibbs update draws a random new sample for this vertex
 * conditioned on the values the neighboring vertices.  Depending on
 * the scheduler, the callback mechanism may or may not be used.
 * Therefore we templatize this function over wheter a callback is
 * used.
 */
template<bool use_callback>
void gibbs_update(gl_types::iscope& scope, 
                  gl_types::icallback& scheduler,
                  gl_types::ishared_data* sdm) {
  // Define types
  typedef factor_graph::factor_type     factor_type;
  typedef factor_type::domain_type      domain_type;
  typedef factor_type::assignment_type  assignment_type;
  typedef factor_graph::variable_type   variable_type;
  typedef std::map<graphlab::vertex_id_t, uint16_t>
    asg_map_type;
  typedef asg_map_type::const_iterator asg_map_iter_type;

  // Get the shared data
  assert(sdm != NULL);
   
  // Get the vertex data  
  vertex_data& vdata = scope.vertex_data();

  // If we have enough samples return
  if(vdata.samples > nsamples) return;
  
  
  // Get the out edges
  graphlab::edge_list out_edges = scope.out_edge_ids();

  // Get all the neighbor vertices assignments
  // ------------------------------------------------------->
  // Build an assignment map
  asg_map_type asg_map;
  for(size_t i = 0; i < out_edges.size(); ++i) {
    graphlab::vertex_id_t target = scope.target(out_edges[i]);
    const vertex_data& neighbor_vdata =
      scope.neighbor_vertex_data(target);
    assert(neighbor_vdata.color != vdata.color);
    asg_map[target] = neighbor_vdata.asg;    
  }
  
  // Condition on all the neighboring vertices
  // ------------------------------------------------------->
  factor_type marginal(vdata.variable);
  factor_type conditional(vdata.variable);
  marginal.uniform();
    
  // Looping over each factor
  for(size_t i = 0; i < vdata.factors.size(); ++i) {
    // Get the factor
    uint32_t factor_id = vdata.factors[i];
    // The factors are all stored as constants in the shared data
    // manager
    const factor_type& factor =
      sdm->get_constant(FACTOR_OFFSET + factor_id).as<factor_type>();

    // Construct the restricted assignment for the factor
    assignment_type asg;
    for(size_t j = 0; j < factor.num_vars(); ++j) {
      const variable_type& other_variable = factor.args().var(j);
      // if the other variable is not this variable
      if(vdata.variable != other_variable) {
        // Get the assignment in the large neighborhood assignment map
        asg_map_iter_type iter = asg_map.find(other_variable.id());
        assert(iter != asg_map.end());
        // union the assignment 
        asg &= assignment_type(other_variable, iter->second);
      }
    }
    conditional.condition(factor, asg);    
    // and multiply into the marginal
    marginal *= conditional;
  }
  // noramlzie the marginal 
  marginal.normalize();

  // Draw a sample
  // ------------------------------------------------------->
  assignment_type sample = marginal.sample();
  // Slight hack: the sample should be one dimensional hence:
  // asg_at(0)
  vdata.asg = sample.asg_at(0);
  // Update the counts
  vdata.counts[vdata.asg]++;
  // Record the belief (this is Rao-Blackwellized estimator)
  vdata.belief += marginal;
  vdata.samples++;
  
  // Reschedule self if necessary
  // ------------------------------------------------------->
  // If the set scheduler is used then a callback is not needed since
  // the scheduler uses a predetermined fixed ordering.  However other
  // schedulers, like fifo rely on the callback to know what the
  // schedule next.
  if(use_callback) {
    if(vdata.samples < nsamples) {
      gl_types::update_function update_function =
        gibbs_update<use_callback>;
      gl_types::update_task task(scope.vertex(),
                                 update_function);      
      scheduler.add_task(task, DEFAULT_RESIDUAL);
    } 
  }
} // end of gibbs_update




// The color update is used to color the clique graph MRF
void color_update(gl_types::iscope& scope, 
                  gl_types::icallback& scheduler,
                  gl_types::ishared_data* sdm) {
  // Get the vertex data
  vertex_data& vdata = scope.vertex_data();
  
  // Get the out edges
  graphlab::edge_list out_edges =  scope.out_edge_ids();
  
  std::set<uint16_t> neighbor_colors;
  // Collect the neighbors colors
  for(size_t i = 0; i < out_edges.size(); ++i) {
    graphlab::vertex_id_t target = scope.target(out_edges[i]);
    const vertex_data& neighbor_vdata =
      scope.neighbor_vertex_data(target);
    neighbor_colors.insert(neighbor_vdata.color);
  }
  
  // Greedily color the vertex
  vdata.color = 0;
  foreach(uint16_t color, neighbor_colors) {
    if(vdata.color < color) break;
    else {
      assert(vdata.color == color);
      vdata.color++;
    }
    // Ensure no wrap around (that would be tragic)
    assert(vdata.color > 0);
  }
}


//! Get the color from the vertex data (used in the sync operation)
uint16_t get_color(const vertex_data& vdata) { return vdata.color; }

uint16_t parallel_graph_color(gl_types::graph& graph,
                              size_t ncpus) {
  gl_types::iengine* engine =
    graphlab::engine_factory::new_engine("threaded",
                                         "multiqueue",
                                         "locked",
                                         graph,
                                         ncpus);

  assert(engine != NULL);
  gl_types::thread_shared_data sdm;

  // Setup a sync to compute the maximum vertex color
  typedef uint16_t accum_type;
  accum_type zero = 0;
  // Use a prefab sync_function that takes the max (uses the function
  // get_color) to read the user vertex data
  gl_types::ishared_data::sync_function_type sync_fun =
    gl_types::sync_ops::max<accum_type, get_color>;
  // Use the prefab identity apply function which simply applies the
  // new value without any additional transformations
  gl_types::ishared_data::apply_function_type apply_fun =
    gl_types::apply_ops::identity<accum_type>;
  // Store the zero in the field initially
  sdm.create_atomic(MAX_COLOR_ID, graphlab::any(zero));
  // Setup the sync 
  sdm.set_sync(MAX_COLOR_ID,
               sync_fun,
               apply_fun,
               graphlab::any(zero));
  // Register the shared data
  engine->set_shared_data_manager(&sdm);
  double initial_priority = 1.0;
  engine->add_task_to_all(color_update, initial_priority);
  engine->start();
  // Run an extra sync at the end to compute the largest color
  sdm.sync(MAX_COLOR_ID);
  uint16_t num_colors = sdm.get(MAX_COLOR_ID).as<uint16_t>() + 1;
  delete engine;
  return num_colors;
} // end of parallel graph color



/**
 * This constructs the Markov Random Field from the factor graph.  The
 * markov random field is not a pairwise MRF but the clique graph
 * corresponding to the factor graph.
 */
void make_gibbs_graph(const factor_graph& fgraph,
                      gl_types::graph& graph) {
  assert(!fgraph.variables().empty());
  assert(!fgraph.factors().empty());
  assert(fgraph.variables().rbegin()->id() == fgraph.variables().size()-1);
  
  typedef factor_graph::factor_type factor_type;
  typedef factor_graph::domain_type domain_type;

  // Create a reverse map which maps from varialbes to the factors
  // that they are member of:
  typedef std::vector< std::set<uint32_t> > reverse_map_type;
  reverse_map_type reverse_map(fgraph.variables().size());

  
  // Find all the factors for each variable
  for(uint32_t i = 0; i < fgraph.factors().size(); ++i) {
    const domain_type& domain = fgraph.factors()[i].args();
    for(uint32_t j = 0; j < domain.num_vars(); ++j) {
      reverse_map[domain.var(j).id()].insert(i);
    }
  }

  vertex_data vdata;
  // Add all the variables
  foreach(factor_graph::variable_type variable, fgraph.variables()) {
    // Set the arity
    vdata.variable = variable;
    vdata.belief.set_args(variable);
    vdata.counts.resize(variable.size());
    // Make a uniformly "zero" factor initially.
    vdata.belief.uniform(-std::numeric_limits<double>::max());
    // Set all the factors
    size_t num_factors = reverse_map[variable.id()].size();
    vdata.factors.resize(num_factors);
    size_t index = 0;    
    foreach(uint32_t factor_id, reverse_map[variable.id()]) {
      vdata.factors[index++] = factor_id;
    }
    assert(num_factors == index);
    graphlab::vertex_id_t vid = graph.add_vertex(vdata);
    assert(vid == variable.id());
  }  
  assert(graph.num_vertices() == fgraph.variables().size());  

  std::set< std::pair<uint32_t, uint32_t> > edge_set;
  // Add all the edges
  foreach(const factor_graph::factor_type& factor, fgraph.factors()) {
    const domain_type& domain = factor.args();
    // Make a clique over the domain
    for(size_t i = 0; i < domain.num_vars(); ++i) {      
      graphlab::vertex_id_t vi = domain.var(i).id();
      for(size_t j = i+1; j < domain.num_vars(); ++j) {
        graphlab::vertex_id_t vj = domain.var(j).id();
        assert(vi != vj);
        if(edge_set.count(std::make_pair(vi,vj)) == 0) {
          graph.add_edge(vi, vj);
          graph.add_edge(vj, vi);
          edge_set.insert(std::make_pair(vi,vj));
        }
      }
    }
  }
  // Finalize the graph
  graph.finalize();
} // End of make Gibbs graph

/**
 * We save all the factors to the shared data instead of in the graph
 * since they will remain constant throughout the exectution.
 */
void fill_shared_data(const factor_graph& fgraph,
                      gl_types::thread_shared_data& sdm) {
  typedef factor_graph::factor_type factor_type;
  // Add all the factors to the shared data manager
  for(size_t i = 0; i < fgraph.factors().size(); ++i) {
    const factor_type& factor = fgraph.factors()[i];
    sdm.set_constant(i + FACTOR_OFFSET, graphlab::any(factor));                     
  }  
} // end of fill shared data


/**
 * Save all the beliefs to a file
 */
void save_beliefs(const std::string& filename,
                  factor_graph& fgraph,
                  const gl_types::graph& graph) {
  // Open the file to store the final belief vectors
  std::ofstream fout(filename.c_str());
  assert(fout.good());
  // Save all the beliefs
  for(size_t i = 0; i < graph.num_vertices(); ++i) {
    const vertex_data& vdata = graph.vertex_data(i);
    fout << fgraph.var_name(i) << " // ";    
    size_t sum = 0;
    // normalize the probabilities
    for(size_t j = 0; j < vdata.counts.size(); ++j) sum += vdata.counts[j];
    assert(sum == vdata.samples);
    // Save the normalized belief estimates to the file
    for(size_t j = 0; j < vdata.counts.size(); ++j) {
      fout << (double(vdata.counts[j])/sum);
      if((j + 1) < vdata.counts.size()) fout << ", ";
    }
    fout << "\n";
  }
  fout.close();   
} // end of save beliefs


// Selector for vertex color
bool select_color(uint32_t color,
                  gl_types::vertex_id_t vertex,
                  const vertex_data& vdata) {
  return vdata.color == color;
}


/* TODO: set scheduler is GONE!

// Set scheduler colored schedule
void basic_color_schedule(gl_types::set_scheduler &sched) { 
  // All sets must be created before scheduling calls
  std::vector<gl_types::ivertex_set*> colorsets;
  // Construct color sets
  for (size_t i = 0; i < total_colors; ++i){
    colorsets.push_back(&sched.attach(gl_types::rvset(gl_types::selector_function(boost::bind(select_color, i, _1, _2)), true), sched.root_set()));
  }

  const bool use_callback = false;
  
  gl_types::update_function update_function = gibbs_update<use_callback>;

  
  // Actually construct color sets
  sched.init();
  for (size_t i = 0; i < total_colors; ++i){
    std::cout <<"Set " << i <<" has " <<colorsets[i]->size() << "vertices\n";
  }
  
  // Actually execut the schedule generating all the samples
  for (size_t iter = 0; iter < nsamples; ++iter) {
    for (size_t i = 0; i < total_colors; ++i){
      sched.execute(*colorsets[i], update_function);
      sched.wait();
    }
    
  }
}

// Set scheduler colored schedule
void planned_color_schedule(gl_types::set_scheduler &sched) {

  
  // All sets must be created before scheduling calls
  std::vector<gl_types::ivertex_set*> colorsets;
  // Construct color sets
  for (size_t i = 0; i < total_colors; ++i){
    colorsets.push_back(&sched.attach(gl_types::rvset(gl_types::selector_function(boost::bind(select_color, i, _1, _2)), true), sched.root_set()));

  }
  
  gl_types::execution_plan eplan;

  const bool use_callback = false;
  gl_types::update_function update_function =
    gibbs_update<use_callback>;

  
  // Actually construct color sets
  sched.init();
  
  for (size_t i = 0; i < total_colors; ++i) {
    std::cout <<"Set " << i <<" has " <<colorsets[i]->size() << "vertices\n";
    eplan.execute(*colorsets[i], update_function);
  }
  
  graphlab::timer ti;
  ti.start();
  eplan.generate_plan(sched.get_graph(), sched.num_cpus());
  std::cout << "Execution Plan took " << ti.current_time() << " to compile\n";


  // Actually execut the schedule generating all the samples
  for (size_t iter = 0; iter < nsamples; ++iter) {
    sched.execute_plan(eplan);
  }
}

*/


// //! color the graph sequentially
// uint16_t sequential_graph_color(gl_types::graph& graph) {
//   uint16_t max_color = 0;

//   // Compute a random order
//   std::vector<uint32_t> order(graph.num_vertices());
//   for(size_t i = 0; i < order.size(); ++i) order[i] = i;
//   std::random_shuffle(order.begin(), order.end());
  
//   //  for(size_t i = 0; i < graph.num_vertices(); ++i) {
//   foreach(uint32_t i, order){
//     vertex_data& vdata = graph.vertex_data(i);
//     std::set<uint16_t> neighbor_colors;
//     // Collect the neighbors colors
//     const std::vector<graphlab::edge_id_t>& out_edges =
//       graph.out_edge_ids(i);
//     for(size_t j = 0; j < out_edges.size(); ++j) {
//       graphlab::vertex_id_t target = graph.target(out_edges[j]);
//       const vertex_data& other = graph.vertex_data(target);
//       neighbor_colors.insert(other.color);
//     }
//     // Greedily color the vertex
//     vdata.color = 0;
//     foreach(uint16_t color, neighbor_colors) {
//       if(vdata.color < color) break;
//       else {
//         assert(vdata.color == color);
//         ++(vdata.color);
//       }
//       // Ensure no wrap around (that would be tragic)
//       assert(vdata.color > 0);
//     }
//     // Determine the max color
//     max_color = std::max(vdata.color, max_color);
//   }
//   return max_color+1;
// }




#include <graphlab/macros_undef.hpp>
#endif
