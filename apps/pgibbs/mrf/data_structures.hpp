#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

/**
 *
 * Parallel blocked gibbs using graphlab
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>

// Including Standard Libraries

#include <iostream>
#include <iomanip>



#include <graphlab.hpp>


// #include "image.hpp"


typedef graphlab::vertex_id_t    vertex_id_t;
typedef graphlab::edge_id_t      edge_id_t;
typedef graphlab::unary_factor   unary_factor;
typedef graphlab::binary_factor  binary_factor;

// Represents a null VID in the tree
const vertex_id_t NULL_VID = -1;

// STRUCTS (Edge and Vertex data) =============================================>
/**
 * The data associated with each variable in the pairwise markov
 * random field
 */
struct vertex_data {

  enum vertex_state {
    AVAILABLE,         // Vertex is completely available [Default]

    CANDIDATE,         // The vertex is currently a candidate
    
    BOUNDARY,          // The vertex is on the boundary of a tree

    TREE_NODE,         // The vertex is in a tree
    
    CALIBRATED         // The vertex has been calibrated and has 
                       // computed the message to the parent 
  };


  uint16_t asg;
  size_t updates;

  //! The node potential associated with the vertex
  unary_factor potential;
  
  //! store the Roa-Blackwell conditional belief estimate
  unary_factor belief;
  

  unary_factor bp_marginal;

  //! the parent in the tree
  vertex_id_t parent;
  vertex_state state;
  size_t height;
  graphlab::atomic<size_t> child_candidates;
  

  
  vertex_data() :
    asg(0), updates(0),
    parent(NULL_VID),
    state(AVAILABLE),
    height(0),
    child_candidates(0) { }

  void print() const {
    std::cout << "Parent: " << parent << ", "
              << "State: " << state << ", "
              << "Height: " << height << std::endl;
  }

  void save(graphlab::oarchive &arc) const {
    arc << asg;
    arc << updates;
    arc << potential;
    arc << belief;
    arc << parent;

    arc.o->write(reinterpret_cast<const char*>(&state), sizeof(vertex_state));
    arc << height;
    size_t value = child_candidates.value;
    arc << value;
  }
  
  void load(graphlab::iarchive &arc) {
    arc >> asg;
    arc >> updates;
    arc >> potential;
    arc >> belief;
    arc >> parent;
    
    arc.i->read(reinterpret_cast<char*>(&state), sizeof(vertex_state));
    arc >> height;
    size_t value;
    arc >> value;
    child_candidates.value = value;
    
  }
}; // End of vertex data



/**
 * The data associated with each directed edge in the pairwise markov
 * random field
 */
struct edge_data { 
  double weight;
  uint32_t factor_id;
  unary_factor message;
  bool   exploring;
  edge_data() : weight(0), factor_id(-1), exploring(false) { }

  void save(graphlab::oarchive &arc) const {
    arc << weight;
    arc << factor_id;
    arc << message;
    arc << exploring;
  }
  
  void load(graphlab::iarchive &arc) {
    arc >> weight;
    arc >> factor_id;
    arc >> message;
    arc >> exploring;
  }
}; 

// define the graph type:
typedef graphlab::graph< vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl;





std::string make_filename(const std::string& base,
                          const std::string& suffix,
                          size_t number) {
  std::stringstream strm;
  strm << base
       << std::setw(10) << std::setfill('0')
       << number
       << suffix;
  std::cout << strm.str() << std::endl;
  return strm.str();
}





/** Save the beliefs stored in the graph */
void save_beliefs(const graph_type& graph,
                  const std::string& filename) {
  std::ofstream fout(filename.c_str());
  fout.precision(8);
  graphlab::unary_factor marginal;
  for(size_t v = 0; v < graph.num_vertices(); ++v) {
    const vertex_data& vdata = graph.vertex_data(v);
    marginal = vdata.belief;
    marginal.normalize();
    fout << vdata.updates << '\t';
    for(size_t asg = 0; asg < marginal.arity(); ++asg) {
      fout << std::exp( marginal.logP(asg) );
      if(asg + 1 < marginal.arity() ) fout << '\t';      
    }
    fout << '\n';
  } 
  fout.close();
} // End of save beliefs


void save_asg(const graph_type& graph,
              const std::string& filename) {
  std::ofstream fout(filename.c_str());
  graphlab::unary_factor marginal;
  for(size_t v = 0; v < graph.num_vertices(); ++v) 
    fout << graph.vertex_data(v).asg << '\n';
  fout.close();
} // End of save beliefs


void save_color(const graph_type& graph,
                const std::string& filename) {
  std::ofstream fout(filename.c_str());
  graphlab::unary_factor marginal;
  for(size_t v = 0; v < graph.num_vertices(); ++v) 
    fout << graph.color(v) << '\n';
  fout.close();
} // End of save beliefs


//! Compute the unormalized likelihood of the current assignment
double unnormalized_likelihood(const graph_type& graph,
                               const gl::ishared_data& shared_data,
                               const size_t EDGE_FACTOR_OFFSET) {
  double sum = 0;

  // Compute all the vertex likelihoods
  for(vertex_id_t i = 0; i < graph.num_vertices(); ++i) {
    const vertex_data& vdata = graph.vertex_data(i);
    sum += vdata.potential.logP(vdata.asg);
  }

  // Compute all edge likelihoods
  for(edge_id_t i = 0; i < graph.num_edges(); ++i) {
    vertex_id_t source = graph.source(i);
    vertex_id_t target = graph.target(i);
    if(source < target) {
      const vertex_data& source_vdata = graph.vertex_data(source);
      const vertex_data& target_vdata = graph.vertex_data(target);
      const edge_data& edata = graph.edge_data(i);
      const binary_factor& edge_factor =
        shared_data.get_constant(EDGE_FACTOR_OFFSET + 
                                  edata.factor_id).as<binary_factor>();
      sum += edge_factor.logP(source, source_vdata.asg,
                              target, target_vdata.asg);
    }   
  }
  return sum;
}






#endif



