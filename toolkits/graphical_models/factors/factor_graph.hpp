/**  
 *  Software submitted by 
 *  Systems & Technology Research / Vision Systems Inc., 2013
 *
 *  Approved for public release; distribution is unlimited. [DISTAR Case #21428]
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



#ifndef VSI_FACTOR_GRAPH_HPP
#define VSI_FACTOR_GRAPH_HPP

#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>

//#include <graphlab/engine/engine_includes.hpp>
#include <graphlab/graph/graph_includes.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/vertex_program/vertex_program_includes.hpp>


#include "table_factor.hpp"
#include "bp_graph_data.h"



//#include <graphlab/macros_def.hpp>
namespace belief_prop {


/** 
 * Defines a factor_graph; i.e., a bipartite graph whose verticies can 
 * be divide into two disjoint sets: a set of variables V and a set of 
 * factors F. Undirected edges connect factors to variables. An edge 
 * exists between f and v if v is a member of the factor's domain.
 * 
 * A variable specifies a unary discrete probabiltiy mass function over 
 * a set of labels. A variable's PMF is assumed to be defined by a 
 * dense_table in which each label has a corresponding probability. 
 * 
 * A factor specifies a discrete joint probability mass function over 
 * a set of variables. The n-D PMF is defined by either a dense_table  
 * or a sparse_table. The values in the table are stored such that the 
 * first variable added to the domain iterates the fastest. 
 *
 * The typical usasge of this interface would consist of 1) adding a set
 * of variables by using add_variable(), 2) defining the prior distribution
 * over each variable by using one of the prior related methods (such as 
 * set_prior_for_variable() ), 3) adding a set of factors by using 
 * add_factor(), 4) constructing the distributed graph using make_bp_graph(), 
 * 5) running belief propagation on the distributed graph to propagate the 
 * evidence across the graph (outside the scope of this interface), and 
 * 6) loading the results using pull_beliefs_for_variables(). 
 * 
 * The resulting belief for a variable can be queried using 
 * belief_for_variable().logP() (once the evidence has been propagated 
 * across the distributed graph and the results pulled back into the  
 * factor_graph using pull_beliefs_for_variables() ). 
 *
 * \author Scott Richardson     10/2012
 *
 */
// NOTE a variable has a domain that spans only itself. a factor has
// a domain that spans its neighboring variables
template<size_t MAX_DIM>
class factor_graph {
  typedef typename graph_type<MAX_DIM>::type graph_type_t;
  typedef vertex_data<MAX_DIM>               vertex_data_t;
  typedef edge_data<MAX_DIM>                 edge_data_t;
  // dense_table is used to define the distributions over variables.  
  // *this shouldnt need to know specifically about sparse_table
  typedef graphlab::dense_table<MAX_DIM>     dense_table_t;
  typedef graphlab::table_base<MAX_DIM>      table_base_t;
  typedef graphlab::table_factor<MAX_DIM>    factor_type;
  typedef graphlab::discrete_variable        variable_type;

// variable related methods
public: 
  factor_graph() : _unique_var_id(0) { }

public: 
  /** Add a new discrete variable to the factor graph */
  // TODO rename to create_variable()
  variable_type add_variable(const size_t n_labels, 
      const std::string& default_var_name="") 
  {
    size_t id = _unique_var_id++;
    // assert that the id can be used to index into _factors
    DCHECK_EQ(id, num_factors());

    // Only store the variable in the local variable map if the variable (key) 
    // does not exist in the graph (map) already
//    var_map_const_iter_type it = _var_map.find(var_name);
//    if(it != _var_map.end()) {
//      std::cout << "WARNING: variable already exists for that name" << std::endl;
//      return it->second;
//    }

    // Create a new variable
    variable_type variable(id, n_labels);

    std::string var_name(default_var_name);
    if(var_name.empty()) {
      std::stringstream ss; ss << id;
      var_name = ss.str();
    } 
    // Save the factor to the factor graph
    vertex_data_t vert = add_vertex(variable, var_name);
    logstream(LOG_INFO) << "var_id=" << id 
                        << " description='" << vert.name << "'" << std::endl;

    return variable; 
  }

  /** 
   * Direct access to the variable's belief distribution.  
   * useful for initialization 
   */
  dense_table_t& belief_for_variable(const variable_type& var) {
    DCHECK_LT(var.id(), num_factors());

    // NOTE the variable's prior distribution is always dense
    dense_table_t* belief = 
        dynamic_cast<dense_table_t*>(factors()[var.id()].belief.table());
    DCHECK_NE(belief, NULL); 

    return *belief;
  }
  dense_table_t& belief_for_var(const variable_type& var) {
    return belief_for_variable(var);
  }
  /** 
   * Direct access to the variable's belief distribution.  
   */
  const dense_table_t& belief_for_variable(const variable_type& var) const {
    DCHECK_LT(var.id(), num_factors());

    // NOTE the variable's prior distribution is always dense
    dense_table_t const *const belief = 
        dynamic_cast<dense_table_t const *const>(
        factors()[var.id()].belief.table());
    DCHECK_NE(belief, NULL); 

    return *belief;
  }
  const dense_table_t& belief_for_var(const variable_type& var) const {
    return belief_for_variable(var);
  }

  /** 
   * Direct access to the variable's prior distribution.  
   * useful for initialization 
   */
  dense_table_t& prior_for_variable(const variable_type& var) {
    DCHECK_LT(var.id(), num_factors());

    // NOTE the variable's prior distribution is always dense
    dense_table_t* potential = 
        dynamic_cast<dense_table_t*>(factors()[var.id()].potential.table());
    DCHECK_NE(potential, NULL); 

    return *potential;
  }
  /** 
   * Direct access to the variable's prior distribution.  
   */
  const dense_table_t& prior_for_variable(const variable_type& var) const {
    DCHECK_LT(var.id(), num_factors());

    // NOTE the variable's prior distribution is always dense
    dense_table_t const *const potential = 
        dynamic_cast<dense_table_t const *const>(factors()[var.id()].potential.table());
    DCHECK_NE(potential, NULL); 

    return *potential;
  }

  // TODO rename to stage_prior_for_variable, etc. 
  void set_prior_for_variable(const variable_type& var, 
      const std::vector<double>& data) 
  {
    DCHECK_LT(var.id(), num_factors());
    factor_type& potential = factors()[var.id()].potential;
    // NOTE the variable's prior distribution is always dense
    dense_table_t* table = dynamic_cast<dense_table_t* >(potential.table());
    DCHECK_NE(table, NULL);

    *table = dense_table_t(table->domain(), data);
  }

  void set_belief_for_variable(const variable_type& var, 
      const std::vector<double>& data) 
  {
    DCHECK_LT(var.id(), num_factors());
    factor_type& belief = factors()[var.id()].belief;
    // NOTE the variable's prior distribution is always dense
    dense_table_t* table = dynamic_cast<dense_table_t* >(belief.table());
    DCHECK_NE(table, NULL);

    *table = dense_table_t(table->domain(), data);
  }

  void set_prior_for_variable(const variable_type& var, 
      const dense_table_t& table) 
  {
    DCHECK_LT(var.id(), num_factors());

    factors()[var.id()].potential = factor_type(factor_type::DENSE_TABLE, table);
  }

  void set_belief_for_variable(const variable_type& var, 
      const dense_table_t& table) 
  {
    DCHECK_LT(var.id(), num_factors());

    factors()[var.id()].belief = factor_type(factor_type::DENSE_TABLE, table);
  }

  variable_type get_variable(const size_t id) {
    DCHECK_LT(id, num_factors());
    DCHECK_EQ(factors()[id].potential.table()->ndims(), 1);
    return factors()[id].potential.table()->var(0);
  }

private:
  // REVIEW prevent a variable from being double added
  const vertex_data_t& add_vertex(const variable_type &variable, 
      const std::string& var_name) 
  {
    // assert that the id can be used to index into _factors
    DCHECK_EQ(variable.id(), num_factors());

    // Define a unary factor (the concept not the class) over the var to 
    // support a prior. 
    // NOTE the variable's prior distribution is always dense
    factor_type prior(factor_type::DENSE_TABLE, dense_table_t(variable));
    // using shift normalization, this is equivalent to uniform()
    prior.table()->zero();

    // NOTE the variable's prior distribution is always dense
    factor_type belief(factor_type::DENSE_TABLE, dense_table_t(variable));
    // using shift normalization, this is equivalent to uniform()
    belief.table()->zero(); 
    // assert that prior and belief have the same domain
    //ASSERT_TRUE(belief.domain() == prior.domain());

    // catalog the variable and associated metadata
    vertex_data_t vertex(prior, belief, true, var_name);
    _factors.push_back(vertex);
    //ASSERT_EQ(_factors[variable.id()].belief, vertex.belief);
    //ASSERT_EQ(_factors[variable.id()].potential, vertex.potential);

    return _factors.back();
  }

  factor_type uniform_factor_from_factor(const factor_type& factor) {
    factor_type belief = factor;
    // using shift normalization, this is equivalent to uniform()
    belief.table()->zero(); 

    return belief;
  }

// factor related methods
public:
  /** Add a new discrete factor to the factor graph */
  // REVIEW prevent a factor from being double added
  void add_factor(const table_base_t& factor, 
      const std::string& default_factor_name = "") 
  {
    size_t id = _unique_var_id++;
    DCHECK_NE(factor.ndims(), 0);
    // assert that the id can be used to index into _factors
    DCHECK_EQ(id, num_factors());

    std::string factor_name; 
    if(default_factor_name.empty()) {
      std::stringstream ss;
      ss << id;
      factor_name = ss.str();
    } else { 
      factor_name = default_factor_name;
    }

    // assert all variables have already been added to the graph
    logstream(LOG_INFO) << "ndims=" << factor.ndims() << " id=" << id 
        << " description='" << factor_name << "'" << std::endl;
    for(size_t i = 0; i < factor.ndims(); ++i) {
      logstream(LOG_INFO) << "  factor.var(" << i << ").id()=" << factor.var(i).id() << std::endl;
      DCHECK_LT(factor.var(i).id(), num_factors());
      DCHECK_EQ(factor.var(i).id(), get_variable(factor.var(i).id()).id());
    }

    factor_type node(factor);
    factor_type uniform = uniform_factor_from_factor(node);
    vertex_data_t vertex(node, uniform, false, factor_name);
    _factors.push_back(vertex);
    //ASSERT_EQ(_factors[id], vertex);
  }

// utils
public:
  size_t num_factors() const {
    return factors().size();
  }

  // FIXME O(n)
  size_t num_variables() const {
    size_t ndims = 0;
    for(typename std::vector<vertex_data_t >::const_iterator factor = factors().begin(); 
        factor != factors().end(); ++factor) {
      // any vertex that has a domain with only a single dimension is a variable ...
      if(factor->potential.table()->ndims() == 1) {
        ndims++;
      }
    }
    return ndims;
  }

  const std::string& name(const size_t id) const {
    DCHECK_LT(id, num_factors());
    return factors()[id].name;
  }

  /**
   * write a dot file which can be loaded into graphviz
   */
  void save_graph_summary(const std::string& filename) {
    DCHECK_EQ(_unique_var_id, num_factors());

    std::ofstream fout(filename.c_str());
    if(fout.is_open() == false) {
      std::cerr << "ERROR: " << filename << " not opened." << std::endl;
      return;
    }

    fout << "graph G {" << std::endl;
    fout << "layout=sfdp;" << std::endl;
    fout << "overlap=false;" << std::endl; 
    //fout << "sccmap;" << std::endl;
    fout << "K=2;" << std::endl;
    //fout << "clusterrank=local;" << std::endl;

    // Iterate all the factors and all the edges. NOTE all variables are also factors
    typename std::vector<vertex_data_t >::const_iterator factor;
    size_t factor_idx = 0;
    for(factor = factors().begin(); 
        factor != factors().end(); ++factor_idx, ++factor) 
    {
      //fout << "subgraph cluster_" << factor_idx << " {" << std::endl;
      //fout << "color=none;" << std::endl;
      // Iterate edges for a factor
      for(size_t i = 0; i < factor->potential.table()->ndims(); ++i) {
        variable_type variable = factor->potential.table()->var(i);
        // all variables are also factors, dont link a variable to itself
        if(variable.id() == factor_idx) continue; 

        fout << "\"" << factor->name << " {" << factor_idx << "}" 
               << "\" -- \"" << name(variable.id()) << "{" << variable.id() << "}" 
               << "\";" << std::endl; 
      }
      //fout << "}" << std::endl;
    }
    fout << "}" << std::endl;
  } // end of save_graph_summary

  /**
   * Construct a belief propagation graph from a factor graph
   */
  // NOTE could rewrite this function to construct the graph in parallel using 
  // the graphlab::distributed_control obj
  // TODO rename to finalize_distributed_graph()
  // REVIEW because bound, damping and regularization are now an attribute of a 
  // factor, perhaps they should be set in add_factor()
  void make_bp_graph(graph_type_t& graph, 
      double bound, double damping, double regularization=0.0) 
  {
    DCHECK_NE(num_factors(), 0);
    DCHECK_EQ(_unique_var_id, num_factors());

    // TODO clear the graph
    
    graphlab::timer timer; 

    if (graph.dc().procid() == 0) 
    {
      // Add all the factors and all the edges. NOTE all variables are also factors
      typename std::vector<vertex_data_t >::iterator factor = factors().begin();
      typename std::vector<vertex_data_t >::const_iterator end = factors().end();
      size_t factor_idx = 0;
      for( ; factor != end; ++factor_idx, ++factor) {
        // Add the factor to the graph
        factor->BOUND = bound; 
        factor->DAMPING = damping;
        factor->REGULARIZATION = regularization;
        graph.add_vertex(factor_idx, *factor);
        // TODO does the order in which i add variables and factors matter? 

        // Attach all the edges
        for(size_t i = 0; i < factor->potential.table()->ndims(); ++i) {
          variable_type variable = factor->potential.table()->var(i);
          // all variables are also factors, dont link a variable to itself
          if(variable.id() == factor_idx) continue; 

          // NOTE from graph::add_edge() - An edge can only be added if both the  
          // source and target vertex id's are already in the graph. Duplicate  
          // edges are not supported and may result in undefined behavior
          // NOTE messages are always dense
          dense_table_t msg(variable);
          msg.zero(); // using shift normalization, this is equivalent to uniform()
          graph.add_edge(factor_idx, variable.id(), edge_data_t(msg));
        }
      }
      graph.dc().cout() << "Loading graph. Finished in " 
          << timer.current_time() << std::endl;
    }
    timer.start();
    graph.finalize();
    graph.dc().cout() << "Finalizing graph. Finished in " 
        << timer.current_time() << std::endl;

    graph.dc().cout() 
      << "================ "
      << "Graph statistics on proc " << graph.dc().procid() << " of " << graph.dc().numprocs()
      << " ==============="
      << "\n Num vertices: " << graph.num_vertices()
      << "\n Num edges: " << graph.num_edges()
      << "\n Num replica: " << graph.num_replicas()
      << "\n Replica to vertex ratio: " 
      << float(graph.num_replicas())/graph.num_vertices()
      << "\n --------------------------------------------" 
      << "\n Num local own vertices: " << graph.num_local_own_vertices()
      << "\n Num local vertices: " << graph.num_local_vertices()
      << "\n Replica to own ratio: " 
      << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
      << "\n Num local edges: " << graph.num_local_edges()
      //<< "\n Begin edge id: " << graph.global_eid(0)
      << "\n Edge balance ratio: " 
      << float(graph.num_local_edges())/graph.num_edges()
      << std::endl;
  } // end of make_bp_graph

  /** 
   * Update the value of the variables' beliefs given the propegated 
   * graphlab distributed-graph. 
   */ 
  // TODO rename to pull_beliefs
  void pull_beliefs_for_variables(graph_type_t& graph) {
    typedef dense_table_t const *const const_ptr;

    // aggregate (reduce) the variable verticies into a vector
    aggregate_vertex_data agg = 
        graph.template map_reduce_vertices<aggregate_vertex_data>(
          aggregate_vertex_data()
        ); // wow

    for(unsigned i=0; i < agg.size(); ++i) {
      const vertex_data_t& ovdata = agg.agg[i];
      DCHECK_EQ(ovdata.belief.table()->ndims(), 1);

      const_ptr other = dynamic_cast<const_ptr>(ovdata.belief.table());
      if(other == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }

      vertex_data_t& vdata = factors()[other->args().var(0).id()];
      DCHECK_EQ(vdata.belief.table()->ndims(), 1);

      const_ptr tbl = dynamic_cast<const_ptr>(vdata.belief.table());
      if(tbl == NULL) {
        std::cout << "ERROR: std::bad_cast" << std::endl;
        // REVIEW should probably raise an exception
        ASSERT_TRUE(false);
      }
      ASSERT_EQ(tbl->args().var(0).id(), other->args().var(0).id());

      vdata.belief = ovdata.belief;
    }
  }

private:
  struct aggregate_vertex_data {
    std::vector<vertex_data_t> agg;

    aggregate_vertex_data() { }

    aggregate_vertex_data operator()(const typename graph_type_t::vertex_type& vertex) {
      const vertex_data_t& ovdata = vertex.data();

      aggregate_vertex_data out;
      // variables are always dense and 1D, although not all 1D vertices are variables...
      if(ovdata.isVariable == true) {
        out.agg.push_back(ovdata);
      }

      return out;
    }

    aggregate_vertex_data& operator+=(const aggregate_vertex_data& other) {
      agg.insert(agg.end(), other.agg.begin(), other.agg.end());
      return *this;
    }

    unsigned size() { return agg.size(); }

    void save(graphlab::oarchive& arc) const { arc << agg; }
    void load(graphlab::iarchive& arc) { arc >> agg; }
  };

public:
  void print_variable(const variable_type& var, const std::vector<double>& labels,
      std::ostream& out = std::cout) 
  {
    const dense_table_t& table = belief_for_variable(var);
    DCHECK_EQ(table.size(), labels.size());

    //std::ios::fmtflags f(out.flags()); // i cant believe this is how you do this

    out << "var_" << var.id() << ": " 
              << factors()[var.id()].name << std::endl;
    out << std::setw(8) << "index" << std::setw(16) << "logP" << std::setw(16) << "label" << std::endl;
    size_t end = table.size();
    for(size_t i=0; i < end; ++i) {
      out << std::setw(8) << i 
          << std::setw(16) << table.logP(i) 
          << std::setw(16) << labels[i] << std::endl;
    }

    //out.flags(f);
  }

// accessors
private:
  const std::vector<vertex_data_t>& factors() const { return _factors; }
  std::vector<vertex_data_t>& factors() { return _factors; }

private:
  // NOTE if ever multithreaded, this requires atomic access 
  size_t _unique_var_id;
  // REVIEW deep-copying data into std::vectors is slow 
  std::vector<vertex_data_t > _factors; 
};

} // end of namespace belief_prop

//#include <graphlab/macros_undef.hpp>
#endif // VSI_FACTOR_GRAPH_HPP

