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


#ifndef VSI_BP_VERTEX_PROGRAM_HPP
#define VSI_BP_VERTEX_PROGRAM_HPP

/**
 * This file defines the max-sum vertex program for beleif propagation.
 *
 * \author Scott Richardson     10/2012
 */

#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include "table_factor.hpp"
#include "bp_graph_data.h"


namespace belief_prop {


/** 
 * The gather_type for the vertex program needs to compute *= in place
 * of += so we create a new type which computes *= for +=.
 */
template<size_t MAX_DIM>
class factor_product {
  typedef graphlab::table_factor<MAX_DIM>    factor_type;

public:
  factor_type factor;
  // REVIEW deep copying the factor around could get expensive.
  //   after profiling, it seems most of this is either negligible or gets
  //   optimized out
  factor_product(const factor_type& factor = factor_type()) : 
      factor(factor) { }
  factor_product& operator+=(const factor_product& other) {
    DCHECK_EQ(factor.table()->numel(), other.factor.table()->numel());
    factor *= other.factor;
    return *this;
  }
  void save(graphlab::oarchive& arc) const { arc << factor; }
  void load(graphlab::iarchive& arc) { arc >> factor; }
}; // end of struct factor product


/** 
 * Belief Propagation Vertex Program. As implemented, this "program" 
 * performs the max-sum algorithm. GraphLab runs this program at 
 * every vertex.
 * 
 * \author Scott Richardson
 */
template<size_t MAX_DIM>
class bp_vertex_program : 
    public graphlab::ivertex_program< typename graph_type<MAX_DIM>::type, 
                                      factor_product<MAX_DIM>,
                                      graphlab::messages::sum_priority >
{
  // unfortunately this is necessary...from C++ Standard 14.6.2/3:
  // "In the definition of a class template or a member of a class template, if a
  // base class of the class template depends on a template-parameter, the base 
  // class scope is not examined during unqualified name lookup either at the 
  // point of definition of the class template or member or during an instantiation 
  // of the class template or member.
  typedef graphlab::ivertex_program< typename graph_type<MAX_DIM>::type, 
                                     factor_product<MAX_DIM>,
                                     graphlab::messages::sum_priority > ivertex_program_t;
  // NOTE there is a bug in GCC < 4.7 which prevents these using declarations from 
  // compiling (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=14258)
  //using typename ivertex_program_t::edge_dir_type;
  //using typename ivertex_program_t::vertex_type;
  //using typename ivertex_program_t::edge_type;
  typedef typename ivertex_program_t::edge_dir_type edge_dir_type;
  typedef typename ivertex_program_t::vertex_type   vertex_type;
  typedef typename ivertex_program_t::edge_type     edge_type;

  typedef vertex_data<MAX_DIM>             vertex_data_t;
  typedef edge_data<MAX_DIM>               edge_data_t;
  typedef graphlab::table_factor<MAX_DIM>  factor_type; // vertex_data_t::factor_type
  typedef graphlab::dense_table<MAX_DIM>   msg_type;    // edge_data_t::msg_type

public:
  //using typename ivertex_program_t::gather_type;
  //using typename ivertex_program_t::icontext_type;
  typedef typename ivertex_program_t::gather_type gather_type;
  typedef typename ivertex_program_t::icontext_type icontext_type;

public:
  bp_vertex_program() { }

  /**
   * Since we are handling edge direction ourselves, we will use all edges for 
   * gather and scatter
   */
  edge_dir_type gather_edges(icontext_type& context,
                       const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /**
   * Update the old message to be the new message and collect the
   * message value.
   */
  gather_type gather(icontext_type& context, 
               const vertex_type& vertex, edge_type& edge) const {
    // NOTE While gather() technically has a non-const reference to the  
    // source and target vertex data through edge.source() and edge.target(), 
    // it should not modify them. The data on the edge (accessible through 
    // edge.data()) is modifiable however.
    const vertex_data_t&  vdata = vertex.data();
    const vertex_type     other_vertex = get_other_vertex(edge, vertex);
    const vertex_data_t&  ovdata = other_vertex.data();
    edge_data_t& edata = edge.data();
    logstream(LOG_DEBUG) << "in bp_vertex_program::gather(): compute message to '" 
        << vdata.name << "' from vertex '" << ovdata.name << "'" << std::endl;

    // Update the old message with the value of the new message. We
    // then receive the old message during gather and then compute the
    // "cavity" during scatter (again using the old message).
    edata.update_old(other_vertex.id(), vertex.id());
    msg_type& msg = edata.old_message(other_vertex.id(), vertex.id());
    logstream(LOG_DEBUG) << "edata=" << msg << std::endl;

    gather_type& rep = repmat(msg, vdata);
    logstream(LOG_EVERYTHING) << "repmat-ed msg=" << rep.factor << std::endl;

    logstream(LOG_DEBUG) << "END bp_vertex_program::gather()" << std::endl; 
    return rep;
  }; // end of gather function

  /**
   * Multiply message product by node potential and update the belief.
   */
  void apply(icontext_type& context, vertex_type& vertex, 
       const gather_type& total) {
    if(vertex.num_in_edges() + vertex.num_out_edges() == 0) return; 

    // factor_type knows it is in log space (so it adds)
    vertex_data_t& vdata = vertex.data();

    logstream(LOG_DEBUG) << "in bp_vertex_program::apply(): vertex = '" 
        << vdata.name << "'" << std::endl;

    //vdata.belief = vdata.potential * total.factor;
    vdata.belief = vdata.potential;
    //vdata.belief.table()->copy_onto(*(vdata.potential.table())); // should be faster...
    vdata.belief *= total.factor;
    logstream(LOG_EVERYTHING) << "vdata.potential=" << vdata.potential << std::endl;
    logstream(LOG_EVERYTHING) << "total.factor=" << total.factor << std::endl;
    logstream(LOG_EVERYTHING) << "vdata.belief=vdata.potential * total.factor = " 
        << vdata.belief << std::endl;
    if(vdata.isVariable == true) {
      logstream(LOG_INFO) << "belief-prop variable state = '" 
          << vdata.name << "' "
          << vdata.belief << std::endl;
    }

    DCHECK_GT(vdata.belief.table()->numel(), 0);
    // Rescale the belief to ensure numerical stability. (This is
    // essentially normalization in log-space.)
    // REVIEW is this needed to match belief_prop
    //vdata.belief.table()->shift_normalize();
    //logstream(LOG_DEBUG) << "vdata.shift_normalized=" << vdata.belief << std::endl;
    logstream(LOG_DEBUG) << "END bp_vertex_program::apply()" << std::endl;
  }; // end of apply

  /**
   * Since we are handling edge direction ourselves, we will use all edges for 
   * gather and scatter
   */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  /**
   * Compute new message value for each edge.
   */
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {
    const vertex_data_t&  vdata = vertex.data();
    const vertex_type     other_vertex = get_other_vertex(edge, vertex);
    const vertex_data_t&  ovdata = other_vertex.data();
    edge_data_t& edata = edge.data();
    logstream(LOG_DEBUG)
        << "in bp_vertex_program::scatter(): compute message from '" << vdata.name 
        << "' to vertex '" << ovdata.name << "'" << std::endl;

    // construct the cavity
    //factor_type cavity = vdata.belief;
    factor_type& cavity = init_factor(vdata.belief).factor; // initilizes the cache factor
    logstream(LOG_EVERYTHING) << "factor=" << cavity << std::endl;
    const msg_type& incoming_message = edata.old_message(other_vertex.id(), vertex.id());
    cavity /= incoming_message;
    logstream(LOG_DEBUG) << "incoming_message=" << incoming_message << std::endl;
    logstream(LOG_EVERYTHING) << "cavity=" << cavity << std::endl;
    //cavity.table()->normalize();

    // compute the new outgoing message
    const msg_type& old_message = 
        edata.old_message(vertex.id(), other_vertex.id());
    msg_type& new_message = 
        edata.message(vertex.id(), other_vertex.id());
    DCHECK_NE(&new_message, &old_message);

    // max-product algorithm 
    cavity.table()->MAP(new_message);
    // sum-product algorithm
    //cavity.table()->marginalize(new_message);


    if(vdata.REGULARIZATION > 0.0) {
      // shift normalize
      new_message.shift_normalize(); 
      //logstream(LOG_DEBUG) << "normalized=" << new_message << std::endl;
      // regularize
      logstream(LOG_DEBUG) << "regularization_value=" << vdata.REGULARIZATION << std::endl;
      msg_type reg(new_message.domain());
      reg.uniform(1.0);
      new_message.damp(reg, vdata.REGULARIZATION);
      logstream(LOG_DEBUG) << "regularized=" << new_message << std::endl;
    }

    // shift normalize
    new_message.shift_normalize(); 
    logstream(LOG_DEBUG) << "normalized=" << new_message << std::endl;
    msg_type raw_message(new_message);

    // dampen
    new_message.damp(old_message, vdata.DAMPING);
    logstream(LOG_DEBUG) << "damped=" << new_message << std::endl;

    // Compute message residual
    //const double residual = new_message.l1_diff(old_message);
    // dec: The literature seems to indicate that the l_inf norm is better, so lets try that.
    const double residual = new_message.linf_diff(old_message);
    logstream(LOG_INFO) << "residual=" << residual << std::endl;

    context.clear_gather_cache(other_vertex);
    // to prevent drift, we may want to only update the new_message   
    // if the residual is greater than the BOUND because we only 
    // signal the neighboring node when this is true. however this
    // doesnt work as well as expected...
    logstream(LOG_INFO) << "belief-prop message from '" 
        << vdata.name << "' to vertex '" << ovdata.name << "':" 
        << " raw=" << new_message // NOTE newlines here can break atomicity...
        << " damped=" << raw_message 
        << std::endl;
    // Schedule the adjacent vertex
    if(residual > vdata.BOUND) {
      context.signal(other_vertex, residual);
    }
//    else {
//      new_message = old_message; 
//    } 

    logstream(LOG_DEBUG) << "END bp_vertex_program::scatter()" << std::endl;
  }; // end of scatter

  /** Save the values to a binary archive */
  // NOTE no need to serialize the contents of cache, although it is not POD, 
  // so we must serialize something
  void save(graphlab::oarchive& arc) const { arc << gather_type(); }

  /** Read the values from a binary archive */
  void load(graphlab::iarchive& arc) { arc >> cache; }  

private:
  /**
   * Initilize 'gather_type cache'---a factor that is re-used to avoid costly 
   * data-structure construction.
   */
  // REVIEW not sure if this would work in async mode
  gather_type& init_factor(const factor_type& other) const {
    if(cache.factor.table_storage() == factor_type::nil) {
      cache.factor = other;
    } else {
      cache.factor.table()->copy_onto(*(other.table()));
    }

    return cache;
  }
  /** 
   * Return msg copied (broadcasted) across the domain defined by vdata
   */
  // you cant multiply (or add) two edge-messages because the intersection of
  // their domains is null (cf. dense_table::operator*() => dense_table::logP(asg) => 
  // discrete_assignment::restrict() fails), so i repmat the message to cover the domain
  gather_type& repmat(const msg_type& msg, const vertex_data_t& vdata) const {
    gather_type& ones = init_factor(vdata.belief);
    //ones = vdata.belief;
    ones.factor.table()->zero();
    // NOTE implicit broadcasting
    // NOTE factor_type knows it is in log space (so it adds)
    ones.factor *= msg;

    return ones;
  }
  /**
   * Return the other vertex
   */
  const vertex_type get_other_vertex(edge_type& edge, 
                               const vertex_type& vertex) const {
    return vertex.id() == edge.source().id() ? edge.target() : edge.source();
  }

private:
  // keeping a cache here has very good storage requirements as opposed to in vector_data. 
  // NOTE no need to serialize.
  mutable gather_type cache;
}; // end of class bp_vertex_program


} // end of namespace belief_prop

#endif // VSI_BP_VERTEX_PROGRAM_HPP
