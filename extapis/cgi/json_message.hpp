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
 * @file json_message.hpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#ifndef GRAPHLAB_JSON_MESSAGE_HPP
#define GRAPHLAB_JSON_MESSAGE_HPP

#include <graphlab.hpp>

#include "rapidjson.hpp"
#include "dispatcher.hpp"

namespace graphlab {

  /**
   * JSON Message (Abstract)
   * Formats, sends, receives, and parses JSON messages.
   * @internal
   */
  class json_message {
  
  protected:
  
    /** JSON Document -> think of it as the root */
    rapidjson::Document mdocument;
    rapidjson::Document::AllocatorType& mallocator;

    json_message();
  
  public:

    typedef char byte;
    
    ///////////////////////////////// I/O //////////////////////////////////////
    
    /**
     * Parses input data as a JSON document
     * @param[in]   data      input data
     * @param[in]   bytes     size of input data
     */
    void parse(const byte *data, std::size_t bytes);

    /**
     * Prints message to output stream
     */
    friend std::ostream& operator<< (std::ostream &out, json_message &message);
  
  };
  
  /**
   * JSON Invocation
   * A specific JSON message for invoking methods.
   * @internal
   */
  class json_invocation : public json_message {
  public:
  
    /**
     * Creates a new message for specified method, with the given updater state.
     * @param[in]   method    method to invoke
     * @param[in]   state     state of invocation owner (the updater)
     */
    json_invocation(const std::string& method, const std::string& state);
    json_invocation(const std::string& method);
    virtual ~json_invocation();
    
    /**
     * Adds a vertex parameter to the invocation
     * @param[in]   vertex    vertex
     */ 
    void add_vertex(const dispatcher::vertex_type& vertex);
    
    /**
     * Adds an edge parameter to the invocation
     * @param[in]   edge      edge
     */
    void add_edge(dispatcher::edge_type& edge);
    
    /**
     * Adds other gather state to the invocation
     * @param[in]   other     gather state
     */
    void add_other(const std::string& other);
    
    /**
     * Adds gather total to the invocation
     * @params[in]  gather    gather total
     */
    void add_gather(const dispatcher::gather_type& gather_total);
    
  private:
    
    /**
     * Initializes vertexv with state and id from vertex
     * @internal
     */
    rapidjson::Value& create_vertex(rapidjson::Value& vertexv, const dispatcher::graph_type::vertex_type& vertex);
    
    /**
     * Initializes edgev with source, target, and state from edge
     * @internal
     */
    rapidjson::Value& create_edge(rapidjson::Value& edgev, dispatcher::graph_type::edge_type& edge);
    
    /**
     * Ensures that the "params" element exists in the document
     */
    void ensure_params_exist();
    
  };
  
  /**
   * JSON Return
   * A specific type of json_message for return values.
   * @internal
   */
  class json_return : public json_message {
  public:
    
    virtual ~json_return();
    
    /** @return program state from return JSON, or null if not available. */
    const char *program() const;
    
    /** @return vertex state from return JSON, or null if not available. */
    const char *vertex() const;
    
    /** @return edge state from return JSON, or null if not available. */
    const char *edge() const;
    
    /** @return edge direction from return JSON, or null if not available. */ 
    const edge_dir_type edge_dir() const;
    
    const char *result() const;
    const char *signal() const;
    
    
  };

};

#endif /* #ifndef GRAPHLAB_JSON_MESSAGE_HPP */
