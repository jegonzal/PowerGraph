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
   * JSON Schedule (extracted from return JSON)
   */
  class json_schedule {
  public:
    enum neighbors_enum {ALL, SOME, NONE, IN_NEIGHBORS, OUT_NEIGHBORS};
  private:
    std::string mupdater;
    neighbors_enum mtargets;
    std::vector<unsigned> mvertices;
  public:
    
    json_schedule(neighbors_enum targets,
      const std::string &updater="");
    json_schedule(neighbors_enum targets,
      const std::string &updater,
      const rapidjson::Value &vertices);
    
    /**
     * Enum representing vertices to schedule
     */
    neighbors_enum targets() const;
    
    /**
     * Returns updater state for new updater (or 'self')
     */
    const std::string &updater() const;
    
    /**
     * Returns array of vertex IDs if these were specified in the JSON
     */
    const std::vector<unsigned> &vertices() const;
    
  };


  /**
   * JSON Message (Abstract)
   * Formats, sends, receives, and parses JSON messages.
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
   */
  class json_invocation : public json_message {
  public:
  
    const static byte VERTEX = 0x01;
    const static byte EDGES = 0x02;
  
    /**
     * Creates a new message for specified method, with the given updater state.
     * @param[in]   method    method to invoke
     * @param[in]   state     state of invocation owner (the updater)
     */
    json_invocation(const std::string& method = "", const std::string& state = "");
    virtual ~json_invocation();
    
    /**
     * Adds a context parameter to the invocation
     * @param[in]   context   invocation context
     * @param[in]   flags     vertex only, or vertex and edges
     */
    void add_context(dispatcher_update::icontext_type& context, byte flags);
    
  private:
    void add_in_edges(dispatcher_update::icontext_type& context, rapidjson::Value& parent);
    void add_out_edges(dispatcher_update::icontext_type& context, rapidjson::Value& parent);
    rapidjson::Value& create_edge(rapidjson::Value& edgev, const dispatcher_update::icontext_type& context, const dispatcher_update::graph_type::edge_type& edge);
    rapidjson::Value& create_vertex(rapidjson::Value& vertexv, const dispatcher_update::graph_type::vertex_id_type vertex_id, const dispatcher_update::graph_type::vertex_data_type& vertex_data);
    
  };
  
  /**
   * JSON Return
   * A specific JSON message for return values.
   */
  class json_return : public json_message {
  public:
    
    virtual ~json_return();
    
    /**
     * @return updater state from return JSON, or null if not available.
     */
    const char *updater() const;
    
    /**
     * @return vertex state from return JSON, or null if not available.
     */
    const char *vertex() const;
    
    /**
     * @return schedule, extracted from return JSON
     */
    const json_schedule schedule() const;
    
  };

};

#endif /* #ifndef GRAPHLAB_JSON_MESSAGE_HPP */
