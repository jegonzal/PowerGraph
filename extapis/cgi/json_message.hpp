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
   * JSON Message
   * Formats, sends, receives, and parses JSON messages.
   */
  class json_message {
  
  private:

    rapidjson::Document mdocument;
  
  public:
  
    typedef char byte;
    
    const static byte VERTEX = 0x01;
    const static byte EDGES = 0x02;
  
    /**
     * Creates a new message for specified method, with the given updater state.
     * @param[in]   method    method to invoke
     * @param[in]   state     state of invocation owner (the updater)
     */
    json_message(const std::string method = "", const std::string state = "");
    virtual ~json_message();
    
    ////////////////////////////// ACCESSORS ///////////////////////////////////
    
    const char *updater();
    const char *vertex();
    
    /**
     * Adds a context parameter to the invocation
     * @param[in]   context   invocation context
     * @param[in]   flags     vertex only, or vertex and edges
     */
    void add_context(dispatcher_update::icontext_type& context, byte flags);
    
    ///////////////////////////////// I/O //////////////////////////////////////
    
    /**
     * Parses input data as a JSON document
     * @param[in]   data      input data
     * @param[in]   bytes     size of input data
     */
    void parse(byte *data, std::size_t bytes);

    /**
     * Prints message to output stream
     */
    friend std::ostream& operator<< (std::ostream &out, json_message &message);
  
  };

};

#endif /* #ifndef GRAPHLAB_JSON_MESSAGE_HPP */
