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

#define JSONSL_STATE_GENERIC
#include <jsonsl.h>
#include <graphlab.hpp>
#include "rapidjson.hpp"

namespace graphlab {

  /**
   * JSON Message
   * Formats, sends, receives, and parses JSON messages.
   */
  class json_message {
  
  private:

    rapidjson::Document mdocument;
    jsonsl_t mjsn;
  
  public:
  
    typedef char byte;
    
    /** Maximum recursion depth for lexer */
    const static int MAX_LEVELS;
  
    /**
     * Create a new message for specified method, with the given updater state.
     */
    json_message(const std::string method = "", const std::string state = "");
    
    virtual ~json_message();
    
    bool feed(byte *data, std::size_t nread);
    
    friend std::ostream& operator<< (std::ostream &out, json_message &message);
  
  };

};

#endif /* #ifndef GRAPHLAB_JSON_MESSAGE_HPP */
