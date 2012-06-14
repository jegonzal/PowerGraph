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
 * @file cgi.hpp
 * Convenience methods for writing a C++ GraphLab CGI program. You don't really
 * need this if you are using your own json parser (i.e. not using rapidjson)
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */

#ifndef GRAPHLAB_CGI_HPP
#define GRAPHLAB_CGI_HPP

#include <sstream>
#include <iostream>
#include "../../rapidjson.hpp"

/** initial read buffer size */
#define INITIAL_LENGTH 256

namespace graphlab {

  /**
   * Interface for all CGI handlers. Every method takes a JSON document and returns
   * a JSON document.
   */
  class icgi_handler {
    public:
      virtual ~icgi_handler() {}
      virtual void gather_edges(rapidjson::Document& invocation, rapidjson::Document& return_json) = 0;
      virtual void scatter_edges(rapidjson::Document& invocation, rapidjson::Document& return_json) = 0;
      virtual void gather(rapidjson::Document& invocation, rapidjson::Document& return_json) = 0;
      virtual void merge(rapidjson::Document& invocation, rapidjson::Document& return_json) = 0;
      virtual void apply(rapidjson::Document& invocation, rapidjson::Document& return_json) = 0;
      virtual void scatter(rapidjson::Document& invocation, rapidjson::Document& return_json) = 0;
      virtual void init_edge(rapidjson::Document& invocation, rapidjson::Document& return_json){}
      virtual void init_vertex(rapidjson::Document& invocation, rapidjson::Document& return_json){}
  };
  
  /**
   * CGI client. Parses STDIN, invokes callbacks on handler, and sends results
   * to STDOUT.
   */
  class cgi {
  private:
    icgi_handler& handler;
    
  public:
    
    /**
     * Creates a new CGI client.
     * @param[in]   cgi_handler     implements icgi_handler and contains most of
     *                              the program logic
     */
    cgi(icgi_handler& cgi_handler) : handler(cgi_handler){}
    
    /**
     * Start listening for JSON invocation. This method blocks until "exit" is
     * received.
     */
    void listen(){
      
      std::string line;
      std::size_t length = 0;
      std::size_t current_length = INITIAL_LENGTH;
      char *buffer = new char[current_length];
  
      // loop until exit is received
      while (true){
    
        // TODO: assume NULL-terminated string for now
        
        // read length, and expand buffer if necessary
        std::cin >> length; std::getline(std::cin, line);
        if (length + 1 > current_length){
          current_length = length + 1;
          delete[] buffer;
          buffer = new char[current_length];
        }
    
        // read message
        std::cin.read(buffer, length);
        buffer[length] = '\0';    // terminate string w. null
        rapidjson::StringBuffer return_buffer;
    
        // invoke corresponding callback method
        std::cerr << buffer << std::endl;
        const char *return_json = handle_invocation(buffer, return_buffer);
        if (!return_json) break;
    
        // return
        std::cout << strlen(return_json) << "\n";
        std::cout << return_json << std::flush;
    
      }
  
      delete[] buffer;
      
    }
    
  private:
  
    const char *handle_invocation(const char *buffer, rapidjson::StringBuffer& return_buffer){

      if (NULL == buffer) return NULL;
    
      rapidjson::Document invocation;
      if (invocation.Parse<0>(buffer).HasParseError()){
        std::ostringstream error("JSON syntax error: ");
        error << buffer;
        throw error.str();
      }
      
      if (!invocation.HasMember("method"))
        throw "Missing method field.";

      if (!strcmp(invocation["method"].GetString(), "exit"))
        return NULL;
      
      rapidjson::Document return_json;
      return_json.SetObject();
      
      const char *method = invocation["method"].GetString();
      
      // invoke appropriate method
      if (!strcmp(method, "gather_edges"))
        handler.gather_edges(invocation, return_json);
      else if (!strcmp(method, "scatter_edges"))
        handler.scatter_edges(invocation, return_json);
      else if (!strcmp(method, "gather"))
        handler.gather(invocation, return_json);
      else if (!strcmp(method, "merge"))
        handler.merge(invocation, return_json);
      else if (!strcmp(method, "apply"))
        handler.apply(invocation, return_json);
      else if (!strcmp(method, "scatter"))
        handler.scatter(invocation, return_json);
      else if (!strcmp(method, "init_edge"))
        handler.init_edge(invocation, return_json);
      else if (!strcmp(method, "init_vertex"))
        handler.init_vertex(invocation, return_json);
      else {
        std::ostringstream error("Unknown method: ");
        error << method;
        throw error.str();
      }
      
      // write return json to stdout
      return_buffer.Clear();
      rapidjson::Writer<rapidjson::StringBuffer> writer(return_buffer);
      return_json.Accept(writer);
      return return_buffer.GetString();
      
    }
    
  };
  
};

#endif
