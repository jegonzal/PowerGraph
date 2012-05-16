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
 * @file dispatcher.hpp
 * @author Jiunn Haur Lim <jiunnhal@cmu.edu>
 */
 
#ifndef GRAPHLAB_DISPATCHER_HPP
#define GRAPHLAB_DISPATCHER_HPP

#include <graphlab.hpp>

////////////////////////////// CGI UPDATE FUNCTOR //////////////////////////////

namespace graphlab {

  struct cgi_vertex {
    std::string state;
  };
  
  struct cgi_edge {
    std::string state;
  };
  
  class json_return;
  
  /** dispatcher update functor */  
  class dispatcher_update : 
    public graphlab::iupdate_functor<graph<cgi_vertex, cgi_edge>, dispatcher_update> {
  public:
    typedef graphlab::graph<cgi_vertex, cgi_edge> graph_type;
  private:
    std::string mstate;
  public:
    dispatcher_update(const std::string& state="");
    dispatcher_update(const dispatcher_update& other);
    inline void operator+=(const dispatcher_update& other);
    void operator()(icontext_type& context);
  private:
    /** schedule updates on vertices based on return values */
    void schedule(icontext_type& context, const json_return& result);
  }; // end of dispatcher update functor
  
};

////////////////////////////////////////////////////////////////////////////////

#endif /* #ifndef GRAPHLAB_DISPATCHER_HPP */
