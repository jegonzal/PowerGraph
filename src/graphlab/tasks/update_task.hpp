/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_UPDATE_TASK_HPP
#define GRAPHLAB_UPDATE_TASK_HPP

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {
  // Predecleration 
  template<typename Graph> class icallback;
  template<typename Graph> class ishared_data;

  
  template<typename Graph>
  class update_task {
  public:

    typedef Graph graph_type;
    typedef iscope<Graph> iscope_type;
    typedef icallback<Graph> callback_type;
    typedef ishared_data<Graph> ishared_data_type;

    /// This is the standard vertex update function
    typedef void(*update_function_type)(iscope_type& scope,
                                        callback_type& callback);


    struct  hash_functor {
      size_t operator()(const update_task& t) const { return t.hash();  }
    };

  private:
    vertex_id_t vertexid;  
    update_function_type func;

  public:
    explicit update_task(vertex_id_t vertexid = -1, 
                         update_function_type func = NULL) :
      vertexid(vertexid), func(func) { }

    ~update_task() {}
 
    // Move into engine to simplify typing of this object
    // void execute(scope_type& scope,
    //              ischeduler_callback &scheduler,
    //              const shared_data_manager* data_manager) {
    //   assert(func != NULL);
    //   func(scope, scheduler, data_manager);
    // }
    
    vertex_id_t vertex() const {
      return vertexid;
    }
    
    update_function_type function() const {
      return func;
    }

    /// Returns true if tasks are identical
    bool operator==(const update_task &i) const{
      return vertexid == i.vertexid && func == i.func ;
    }
    
    /// comparator
    bool operator<(const update_task &i) const{
      return (vertexid < i.vertexid) || 
        (vertexid == i.vertexid && func < i.func);
    }
    
    size_t hash() const {
      // TODO: this was arbitrarily decided. need something better here
      return vertexid ^ (size_t)(void*)(func);
    }
    
    void save(oarchive &oarc) const {
      oarc << vertexid;
      oarc << reinterpret_cast<size_t>(func);
    }
    
    void load(iarchive &iarc)  {
      iarc >> vertexid;
      size_t funcptr;
      iarc >> funcptr;
      func = reinterpret_cast<update_function_type>(funcptr);
    }
  
  };
   
  
}
#include <graphlab/macros_undef.hpp>
#endif
