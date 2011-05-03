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

#ifndef BINARY_VERTEX_TASK_SET
#define BINARY_VERTEX_TASK_SET


/**
  * \Simple task set that only keeps record of whether 
  *  function was scheduled for given vertex.
  **/


#include <vector>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/parallel/atomic.hpp>

namespace graphlab {
  
  template<typename Graph>
  class binary_vertex_task_set {    
  public:
    typedef Graph graph_type;
    typedef update_task<Graph> update_task_type;
    typedef typename update_task_type::update_function_type 
    update_function_type;



    std::vector<size_t> vertexbits;
    std::vector<update_function_type> updatefuncs;
    size_t num_of_updatefunctions;
    spinlock updflock;

    const static size_t MAX_UPDATEFUNCTIONS;

  public:    
    
    binary_vertex_task_set(size_t numvertices) :
      vertexbits(numvertices, 0),
      updatefuncs(MAX_UPDATEFUNCTIONS, 0),
      num_of_updatefunctions(0) {
      // std::cout << "Constructed for " << numvertices
      //           << " vertices" << std::endl;
    }   
   
    
    bool get(const update_task_type& task) {
      update_function_type func = task.function();
      vertex_id_t vid = task.vertex();
      size_t mask = get_update_func_mask(func);
      return ((vertexbits[vid] & mask) != 0);
    }
    
    bool add(const update_task_type& task) {
      update_function_type func = task.function();
      vertex_id_t vid = task.vertex();
      size_t mask = get_update_func_mask(func);
      if ((vertexbits[vid] & mask) == 0) {
        size_t before = __sync_fetch_and_or(&vertexbits[vid], mask);
        return (before & mask) == 0;	// did someone else already set the bit 
      } 
      return false;
    }
    
    void remove(const update_task_type& task) {
      update_function_type func = task.function();
      vertex_id_t vid = task.vertex();
      size_t mask = get_update_func_mask(func);

      // Set bit for this function to 0
      __sync_fetch_and_and(&vertexbits[vid], ~mask);
    }
    
    size_t pop_all_tasks(size_t vid) {
      size_t val = __sync_lock_test_and_set(&vertexbits[vid], 0);
      return val;
    }

  private:
    size_t get_update_func_mask(update_function_type upf) {
      for(size_t i = 0; i < num_of_updatefunctions; i++ ){
        if (updatefuncs[i] == upf) return 1<<i;
      }
      
      // Ok not found, have to make it. Now we need to lock.
      updflock.lock();
      
      // Check once more
      for(size_t i = 0; i < num_of_updatefunctions; i++ ){
        // Update functions resized by a separate thread while
        // grabbing the lock
        if (updatefuncs[i] == upf) {
          updflock.unlock();
          return 1<<i;
        }
      }
      assert(num_of_updatefunctions < MAX_UPDATEFUNCTIONS);
      updatefuncs[num_of_updatefunctions] = upf;
      int newid = 1 << num_of_updatefunctions;
      num_of_updatefunctions++;
      updflock.unlock();
      return newid;
    }
    
  };

  
  template <typename Graph>
  const size_t binary_vertex_task_set<Graph>::MAX_UPDATEFUNCTIONS =
    sizeof(size_t) * 8; 

}

#endif
