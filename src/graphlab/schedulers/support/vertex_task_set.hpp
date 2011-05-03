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


/**
 * Maintain a set of tasks associated with a vertex.
 **/
#ifndef VERTEX_TASK_SET_HPP
#define VERTEX_TASK_SET_HPP

#include <queue>

#include <graphlab/tasks/update_task.hpp>
#include <graphlab/parallel/pthread_tools.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  template<typename Graph>
  class vertex_task_set {
  public:
    typedef Graph graph_type;
    typedef update_task<Graph> update_task_type;
    typedef typename update_task_type::update_function_type 
    update_function_type;

    typedef std::pair< update_function_type, double > fun_priority_pair;
    typedef std::vector< fun_priority_pair > vertex_fun_set;

  private:
    /**
     * taskset[vertexid] contains a list of updatefunctions pending
     * on this vertex
     */
    std::vector< vertex_fun_set > task_set;    
    
    /// The accompanying set of locks for the task_set
    std::vector<spinlock> locks;

    /**
     * Examine a task set for an update function returning true if
     * found.
     */
    bool vset_find(const vertex_fun_set& fun_set,
                   const update_function_type& fun,
                   size_t& ret_index) const {      
      for(size_t i = 0; i < fun_set.size(); ++i) {
        if(fun_set[i].first == fun) {
          ret_index = i;
          return true;
        }
      }
      return false;
    } // end of find

    /** Get the top element index  */
    size_t vset_top(const vertex_fun_set& fun_set) const {
      assert(!fun_set.empty());
      size_t top_index = 0;
      for(size_t index = 0; index < fun_set.size(); ++index) {
        if(fun_set[index].second > fun_set[top_index].second)
          top_index = index;
      }
      return top_index;
    } // end of top

    /**
     * get the total priority of a vset
     */
    double vset_total(const vertex_fun_set& fun_set) const {
      double total = 0;
      for(size_t i = 0; i < fun_set.size(); ++i) 
        total += fun_set[i].second;
      return total;
    } // end of total

    /**
     * Remove an element form a set quickly
     */
    void vset_remove(vertex_fun_set& fun_set,
                     size_t index) {
      assert(index < fun_set.size());
      // Here we do an O(1) removal
      // Place the last element at the current elements location
      fun_set[index] = fun_set[fun_set.size() - 1];
      // Erase the last element 
      fun_set.erase(fun_set.end() - 1);
    } // end of vset_remove
    
    
  public:
    /** Initialize the per vertex task set */
    vertex_task_set(size_t numvertices) :
      task_set(numvertices),
      locks(numvertices) { }

    /**
     * Resize the internal locks for a different graph
     */
    void resize(size_t num_vertices) {
      task_set.resize(num_vertices);
      locks.resize(num_vertices);
    }

    
    /** Find and remove the task from the set */
    bool remove(const update_task_type &task) {
      assert(task.vertex() < locks.size());
      bool success = false;
      // Grab the lock on the vertex
      locks[task.vertex()].lock();
      // Get the update functions associated with the vertex
      vertex_fun_set& fun_set = task_set[task.vertex()];
      // Try and find the task in the tasks in the local vertex "set"
      size_t index = 0;
      // only if the task was found is there a need to erase
      if(vset_find(fun_set, task.function(), index)) {
        // Remove the element at index
        vset_remove(fun_set, index);
        success = true;
      }
      // release the lock
      locks[task.vertex()].unlock();
      return success;
    }

    /** Add a task to the set returning false if the task was already
        present. Promote task to max(old priority, new priority) */
    bool add(const update_task_type& task, double priority) {
      assert(task.vertex() < locks.size());
      // Get the update functions associated with the vertex
      vertex_fun_set& fun_set = task_set[task.vertex()];
      // grab the lock on the vertex
      locks[task.vertex()].lock();
      // Try and find the task in the tasks in the local vertex "set"
      size_t index = 0;
      bool already_present = 
        vset_find(fun_set, task.function(), index);
      if(!already_present) {
        // If not found then this is the first add
        fun_set.push_back(std::make_pair(task.function(), priority));
      } else { 
        // it was already present so we should promote the priority to
        // the max of the two
        fun_set[index].second = std::max(fun_set[index].second, priority);
      }
      // release the lock
      locks[task.vertex()].unlock();
      // Return the result
      return !already_present;
    } // end of add task to set 
    
    
    /** Add a task to the set returning false if the task was already
        present. */
    bool add(const update_task_type& task) {
      assert(task.vertex() < locks.size());
      // Get the update functions associated with the vertex
      vertex_fun_set& fun_set = task_set[task.vertex()];
      // grab the lock on the vertex
      locks[task.vertex()].lock();
      // Try and find the task in the tasks in the local vertex "set"
      size_t index = 0;
      bool already_present = 
        vset_find(fun_set, task.function(), index);
      if(!already_present) {
        // If not found then this is the first add
        fun_set.push_back(std::make_pair(task.function(), 0.0));
      } 
      // release the lock
      locks[task.vertex()].unlock();
      // Return the result
      return !already_present;
    } // end of add task to set 


    
    /**
     * Look at the top task associated with the vertex.  If there are
     * no associated tasks then return false;
     */
    bool top(vertex_id_t vertex_id,
             update_task_type& ret_task,
             double& ret_priority) {
      assert(vertex_id < locks.size());
      // Flag to record if the unique add was successful
      bool successful(false);
      // grab the lock on the vertex
      locks[vertex_id].lock();
      // Get the vertex functions
      vertex_fun_set& fun_set = task_set[vertex_id];      
      // Determine if there is an element to remove
      if(!fun_set.empty()) {
        size_t top_index = vset_top(fun_set);
        // Get the last task
        ret_task  = update_task_type(vertex_id, fun_set[top_index].first);
        ret_priority = fun_set[top_index].second;
        successful = true;
      } else {
        successful = false;
      }
      // release the lock
      locks[vertex_id].unlock();
      // Return the result
      return successful;
    } // end of top


    /**
     * Get the top priority for a vertex
     */
    double top_priority(vertex_id_t vertex_id) {
      assert(vertex_id < locks.size());
      double priority(0);
      // grab the lock on the vertex
      locks[vertex_id].lock();
      // Get the vertex functions
      vertex_fun_set& fun_set = task_set[vertex_id];      
      // Determine if there is an element to remove
      if(!fun_set.empty()) {
        size_t top_index = vset_top(fun_set);
        priority = fun_set[top_index].second;        
      } 
      // release the lock
      locks[vertex_id].unlock();
      // Return the result
      return priority;
    } // get the top priority

    /**
     * Get the top priority for a vertex
     */
    double total_priority(vertex_id_t vertex_id) {
      assert(vertex_id < locks.size());
      double priority(0);
      // grab the lock on the vertex
      locks[vertex_id].lock();
      // Get the vertex functions
      vertex_fun_set& fun_set = task_set[vertex_id];      
      // Determine if there is an element to remove
      if(!fun_set.empty()) {
        priority = vset_total(fun_set);
      } 
      // release the lock
      locks[vertex_id].unlock();
      // Return the result
      return priority;
    } // get the total priority

    
    
    /**
     * Look at the top task associated with the vertex.  If there are
     * no associated tasks then return false;
     */
    bool pop(vertex_id_t vertex_id,
             update_task_type& ret_task,
             double& ret_priority) {
      assert(vertex_id < locks.size());
      // Flag to record if the unique add was successful
      bool successful(false);
      // grab the lock on the vertex
      locks[vertex_id].lock();
      // Get the vertex functions
      vertex_fun_set& fun_set = task_set[vertex_id];      
      // Determine if there is an element to remove
      if(!fun_set.empty()) {
        size_t top_index = vset_top(fun_set);
        // Get the top task
        ret_task  = update_task_type(vertex_id, fun_set[top_index].first);
        ret_priority = fun_set[top_index].second;
        // Same as top except we remove the top
        vset_remove(fun_set, top_index);        
        successful = true;
      } else {
        successful = false;
      }
      // release the lock
      locks[vertex_id].unlock();
      // Return the result
      return successful;
    } // end of top

    
    
    
    
  };

}
#include <graphlab/macros_undef.hpp>

#endif
