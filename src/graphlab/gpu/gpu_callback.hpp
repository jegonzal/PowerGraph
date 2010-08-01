#ifndef GRAPHLAB_GPU_CALLBACK_HPP
#define GRAPHLAB_GPU_CALLBACK_HPP

#include <vector>


#include <graphlab/graph/graph.hpp>


namespace graphlab {

  class gpu_callback {
  public:

    /**
     * Adds a task to execute the update function on the vertex with
     * the given priority.
     */
    void add_task(vertex_id_t vertex, double priority) { }
                         
  };

}; //end graphlab namespace

#endif
