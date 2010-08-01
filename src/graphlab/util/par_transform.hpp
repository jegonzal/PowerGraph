#ifndef GRAPHLAB_PAR_TRANSFORM_HPP
#define GRAPHLAB_PAR_TRANSFORM_HPP

#include <iterator>
#include <graphlab/core.hpp>

namespace graphlab {

  namespace partransform_impl {

    template <typename CoreType, typename TransformF, typename T>
    void transform_update_function(typename CoreType::types::iscope& scope,
                                   typename CoreType::types::icallback& scheduler,
                                   typename CoreType::types::ishared_data* data_manager) {
      // get the function

      any a = data_manager->get_constant(0);
      TransformF* transformfunction = (TransformF*)(a.as<void*>());
      (*transformfunction)(*(scope.vertex_data()));
    }

  }

  template <typename Iterator, typename TransformF>
  void par_transform(Iterator start, Iterator end, const TransformF& f, size_t ncpus = 0){


    typedef typename Iterator::value_type value_type;
    typedef core<value_type*, char> core_type;
  
    core_type glcore;
    // create alot of vertices
    glcore.set_scope_type("null");
    glcore.set_engine_type("async");
    glcore.set_scheduler_type("sweep");
    if (ncpus > 0) glcore.set_ncpus(ncpus);

    Iterator i = start;
    while (i!=end) {
      glcore.graph().add_vertex(&(*i));
      ++i;
    }
  
    glcore.shared_data().set_constant(0, (void*)(&f));
    glcore.add_task_to_all(partransform_impl::
                           transform_update_function<core_type, TransformF, value_type>, 1.0);
    glcore.start();
  }



  template <typename Iterator, typename TransformF>
  void par_transform(Iterator start, Iterator end, const TransformF& f,
                     const engine_options& eopts){


    typedef typename Iterator::value_type value_type;
    typedef core<value_type*, char> core_type;
  
    core_type glcore;
    // create alot of vertices
    glcore.set_engine_options(eopts);
    

    Iterator i = start;
    while (i!=end) {
      glcore.graph().add_vertex(&(*i));
      ++i;
    }
  
    glcore.shared_data().set_constant(0, (void*)(&f));
    glcore.add_task_to_all(partransform_impl::
                           transform_update_function<core_type, TransformF, value_type>, 1.0);
    glcore.start();
  }


}
#endif
