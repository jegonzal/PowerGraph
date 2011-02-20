#ifndef GRAPHLAB_DISTRIBUTED_GLSHARED_BASE
#define GRAPHLAB_DISTRIBUTED_GLSHARED_BASE

#include <vector>
#include <boost/shared_ptr.hpp>
#include <graphlab/shared_data/glshared.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc_types.hpp>

namespace graphlab {

class distributed_glshared_base;
class distributed_glshared_manager;

namespace distgl_impl {
  /** the actual registration table is stored as a static variable
  * in a function to avoid the "global variable construction order" problem.
  */
  std::vector<distributed_glshared_base*>& get_global_dist_glshared_registry();

  /**
  * Registers a distributed_glshared object into the global registry
  */
  void register_dist_glshared(distributed_glshared_base* glsharedobj);
} // namespace distgl_impl

/**
 distributed glshared base class. Allows the glshared manager
 to manage it
*/
class distributed_glshared_base: public glshared_base {
 protected:
  mutable distributed_glshared_manager* manager;
  size_t id;
  mutable bool invalidated;
  
  friend class distributed_glshared_manager;

 public:
 
  typedef glshared_base::apply_function_type apply_function_type;

  distributed_glshared_base() : manager(NULL) {
    distgl_impl::register_dist_glshared(this);
    invalidated = true;
  }
  
  virtual any get_any() const = 0;
  virtual void set_any(const any&) = 0;
  virtual void apply(apply_function_type fun,
                     const any& srcd) = 0;                       
  virtual bool is_unique() const = 0;
  
  virtual void save(oarchive &oarc) const = 0;
  virtual void load(iarchive &iarc) = 0;
  virtual const char* type_name() const = 0;

  // return the machine on which operations perform the fastest
  virtual procid_t preferred_machine() const = 0;
};

}

#endif
