#ifndef GRAPHLAB_DC_DIST_OBJECT_BASE_HPP
#define GRAPHLAB_DC_DIST_OBJECT_BASE_HPP
#include <vector>
#include <graphlab/rpc/dc_internal_types.hpp>
namespace graphlab {

namespace dc_impl {
/**
 * \ingroup rpc
Provides an interface for extracting and updating counters from dc_dist_objects
*/
class dc_dist_object_base{
 public:
  /// Increment the number of calls sent from this object
  virtual void inc_calls_sent(procid_t source) = 0;
  /// Increment the number of calls received by this object
  virtual void inc_calls_received(procid_t dest) = 0;
  
  /// Increment the number of bytes sent from this object
  virtual void inc_bytes_sent(procid_t target, size_t bytes) = 0;

  /// Return the number of calls received by this object
  virtual size_t calls_received() const = 0;
  /// Return the number of calls sent from this object
  virtual size_t calls_sent() const = 0;
};

}
}

#endif
