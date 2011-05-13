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

