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

#ifndef GRAPHLAB_NET_UTIL_HPP
#define GRAPHLAB_NET_UTIL_HPP
#include <string>
#include <stdint.h>

namespace graphlab {
  /**
  * \ingroup util
  * Returns the first non-localhost ipv4 address 
  */
  uint32_t get_local_ip();

  /**
  * \ingroup util
  * Returns the first non-localhost ipv4 address as a standard dot delimited string
  */
  std::string get_local_ip_as_str();
  /** \ingroup util 
   * Find a free tcp port and binds it. Caller must release the port.
   * Returns a pair of [port, socket handle]
   */
  std::pair<size_t, int> get_free_tcp_port();
};

#endif

