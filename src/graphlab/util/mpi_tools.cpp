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

#include <graphlab/util/net_util.hpp>
#include <graphlab/util/mpi_tools.hpp>

namespace graphlab {
namespace mpi_tools {
  
  

void get_master_ranks(std::set<size_t>& master_ranks) {
  uint32_t local_ip = get_local_ip();
  std::vector<uint32_t> all_ips;
  all_gather(local_ip, all_ips);
  std::set<uint32_t> visited_ips;
  master_ranks.clear();
  for(size_t i = 0; i < all_ips.size(); ++i) {
    if(visited_ips.count(all_ips[i]) == 0) {
      visited_ips.insert(all_ips[i]);
      master_ranks.insert(i);
    }
  }
}

} 
}
