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