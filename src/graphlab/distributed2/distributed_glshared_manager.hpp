#ifndef DISTRIBUTED_GLSHARED_MANAGER_HPP
#define DISTRIBUTED_GLSHARED_MANAGER_HPP
#include <vector>
#include <string>
#include <map>
#include <boost/bind.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/coherent_dht.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
// forward declaration of glshared
template <typename T>
class distributed_glshared;
/**
The manager which provides the distributed glshared variables
with synchronization and communication capabilities
*/
class distributed_glshared_manager {
 private:
  dc_dist_object<distributed_glshared_manager> rmi;
  // a list of all the objects attached
  std::vector<distributed_glshared_base*> glsharedobjs;
  // a reverse map 
  std::map<distributed_glshared_base*, size_t> objrevmap;
  // the DHT used to synchronize everyone
  coherent_dht<size_t, std::string> dht;
  
 public:
  distributed_glshared_manager(distributed_control &dc);  
  /*
  completes an atomic exchange of an entry
  */
  std::string exchange(size_t entry, const std::string &val) ;
  
  void invalidate(size_t entry, const std::string& value,bool incache) ;


  /**
  Synchronize variable with index i.
  Call
  */
  void write_synchronize(size_t entry, bool async = false);

  void write_synchronize(distributed_glshared_base* obj, bool async = false);
  
  
  void read_synchronize(size_t entry, bool async = false);
  
  void read_synchronize(distributed_glshared_base* obj, bool async = false);
};

} // namespace graphlab
#endif  
