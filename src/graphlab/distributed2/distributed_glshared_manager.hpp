#ifndef DISTRIBUTED_GLSHARED_MANAGER_HPP
#define DISTRIBUTED_GLSHARED_MANAGER_HPP
#include <vector>
#include <string>
#include <map>
#include <boost/bind.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/coherent_dht.hpp>
#include <graphlab/distributed2/distributed_glshared_base.hpp>
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
  typedef distributed_glshared_base::apply_function_type apply_function_type;
  
  distributed_glshared_manager(distributed_control &dc);  
  ~distributed_glshared_manager();
  /*
  completes an atomic exchange of an entry
  */
  std::string exchange(size_t entry, const std::string &val) ;

  template <typename T>
  void apply_from_remote(size_t entry, size_t fun, const any& param) {
    apply<T>(entry, reinterpret_cast<apply_function_type>(fun), param);
  }
  
  template <typename T>
  void apply(size_t entry, apply_function_type fun, const any& param) {
    if (dht.owning_machine(entry) == rmi.procid()) {
      std::string& valref = dht.begin_critical_section(entry);
      // deserialize the entry from the DHT
      std::stringstream strm(valref);
      iarchive iarc(strm);
      T curval;
      iarc >> curval;
      // put it in an any
      any curany = curval;
      // call the function
      fun(curany, param);
      // serialize it and put it back
      std::stringstream ostrm;
      oarchive oarc(ostrm);
      oarc << (curany.as<T>());
      // store back the value
      valref = ostrm.str();
      dht.end_critical_section(entry);
      dht.push_changes(entry, false, procid_t(-1));
    }
    else {
      rmi.remote_request(dht.owning_machine(entry),
                          &distributed_glshared_manager::apply_from_remote<T>,
                          entry,
                          reinterpret_cast<size_t>(fun),
                          param);
    }
    dht.invalidate(entry);
  }
  
  void invalidate(size_t entry, const std::string& value,bool incache) ;
  
  inline procid_t preferred_machine(size_t entry) {
    return dht.owning_machine(entry);
  }

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
