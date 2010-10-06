#include <graphlab/distributed/metrics/distributed_metrics.hpp>
#include <graphlab/distributed/distributed_control.hpp> 

namespace graphlab {

mutex singletonlock;

distributed_metrics * distributed_metrics_receive_target = NULL;

distributed_metrics * distributed_metrics::instance(distributed_control * dc) {
    singletonlock.lock();
    if (distributed_metrics_receive_target == NULL) {
         distributed_metrics_receive_target = new distributed_metrics(dc);
    }
    singletonlock.unlock();
    return distributed_metrics_receive_target;
}
 
  // Remote call
void distributed_metrics::remote_set_value(distributed_control& _dc, size_t source, void *ptr, size_t len,
        std::string key, double value) {
     instance(&_dc)->set_value(source, key, value);
}

void distributed_metrics::remote_start_bandwidth_test(distributed_control& _dc, size_t source, void *ptr, size_t len) {
    instance(&_dc)->start_bandwidth_test(source);
}

void distributed_metrics::receive_bandwidth_test_payload(distributed_control& _dc, size_t source, void *ptr, size_t len) {
    instance(&_dc)->end_bandwith_test(source, len);
}

  

}