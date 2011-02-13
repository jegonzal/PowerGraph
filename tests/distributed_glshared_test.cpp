#include <string>
#include <algorithm>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/distributed2/distributed_glshared.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>

using namespace graphlab;


distributed_glshared<std::string> astring;
distributed_glshared<size_t> anumber;

void atomic_test(distributed_control &dc) {
  std::vector<size_t> seennumbers;
  for (size_t i = 0; i < 100; ++i) {
    size_t num = 1 + i + size_t(100) * dc.procid();
    anumber.exchange(num);
    seennumbers.push_back(num);
  }
  std::vector<std::vector<size_t> > gather;
  gather.resize(dc.numprocs());
  gather[dc.procid()] = seennumbers;
  
  dc.gather(gather, 0);
  
  if (dc.procid() == 0) {
    // combine them all into one big vector
    std::vector<size_t> ret;
    for (size_t i = 0;i < gather.size(); ++i) {
      std::copy(gather[i].begin(), gather[i].end(), std::back_inserter(ret));
    }
    
    // ok. I will be missing one element. That will be the element currently in the number
    ret.push_back(anumber.get_val());
    
    std::sort(ret.begin(), ret.end());
    
    std::vector<size_t> uniqueret = ret;
    std::vector<size_t>::iterator newend = std::unique(uniqueret.begin(), uniqueret.end());
    uniqueret.resize(newend - uniqueret.begin());
    ASSERT_EQ(ret.size(), uniqueret.size());
    ASSERT_EQ(ret[0], 0);
    ASSERT_EQ(ret[ret.size() - 1], (size_t)(100) * (dc.numprocs() - 1) + 100);
  }
}

int main(int argc, char** argv) {
  dc_init_param param;
  assert(init_param_from_env(param));
  global_logger().set_log_level(LOG_DEBUG);

  distributed_control dc(param);
  distributed_glshared_manager glmanager(dc);
  if (dc.procid() == 0) {
    anumber.set(10);
    astring.set("hello");
  }
  dc.full_barrier();
  ASSERT_EQ(anumber.get_val(), 10);
  ASSERT_EQ(astring.get_val(), std::string("hello"));
  dc.barrier();
  if (dc.procid() == 1) {
    anumber.set(0);
    astring.set("pika");
  }
  dc.full_barrier();
  ASSERT_EQ(anumber.get_val(), 0);
  ASSERT_EQ(astring.get_val(), std::string("pika"));
  dc.barrier();
  atomic_test(dc);
  dc.barrier();
}
