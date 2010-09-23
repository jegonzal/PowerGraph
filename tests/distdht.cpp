#include <graphlab/logger/logger.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_hash_table.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
using namespace graphlab;

struct data{
  size_t a;
  double b;
  int crap[32];
  void save(oarchive &oarc) const {
    oarc << a;
    oarc << b;
    serialize(oarc, crap, sizeof(int)*32);
  }
  void load(iarchive &iarc) {
    iarc >> a;
    iarc >> b;
    deserialize(iarc, crap, sizeof(int)*32);
  }
};

int main(int argc, char** argv) {
  size_t max32bit = 4294967296;
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  distributed_control dc(&argc, &argv);
  dc.init_message_processing(1);
  distributed_hash_table blobmap(dc, 128);
  distributed_hash_table blobmap2(dc, 128);
  dc.barrier();

  // test basic insertion
  if (dc.procid() == 0) {
    // write
    for (size_t i = 0; i < 1000; ++i) {
      data val;
      val.a = i; val.b = i;
      blobmap.set(i,val);
      val.a++; val.b++;
      blobmap2.set(i,val);
    }
    for (size_t i = max32bit; i < max32bit+1000; ++i) {
      data val;
      val.a = i; val.b = i;
      blobmap.set(i,val);
      val.a++; val.b++;
      blobmap2.set(i,val);
    }
    //read
    for (size_t i = 0; i < 1000; ++i) {
      data val, cacheval;
      blobmap.get(i,val);
      blobmap.get_cached(i,cacheval);
      ASSERT_EQ(val.a, i);
      ASSERT_EQ(val.b, i);
      ASSERT_EQ(cacheval.a, i);
      ASSERT_EQ(cacheval.b, i);
    }
    for (size_t i = max32bit; i < max32bit+1000; ++i) {
      data val, cacheval;
      blobmap.get(i,val);
      blobmap.get_cached(i,cacheval);
      ASSERT_EQ(val.a, i);
      ASSERT_EQ(val.b, i);
      ASSERT_EQ(cacheval.a, i);
      ASSERT_EQ(cacheval.b, i);
    }
  }

  // test invalidation
  if (dc.procid() == 0) {
    for (size_t i = 0; i < 1000; ++i) {
      blobmap.invalidate(i);
    }
    for (size_t i = 0; i < 1000; ++i) {
      data val;
      ASSERT_TRUE(blobmap.get_cached(i,val));
      ASSERT_EQ(val.a, i);
      ASSERT_EQ(val.b, i);
    }
  }
  dc.barrier();
  // test a different processor
  if (dc.procid() != 1) {
    for (size_t i = 0; i < 1000; ++i) {
      data val;
      
      ASSERT_TRUE(blobmap.get_cached(i,val));
      ASSERT_EQ(val.a, i); ASSERT_EQ(val.b, i);
      ASSERT_TRUE(blobmap.get(max32bit+i,val));
      ASSERT_EQ(val.a, max32bit+i); ASSERT_EQ(val.b, max32bit+i);
    }
  }
  dc.barrier();

  // test caching
    // test a different processor
  for (size_t i = 0; i < 128; ++i) {
    data val;
    data val2;
    blobmap.get_cached(i, val);
    blobmap2.get_cached(i, val2);
    ASSERT_EQ(val.a, i);
    ASSERT_EQ(val.a+1, val2.a);
  }
  for (size_t i = 0; i < 128; ++i) {
    data val;
    data val2;
    blobmap.get_cached(i, val);
    blobmap2.get_cached(i, val2);
    ASSERT_EQ(val.a, i);
    ASSERT_EQ(val.a+1, val2.a);
  }
  dc.barrier();
  logstream(LOG_INFO) << "Gets: " << blobmap.num_gets() << "\tMisses: " << blobmap.num_misses() << std::endl;
  logstream(LOG_INFO) << "Gets2: " << blobmap2.num_gets() << "\tMisses2: " << blobmap2.num_misses() << std::endl;
  logstream(LOG_INFO) << "Cache Size: " << blobmap.cache_size() << std::endl;
  logstream(LOG_INFO) << "Cache Miss Rate: " << blobmap.cache_miss_rate() << std::endl;

// test caching 2
  // read the data, change it on another processor
  // second time I read it from the cache it should stay the same
  // then after I invalidate it and read it again, it should be correct
  // set initial conditions. invalidate on 0, set the data
  if (dc.numprocs() > 1) {
    if (dc.procid() == 0) {
      blobmap.invalidate(1);
    }
    else {
      data val;
      val.a = 1000; val.b = 1000;
      blobmap.set(1,val);
    }
    dc.barrier();
    // read it, it should return 1000
    data val;
    if (dc.procid() == 0) {
      blobmap.get_cached(1, val);
      ASSERT_EQ(val.a, 1000);
    }

    dc.barrier();
    // change it on a different processor
    if (dc.procid() == 1) {
      data val;
      val.a = 1001; val.b = 1001;
      blobmap.set(1,val);
    }
    dc.barrier();
    blobmap.invalidate(1);
    // read it again. It should be correct now
    if (dc.procid() == 0) {
      blobmap.get_cached(1, val);
      ASSERT_EQ(val.a, 1001);
    }
  }
  else {
    logstream(LOG_WARNING) << "Skipping cache correctness test" << std::endl;
  }
  dc.barrier();
  // bandwidth
  if (dc.procid() == 1) {
    logger(LOG_INFO, "Processor 1 Sleeping for 3 seconds to test background poll");
    sleep(3);
  }
  std::vector<size_t> nums;
  const size_t numel = 100000;
  for (size_t i = 0;i < numel; ++i) nums.push_back(rand() % 1000);
  logstream(LOG_INFO) << dc.procid() <<": Bandwidth test started" << std::endl;
  timer ti;
  ti.start();
  for (size_t i = 0;i < nums.size(); ++i) {
    data val;
    blobmap.get(nums[i], val);
  }
  double curtime = ti.current_time();
  double mbps = double(numel * sizeof(data)) / curtime / (1024 * 1024);
  logstream(LOG_INFO) << dc.procid() << ": "<< mbps <<" MBps" << std::endl;
  logstream(LOG_INFO) << numel << ": "<< " requests completed in " << curtime << "seconds" << std::endl;
  logstream(LOG_INFO) << curtime/numel * 1000.0 << " ms RTT" << std::endl;
  dc.barrier(); 
}
