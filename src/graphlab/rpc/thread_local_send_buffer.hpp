#ifndef GRAPHLAB_RPC_THREAD_LOCAL_SEND_BUFFER_HPP
#define GRAPHLAB_RPC_THREAD_LOCAL_SEND_BUFFER_HPP
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/rpc/dc_compile_parameters.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
namespace graphlab {
class distributed_control;


namespace dc_impl {
struct thread_local_buffer {
  std::vector<std::vector<oarchive> > oarc;
  std::vector<mutex> locks;
  std::vector<size_t> bytes_sent;
  size_t prev_acquire_archive_size;
  procid_t procid;
  distributed_control* dc;

  thread_local_buffer();
  ~thread_local_buffer();

  /**
   * Must be called from within the thread owning this buffer.
   * Acquires a buffer to write to
   */
  inline oarchive* acquire(procid_t target) {
    locks[target].lock();
    // need a new archive, or existing one at risk of being resized
    if (oarc[target].size() == 0 || oarc[target].back().off >= INITIAL_BUFFER_SIZE) {
      // allocate a new one
      oarc[target].push_back(oarchive());
    }
    if (oarc[target].back().buf == NULL) {
      oarc[target].back().buf = (char*)malloc(INITIAL_BUFFER_SIZE);
      oarc[target].back().off = 0;
      oarc[target].back().len = INITIAL_BUFFER_SIZE;
    }
    prev_acquire_archive_size = oarc[target].back().off;
    return &oarc[target].back();
  }

  inline size_t get_bytes_sent(procid_t target) {
    return bytes_sent[target];
  }
  /**
   * Must be called from within the thread owning this buffer.
   * Releases a buffer previously acquired with acquire
   */
  inline void release(procid_t target, bool do_not_count_bytes_sent) {
    if (!do_not_count_bytes_sent) {
      bytes_sent[target] += oarc[target].back().off - prev_acquire_archive_size - sizeof(packet_hdr);
      inc_calls_sent(target);
    }
    locks[target].unlock();
  }

  inline void write(procid_t target, char* c, size_t len, bool do_not_count_bytes_sent) {
    if (!do_not_count_bytes_sent) {
      bytes_sent[target] += len;
      inc_calls_sent(target);
    }
    locks[target].lock();
    oarc[target].push_back(oarchive());
    oarc[target].back().buf = c;
    oarc[target].back().off = len;
    oarc[target].back().len = len;
    locks[target].unlock();
  }

  /**
   * Must be called from within the thread owning this buffer.
   * Flushes the buffer to the sender. This should really only be used
   * when the thread is dying since this incurs a large performance penalty by
   * locking up the sender.
   */
  void push_flush();


  /**
   * Can be called anywhere.
   * Flushes the buffer to the sender. This function blocks until all 
   * buffers have been flushed. Equivalent to calling distributed_control::flush()
   */
  void pull_flush();

  /**
   * Can be called anywhere.
   * Flushes the buffer to the sender. This function blocks until all 
   * buffers have been flushed. Equivalent to calling distributed_control::flush()
   */
  void pull_flush(procid_t p);

  /**
   * Can be called anywhere.
   * Flushes the buffer to the sender. This function requests a flush to happen
   * soon. Equivalent to calling distributed_control::flush()
   */
  void pull_flush_soon();




  /**
   * Can be called anywhere.
   * Flushes the buffer to the sender. This function requests a flush to happen
   * soon. Equivalent to calling distributed_control::flush()
   */
  void pull_flush_soon(procid_t p);



  /**
   * Extracts the buffer going to a given target.
   */
  inline std::vector<std::pair<char*, size_t> > extract(procid_t target) {
    std::vector<std::pair<char*, size_t> > ret;
    std::vector<oarchive> arcs;
    if (oarc[target].size() > 0) {
      locks[target].lock();
      std::swap(arcs, oarc[target]);
      locks[target].unlock();

      ret.resize(arcs.size());
      for (size_t i = 0;i < arcs.size(); ++i) {
        ret[i].first = arcs[i].buf;
        ret[i].second = arcs[i].off;
      }
    }
    return ret; 
  }

  void inc_calls_sent(procid_t target);
};
}
}
#endif
