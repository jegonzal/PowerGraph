#ifndef GRAPHLAB_RPC_THREAD_LOCAL_SEND_BUFFER_HPP
#define GRAPHLAB_RPC_THREAD_LOCAL_SEND_BUFFER_HPP
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/rpc/dc_compile_parameters.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/util/dense_bitset.hpp>
namespace graphlab {
class distributed_control;


namespace dc_impl {
struct thread_local_buffer {
  std::vector<std::vector<std::pair<char*, size_t> > > outbuf;
  std::vector<mutex> outbuf_locks;
  std::vector<size_t> bytes_sent;

  fixed_dense_bitset<RPC_MAX_N_PROCS> contended;

  std::vector<mutex> archive_locks;
  std::vector<oarchive> current_archive;
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
    archive_locks[target].lock();
    // need a new archive, or existing one at risk of being resized
    if (current_archive[target].buf == NULL) {
      current_archive[target].buf = (char*)malloc(INITIAL_BUFFER_SIZE);
      current_archive[target].off = 0;
      current_archive[target].len = INITIAL_BUFFER_SIZE;
    }
    prev_acquire_archive_size = current_archive[target].off;
    return &current_archive[target];
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
      bytes_sent[target] += current_archive[target].off - prev_acquire_archive_size - sizeof(packet_hdr);
      inc_calls_sent(target);
    }

    if (current_archive[target].off >= INITIAL_BUFFER_SIZE || contended.get(target)) {
      // shift the buffer into outbuf
      char* ptr = current_archive[target].buf;
      size_t len = current_archive[target].off;
      current_archive[target].buf = NULL; 
      current_archive[target].off = 0;
      archive_locks[target].unlock();

      outbuf_locks[target].lock();
      outbuf[target].push_back(std::make_pair(ptr, len));
      outbuf_locks[target].unlock();
      contended.clear_bit(target);
    } else {
      archive_locks[target].unlock();
    }
  }

  inline void write(procid_t target, char* c, size_t len, bool do_not_count_bytes_sent) {
    if (!do_not_count_bytes_sent) {
      bytes_sent[target] += len;
      inc_calls_sent(target);
    }
    // make sure that messsages sent before this write are sent before this write
    if (current_archive[target].off) {
      archive_locks[target].lock();
      outbuf_locks[target].lock();
      outbuf[target].push_back(std::make_pair(current_archive[target].buf,
                                              current_archive[target].off));
      outbuf_locks[target].unlock();
      current_archive[target].buf = NULL; 
      current_archive[target].off = 0;
      contended.clear_bit(target);
      archive_locks[target].unlock();
    }
    outbuf_locks[target].lock();
    outbuf[target].push_back(std::make_pair(c, len));
    outbuf_locks[target].unlock();
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
    if (outbuf[target].size() > 0) {
      outbuf_locks[target].lock();
      std::swap(ret, outbuf[target]);
      outbuf_locks[target].unlock();
    }
    if (current_archive[target].off > 0 ) {
      if (archive_locks[target].try_lock()) {
        char* ptr = current_archive[target].buf;
        size_t len = current_archive[target].off;
        if (len > 0) {
          current_archive[target].buf = NULL;
          current_archive[target].off = 0;
        }
        archive_locks[target].unlock();
        if (len > 0) ret.push_back(std::make_pair(ptr, len));
        contended.clear_bit(target);
      } else {
        contended.set_bit(target);
      } 
    } 
    return ret;
  }

  void inc_calls_sent(procid_t target);
};
}
}
#endif
