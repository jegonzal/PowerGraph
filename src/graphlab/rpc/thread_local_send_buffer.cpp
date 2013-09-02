#include <graphlab/rpc/thread_local_send_buffer.hpp>
#include <graphlab/rpc/dc.hpp>

namespace graphlab {
namespace dc_impl {

thread_local_buffer::thread_local_buffer() {
  // allocate the buffers
  dc = distributed_control::get_instance();
  size_t nprocs = dc->numprocs(); 

  outbuf.resize(nprocs); 
  current_archive.resize(nprocs); 

  outbuf_locks.resize(nprocs);
  archive_locks.resize(nprocs);

  contended.clear();

  bytes_sent.resize(nprocs, 0);
  dc->register_send_buffer(this);
  procid = dc->procid();
}


thread_local_buffer::~thread_local_buffer() {
  dc->unregister_send_buffer(this);
  push_flush();
  // deallocate the buffers
  for (size_t i = 0; i < current_archive.size(); ++i) {
    if (current_archive[i].buf) {
      free(current_archive[i].buf);
      current_archive[i].buf = NULL;
    }
  }

  for (size_t i = 0; i < outbuf.size(); ++i) {
    for (size_t j = 0; j < outbuf[i].size(); ++j) {
      free(outbuf[i][j].first);
    }
  }
}

void thread_local_buffer::inc_calls_sent(procid_t target) {
  dc->inc_calls_sent(target);
}


void thread_local_buffer::push_flush() {
  for (size_t i = 0; i < outbuf.size(); ++i) {
    std::vector<std::pair<char*, size_t> >  buf = extract(i);
    for (size_t j = 0;j < buf.size(); ++j) {
      dc->write_to_buffer(i, buf[j].first, buf[j].second);
    }
    if (buf.size()) dc->flush_soon(i);
  }
}


void thread_local_buffer::pull_flush() {
  dc->flush();
}

void thread_local_buffer::pull_flush(procid_t p) {
  dc->flush(p);
}


void thread_local_buffer::pull_flush_soon() {
  dc->flush_soon();
}


void thread_local_buffer::pull_flush_soon(procid_t p) {
  dc->flush_soon(p);
}

oarchive* thread_local_buffer::acquire(procid_t target) {
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

void thread_local_buffer::release(procid_t target, bool do_not_count_bytes_sent) {
  if (!do_not_count_bytes_sent) {
    bytes_sent[target] += current_archive[target].off - prev_acquire_archive_size - sizeof(packet_hdr);
    inc_calls_sent(target);
  }

  if (current_archive[target].off >= FULL_BUFFER_SIZE_LIMIT || contended.get(target)) {
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
    if (outbuf[target].size() >= NUM_FULL_BUFFER_LIMIT) {
      pull_flush_soon(target);
    }
  } else {
    archive_locks[target].unlock();
  }
}


void thread_local_buffer::write(procid_t target, char* c, size_t len, 
                                bool do_not_count_bytes_sent) {
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
  if (outbuf[target].size() >= NUM_FULL_BUFFER_LIMIT) {
    pull_flush_soon(target);
  }
}



std::vector<std::pair<char*, size_t> > thread_local_buffer::extract(procid_t target) {
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



} // dc_impl
} // graphlab
