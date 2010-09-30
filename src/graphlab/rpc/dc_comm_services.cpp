#include <boost/iostreams/stream.hpp>
#include <sstream>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_comm_services.hpp>

#include <graphlab/serialization/serialization_includes.hpp>
namespace graphlab {
namespace dc_impl {
  
const char CHILD_TO_PARENT = 0;
const char PARENT_TO_CHILD = 1;

void dc_comm_services::__child_to_parent_barrier_trigger(procid_t source,
                                                         char barrierval) {
  barrier_mut.lock();
  // wrong barrier value
  if (barrierval != barrier_sense) {
    barrier_mut.unlock();
    return;
  }
  if (source == child[0]) {
    child_barrier[0] = barrier_sense;
  }
  else {
    child_barrier[1] = barrier_sense;
  }
  barrier_cond.signal();
  barrier_mut.unlock();
}

void dc_comm_services::__parent_to_child_barrier_release(char releaseval) {
  // send the release downwards
  if (child[0] < comm->numprocs()) {
//    rpc.fast_remote_call(child[0],
//                        &dc_services::__parent_to_child_barrier_release,
//                        releaseval);
  }
  if (child[1] < comm->numprocs()) {
//    rpc.fast_remote_call(child[1],
//                        &dc_services::__parent_to_child_barrier_release,
//                        releaseval);
  }
  barrier_mut.lock();
  barrier_release = releaseval;    
  barrier_cond.signal();
  barrier_mut.unlock();
}



void dc_comm_services::recv(procid_t source, char *c, size_t len) {
  boost::iostreams::stream<boost::iostreams::array_source> strm(c, len);
  iarchive iarc(strm);
  size_t instr;
  iarc >> instr;
  
  char releaseval;    
  
  switch(instr) {
    case CHILD_TO_PARENT:
      iarc >> releaseval;
      __child_to_parent_barrier_trigger(source, releaseval);
      break;
    case PARENT_TO_CHILD:
      iarc >> releaseval;
      __parent_to_child_barrier_release(releaseval);
      break;
    default:
      logger(LOG_FATAL, "received unexpected packet");
  }
}




void dc_comm_services::barrier() {
  // upward message
  char barrier_val = barrier_sense;      
  barrier_mut.lock();
  // create the child to parent message
  std::string childtoparent_message;
  {
    std::stringstream strm;
    oarchive oarc(strm);
    oarc << CHILD_TO_PARENT;
    oarc << barrier_sense;
    childtoparent_message = strm.str();
  }
  
  
  while(1) {
    if ((child_barrier[0] == barrier_sense || child[0] >= comm->numprocs()) &&
        (child_barrier[1] == barrier_sense || child[1] >= comm->numprocs())) {
      // flip the barrier sense
      barrier_sense = ! barrier_sense;
      // call child to parent in parent
      barrier_mut.unlock();
      break;
    }
    barrier_cond.wait(barrier_mut);
  }
  logger(LOG_INFO, "barrier phase 1 complete");
  // I am root. send the barrier releae downwards
  if (comm->procid() == 0) {
    barrier_release = barrier_val;
    if (child[0] < comm->numprocs()) {
//      rpc.fast_remote_call(child[0],
//                          &dc_services::__parent_to_child_barrier_release,
//                          barrier_val);
    }
    if (child[1] < comm->numprocs()) {
//      rpc.fast_remote_call(child[1],
//                          &dc_services::__parent_to_child_barrier_release,
//                          barrier_val);
    } 
  }
  // wait for the downward message releasing the barrier
  barrier_mut.lock();
  while(1) {
    if (barrier_release == barrier_val) break;
    // send the upward packet
    if (comm->procid() != 0) {
      comm->send(parent, childtoparent_message.c_str(), 
                 childtoparent_message.length());
    }
    // wait for 0.1s
    barrier_cond.timedwait_ns(barrier_mut, 100000000);
  }
  barrier_mut.unlock();

  logger(LOG_INFO, "barrier phase 2 complete");
}


} // dc_impl
} // graphlab