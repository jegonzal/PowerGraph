#include <graphlab/rpc/async_consensus.hpp>

namespace graphlab {
async_consensus::async_consensus(distributed_control &dc,
                                const dc_impl::dc_dist_object_base *attach)
          :rmi(dc, this), attachedobj(attach),
           last_calls_sent(0), last_calls_received(0),
           waiting_on_done(false), cancelled(false),
           complete(false), hastoken(dc.procid() == 0) {

  cur_token.total_calls_sent = 0;
  cur_token.total_calls_received = 0;
  cur_token.last_change = rmi.numprocs() - 1;
}

bool async_consensus::done() {
  mut.lock();
  
  while (!complete) {
    waiting_on_done = true;
    if (hastoken) pass_the_token();
    cond.wait(mut);
    // I got woken up!
    // if complete or cancelled, leave now
    if (complete || cancelled) {
      cancelled = false;
      waiting_on_done = false;
      break;
    }
    // otherwise this is a spurious wake up.
    // continue waiting
  }
  mut.unlock();
  return complete;
}

void async_consensus::cancel() {
  mut.lock();
  if (waiting_on_done) {
    cancelled = true;
    cond.signal();
  }
  mut.unlock();
}

void async_consensus::receive_the_token(token &tok) {
  mut.lock();
  // save the token
  hastoken = true;
  cur_token = tok;
  // if I am waiting on done, pass the token.
  if (waiting_on_done) {
    pass_the_token();
  }
  mut.unlock();
}

void async_consensus::pass_the_token() {
  // note that this function does not acquire the token lock
  // the caller must acquire it 
  assert(hastoken);
  // first check if we are done
  if (cur_token.last_change == rmi.procid() && cur_token.total_calls_received == cur_token.total_calls_sent) {
    // we have completed a loop around!
    // broadcast a completion
    for (size_t i = 0;i < rmi.numprocs(); ++i) {
      if (i != rmi.procid()) {
        rmi.control_call(i,
                        &async_consensus::consensus);
      }
    }
    // myself
    complete = true;
    cond.signal();
  }
  else {
    // update the token
    size_t callsrecv;
    size_t callssent;
    
    if (attachedobj) {
      callsrecv = attachedobj->calls_received();
      callssent = attachedobj->calls_sent();
    }
    else {
      callsrecv = rmi.dc().calls_received();
      callssent = rmi.dc().calls_sent();
    }

    if (callssent != last_calls_sent ||
        callsrecv != last_calls_received) {
      cur_token.total_calls_sent += callssent - last_calls_sent;
      cur_token.total_calls_received += callsrecv - last_calls_received;
      cur_token.last_change = rmi.procid();
    }
    //std::cout << "Sending token: (" << cur_token.total_calls_sent << ", " << cur_token.total_calls_received << ")" << std::endl;

    last_calls_sent = callssent;
    last_calls_received = callsrecv;
    // send it along.
    hastoken = false;
    rmi.control_call((rmi.procid() + 1) % rmi.numprocs(),
                    &async_consensus::receive_the_token,
                    cur_token);
  }
}


void async_consensus::consensus() {
  mut.lock();
  complete = true;
  cond.signal();
  mut.unlock();
}
}