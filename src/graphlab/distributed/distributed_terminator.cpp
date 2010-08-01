#include <vector>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_terminator.hpp>

namespace graphlab {
std::vector<distributed_terminator*> distributed_terminator::distributed_term_target;

distributed_terminator::distributed_terminator(distributed_control &dc) :dc(dc) {
  curtoken.totaltransmitted = 0;
  curtoken.totalreceived = 0;
  curtoken.terminatenow = false;

  my_total_transmitted = 0;
  my_total_received = 0;
  delta_transmitted.value = 0;
  delta_received.value = 0;
  quit = false;
  use_control_packets = false;
  if (dc.procid() == 0)  {
    hastoken = true;
  }
  else {
    hastoken = false;
  }
  termid = distributed_term_target.size();
  distributed_term_target.push_back(this);
}



void distributed_terminator::distributed_terminator_handler(distributed_control &dc,
                                    procid_t source,
                                    void * ptr,    //unused
                                    size_t len,   // unused
                                    size_t id,
                                    distributed_termination_token token) {
  // ok I just received the new token!
  // check if it is the same as the one I received last looparound
  // and if we can terminate
  if (token.terminatenow) {
    // note that the ordering here is important
    bool oldterminate = distributed_term_target[id]->curtoken.terminatenow;
    distributed_term_target[id]->curtoken = token;
    // first time I saw the termination. pass it along
    if (!oldterminate) {
    //  logger(LOG_INFO, "%d: Received Terminating Token", dc.procid());
      // let done() forward the token.
      distributed_term_target[id]->hastoken = true;
    }
    else {
      // loop around.
      distributed_term_target[id]->hastoken = true;
    }
  }
  else if (distributed_term_target[id]->curtoken.totaltransmitted == token.totaltransmitted &&
    distributed_term_target[id]->curtoken.totalreceived == token.totalreceived &&
    token.totaltransmitted == token.totalreceived) {
      // check if we may terminate
      // every looks right. we can terminate
     // logger(LOG_INFO, "%d: Creating Terminating Token", dc.procid());
      distributed_term_target[id]->curtoken = token;
      distributed_term_target[id]->curtoken.terminatenow = true;
      distributed_term_target[id]->hastoken = true;
      // let done() forward the token
  }
  else {
    // keep the token;
    distributed_term_target[id]->curtoken = token;
    distributed_term_target[id]->hastoken = true;
  }
}


bool distributed_terminator::done(size_t transmitted, size_t received) {
  if (quit) {
    return true;
  }
  my_total_lock.lock();

  delta_transmitted.inc(transmitted - my_total_transmitted);
  my_total_transmitted = transmitted;

  delta_received.inc(received - my_total_received);
  my_total_received = received;

  // now check if I have the token
  if (hastoken) {
    // we need to pass it on
    send_token();
  }
  my_total_lock.unlock();
  return quit;
}


void distributed_terminator::send_token() {
  ASSERT_TRUE(hastoken);

  size_t transmittedval = delta_transmitted.value;
  delta_transmitted.dec(transmittedval);
  curtoken.totaltransmitted += transmittedval;

  size_t receivedval = delta_received.value;
  delta_received.dec(receivedval);
  curtoken.totalreceived += receivedval;

  if (curtoken.terminatenow == false) {
   // logger(LOG_INFO, "%d: Pushing Token %d %d", dc.procid(),
    //                                            curtoken.totaltransmitted,
     //                                           curtoken.totalreceived);
  }
  else {
 //   logger(LOG_INFO, "%d: Pushing Terminating Token %d %d", dc.procid(),
  //                                              curtoken.totaltransmitted,
   //                                             curtoken.totalreceived);

    quit = true;
  }
  hastoken = false;
  if (use_control_packets == false) {
    dc.remote_callx((dc.procid() + 1) % dc.numprocs(),
                    distributed_terminator::distributed_terminator_handler,
                    NULL,
                    0,
                    termid,
                    curtoken);
  }
  else {
    dc.remote_call_control((dc.procid() + 1) % dc.numprocs(),
                        distributed_terminator::distributed_terminator_handler,
                        NULL,
                        0,
                        termid,
                        curtoken);
  }
}



}