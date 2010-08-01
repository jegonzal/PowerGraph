#ifndef DISTRIBUTED_TERMINATOR_HPP
#define DISTRIBUTED_TERMINATOR_HPP
#include <vector>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/distributed/distributed_control_types.hpp>

namespace graphlab {

  class distributed_control;
  struct distributed_termination_token {
    size_t totaltransmitted;
    size_t totalreceived;
    bool terminatenow;
  };

  class distributed_terminator{
   private:
    distributed_control &dc;
    distributed_termination_token curtoken;

    // the total I have contributed to the token thus far
    // this is protected by the mutex
    mutex my_total_lock;
    size_t my_total_transmitted;
    size_t my_total_received;

    // increments to add to the token next time I send it
    // protected by atomics
    atomic<size_t> delta_transmitted;
    atomic<size_t> delta_received;

    bool hastoken;
    bool quit;

    bool use_control_packets;

    size_t termid;
    static std::vector<distributed_terminator*> distributed_term_target;
    
   public:
    distributed_terminator(distributed_control &dc);

    inline bool terminatenow() {
      return quit;
    }

    inline void reset() {
      quit = false;
      curtoken.terminatenow = false;
    }

    void set_use_control_packets(bool ctrl) {
      use_control_packets = ctrl;
    }
    // if done returns true, distributed termination is complete
    // even when locally done, this function still needs to be called
    // ocassionally. (possibly with sched_yields outside)
    bool done(size_t transmitted, size_t received);

    // message handler
    static void distributed_terminator_handler(distributed_control &dc,
                                        procid_t source,
                                        void * ptr,    //unused
                                        size_t len,   // unused
                                        size_t id,
                                        distributed_termination_token token);

   private:
    // send the token. Caller must be sure that I have the token
    void send_token();

  };

}
#endif
