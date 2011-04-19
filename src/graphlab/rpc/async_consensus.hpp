#ifndef ASYNC_TERMINATOR_HPP
#define ASYNC_TERMINATOR_HPP

#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object_base.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>


namespace graphlab {
  /**
   * This implements a distributed consensus algorithm attached on an
   * object
   */
  class async_consensus {
  public:
    /// attaches to the "attach" object if provided. if NULL, attaches to the main DC
    async_consensus(distributed_control &dc, size_t required_threads_in_done = 1,
                    const dc_impl::dc_dist_object_base* attach = NULL);

    /**
     * Done blocks. It only returns if a cancellation
     * or a completion occurs.
     * Done may be called by many threads.
     * Returns true if complete.
     * Returns false if cancelled.
     */
    bool done();

    void begin_done_critical_section();

    bool end_done_critical_section(bool done);

    void force_done();
  
    bool done_noblock() {
      return complete;
    }

    /**
     * Cancels a "done" call. done() will immediately return
     * false after this. Note that this function is rather costly
     * and could "serialize" your execution so it should not be
     * called too frequently
     */
    void cancel();

    void cancel_one();

    struct token {
      size_t total_calls_sent;
      size_t total_calls_received;
      procid_t last_change;
      void save(oarchive &oarc) const {
        oarc << total_calls_sent << total_calls_received << last_change;
      }

      void load(iarchive &iarc) {
        iarc >> total_calls_sent >> total_calls_received >> last_change;
      }
    };
  
  private:
    dc_dist_object<async_consensus> rmi;
    const dc_impl::dc_dist_object_base* attachedobj;
  
    size_t last_calls_sent;
    size_t last_calls_received;

    size_t required_threads_in_done;
  
    size_t threads_in_done;
    /// set if a cancellation occurs while a thread is waiting in done()
    size_t cancelled;
    /// set when everyone is done
    bool complete;
    /// whether I currently have the token
    bool hastoken;
    /// If I have the token, the value of the token
    token cur_token;

    mutex mut;
    conditional cond;
  

    void receive_the_token(token &tok);
    void pass_the_token();
    void consensus();
  };

}
#endif



