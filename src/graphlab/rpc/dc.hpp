#ifndef GRAPHLAB_DC_HPP
#define GRAPHLAB_DC_HPP
#include <iostream>
#include <boost/iostreams/stream.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/resizing_array_sink.hpp>
#include <graphlab/util/blocking_queue.hpp>
#include <graphlab/util/multi_blocking_queue.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/metrics/metrics.hpp>

#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>

#include <graphlab/rpc/dc_receive.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/rpc/dc_comm_base.hpp>
#include <graphlab/rpc/dc_dist_object_base.hpp>

#include <graphlab/rpc/is_rpc_call.hpp>
#include <graphlab/rpc/portable_dispatch.hpp>
#include <graphlab/rpc/portable_issue.hpp>
#include <graphlab/rpc/function_call_issue.hpp>
#include <graphlab/rpc/request_issue.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/function_ret_type.hpp>

#include <boost/preprocessor.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>

namespace graphlab {


/**
Struct form of the distributed control constructor parameters.
This will allow curious initialization functions to be written.
For instance:
distributed_control dc(dc_MPI_init(argc, argv));
*/
struct dc_init_param{
  std::vector<std::string> machines;
  std::string initstring;
  procid_t curmachineid;
  size_t numhandlerthreads;
  dc_comm_type commtype;
};

#define DEFAULT_NUMHANDLERTHREADS 8
#define DEFAULT_COMMTYPE TCP_COMM

// forward declaration for dc services
class dc_services;
class distributed_control;


/**
The primary distributed RPC object.
The basic operation goes like this:
the distributed_control will initialize the communication type using the template
argument provided. It will then spawn a receive multiplxer and a send multiplexer.
The receive multiplexer is attached the communication object using the call back
interface, while the send multiplexer is provided directly access to the comm object.

*/
class distributed_control{
  public:
        /**  Each element of the function call queue is a data/len pair */
    struct function_call_block{
      function_call_block() {}
      function_call_block(procid_t source, const dc_impl::packet_hdr& hdr, 
                          char* data, size_t len): 
                          source(source), hdr(hdr), 
                          data(data), len(len) {}
      procid_t source;
      dc_impl::packet_hdr hdr;
      char* data;
      size_t len;
    };
  private:
   /// initialize receiver threads. private form of the constructor
   void init(const std::vector<std::string> &machines,
             const std::string &initstring,
             procid_t curmachineid,
             size_t numhandlerthreads,
             dc_comm_type commtype = DEFAULT_COMMTYPE);
   
  /// a pointer to the communications subsystem
  dc_impl::dc_comm_base* comm; 
 
  /// senders and receivers to all machines
  std::vector<dc_impl::dc_receive*> receivers;
  std::vector<dc_impl::dc_send*> senders;
  
  /// A thread group of function call handlers
  thread_group fcallhandlers;
  
  /// a queue of functions to be executed
  multi_blocking_queue<function_call_block> fcallqueue;
  
  /// A map of function name to dispatch function. Used for "portable" calls
  dc_impl::dispatch_map_type portable_dispatch_call_map;
  dc_impl::dispatch_map_type portable_dispatch_request_map;

  
  /// object registrations;
  std::vector<void*> registered_objects;
  std::vector<dc_impl::dc_dist_object_base*> registered_rmi_instance;

  /// For convenience, we provide a instance of dc_services 
  dc_services* distributed_services;

  /// ID of the local machine
  procid_t localprocid;

  
  /// Number of machines
  procid_t localnumprocs;
  
  std::vector<atomic<size_t> > global_calls_sent;
  std::vector<atomic<size_t> > global_calls_received;
  
  bool single_sender;
  
  /// the callback given to the comms class. Called when data is inbound
  friend void dc_recv_callback(void* tag, procid_t src, const char* buf, size_t len);
  
  
  /// the callback given to the comms class. Called when data is inbound
  template <typename T> friend class dc_dist_object;
  
  
  /// disable the operator= by placing it in private 
  distributed_control& operator=(const distributed_control& dc) { return *this; }


  std::map<std::string, std::string> parse_options(std::string initstring);
  
  volatile inline size_t num_registered_objects() {
    return registered_objects.size();
  }
  

  void compute_master_ranks();
  procid_t masterid;
  
  metrics rpc_metrics;
  
 public:
   

  
  distributed_control(dc_init_param initparam) {
    init(initparam.machines, 
         initparam.initstring, 
         initparam.curmachineid, 
         initparam.numhandlerthreads,
         initparam.commtype);
  }

  distributed_control(const std::vector<std::string> &machines,
                      const std::string &initstring, 
                      procid_t curmachineid, 
                      size_t numhandlerthreads = DEFAULT_NUMHANDLERTHREADS,
                      dc_comm_type commtype = DEFAULT_COMMTYPE) {
    init(machines, initstring, curmachineid, numhandlerthreads, commtype);
  }

  ~distributed_control();

  /// returns the id of the current processor
  inline procid_t procid() const {
    return localprocid;
  }
  
  /// returns the number of processors in total.
  inline procid_t numprocs() const {
    return localnumprocs;
  }
  
  /**
  Sets the sequentialization key to the new value, returning the previous value.
  User should 
  oldval = set_sequentialization_key(newval)
  ...
  ... do stuff
  ...
  set_sequentialization_key(oldval)
  */
  static unsigned char set_sequentialization_key(unsigned char newkey);
  
  /**
  Creates a new sequentialization key, returning the old value
  User should 
  oldval = new_sequentialization_key()
  ...
  ... do stuff
  ...
  set_sequentialization_key(oldval)
  */
  static unsigned char new_sequentialization_key();
  static unsigned char get_sequentialization_key();
  
  /**
  This generates the interface functions for the standard calls, basic calls, and fast calls
  The generated code looks like this:
  
  template<typename F , typename T0> void remote_call (procid_t target, F remote_function , const T0 &i0 )
  {
    ASSERT_LT(target, senders.size());
    dc_impl::remote_call_issue1 <F , T0> ::exec(senders[target], 
                                                STANDARD_CALL, 
                                                target, 
                                                remote_function , 
                                                i0 );                                
  }
  The arguments passed to the RPC_INTERFACE_GENERATOR ARE: (interface name, issue processor name, flags)

  */
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
  #define GENI(Z,N,_) BOOST_PP_CAT(i, N)
  #define GENT(Z,N,_) BOOST_PP_CAT(T, N)
  #define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);

  #define RPC_INTERFACE_GENERATOR(Z,N,FNAME_AND_CALL) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  void  BOOST_PP_TUPLE_ELEM(3,0,FNAME_AND_CALL) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, senders.size()); \
    BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,FNAME_AND_CALL),N) \
        <F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(senders[target],  BOOST_PP_TUPLE_ELEM(3,2,FNAME_AND_CALL), target, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \
  
  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (remote_call, dc_impl::remote_call_issue, STANDARD_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (fast_remote_call,dc_impl::remote_call_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, RPC_INTERFACE_GENERATOR, (control_call, dc_impl::remote_call_issue, (FAST_CALL | CONTROL_PACKET)) )
 

  #define REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
  template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
    BOOST_PP_TUPLE_ELEM(3,0,ARGS) (procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    ASSERT_LT(target, senders.size()); \
    return BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(3,1,ARGS),N) \
        <F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> \
          ::exec(senders[target],  BOOST_PP_TUPLE_ELEM(3,2,ARGS), target, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  }   \

  /**
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
 // BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type remote_request, dc_impl::remote_request_issue, 0) )
   BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type remote_request, dc_impl::remote_request_issue, STANDARD_CALL) )

  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type fast_remote_request, dc_impl::remote_request_issue, FAST_CALL) )
  BOOST_PP_REPEAT(6, REQUEST_INTERFACE_GENERATOR, (typename dc_impl::function_ret_type<FRESULT>::type control_request, dc_impl::remote_request_issue, (FAST_CALL | CONTROL_PACKET)) )
 

  
  #undef RPC_INTERFACE_GENERATOR
  #undef REQUEST_INTERFACE_GENERATOR
  #undef GENARC
  #undef GENT
  #undef GENI
  #undef GENARGS
  
  /**
  Immediately calls the function described by the data
  inside the buffer. This should not be called directly.
  */
  void exec_function_call(procid_t source, const dc_impl::packet_hdr& hdr, std::istream &istrm);
  
  
  
  /**
  Performs a deferred function call using the information
  inside the buffer. This function will take over ownership of 
  the buffer and will free it when done
  */
  void deferred_function_call(procid_t source, const dc_impl::packet_hdr& hdr, 
                              char* buf, size_t len);
  

  /**
  This is called by the function handler threads
  */
  void fcallhandler_loop(size_t id);
  
  inline void inc_calls_sent(procid_t procid) {
//    ASSERT_FALSE(full_barrier_in_effect);
    global_calls_sent[procid].inc();
  }

  inline void inc_calls_received(procid_t procid) {
    
    if (!full_barrier_in_effect) {
        global_calls_received[procid].inc();
    }
    else {
      //check the proc I just incremented.
      // If I just exceeded the required size, I need
      // to decrement the full barrier counter
      if (global_calls_received[procid].inc() == calls_to_receive[procid]) {
        // if it was me who set the bit
        if (procs_complete.set_bit(procid) == false) {
          // then decrement the incomplete count.
          // if it was me to decreased it to 0
          // lock and signal
          full_barrier_lock.lock();
          if (num_proc_recvs_incomplete.dec() == 0) {
            full_barrier_cond.signal();
          }
          full_barrier_lock.unlock();
        }
      }
    }

  }

  inline size_t calls_sent() const {
    size_t ctr = 0;
    for (size_t i = 0;i < numprocs(); ++i) {
      ctr += global_calls_sent[i].value;
    }
    return ctr;
  }

  inline size_t calls_received() const {
    size_t ctr = 0;
    for (size_t i = 0;i < numprocs(); ++i) {
      ctr += global_calls_received[i].value;
    }
    return ctr;
  }

  inline size_t bytes_sent() const {
    if (single_sender) {
      return senders[0]->bytes_sent();
    }
    else {
      size_t ret = 0;
      for (size_t i = 0;i < senders.size(); ++i) ret += senders[i]->bytes_sent();
      return ret;
    }
  }  
  
  
  inline size_t network_bytes_sent() const {
    return comm->network_bytes_sent();
  }  
  
  inline size_t bytes_received() const {
    size_t ret = 0;
    for (size_t i = 0;i < receivers.size(); ++i) ret += receivers[i]->bytes_received();
    return ret;
  }  

  /**
    Returns true if this is the process with the lowest ID
    currently running on this machine in this working directory
  */
  inline bool is_master_rank() const {
    return masterid == procid();
  }


  /**
    Returns the lowest ID of all the processes
    currently running on this machine in this working directory
  */  
  inline procid_t master_rank() const {
    return masterid;
  }

  /**
    Instantiates a find_dispatch with the right arguments,
    and store the dispatch function in the hash map.
  */
  template <typename F, F f>
  void register_rpc(std::string c) {
    portable_dispatch_request_map[c] = (dc_impl::dispatch_type)
              dc_impl::portable_detail::find_dispatcher<F,        // function type
                              FRESULT,                            // result
                              boost::function_traits<               
                                    typename boost::remove_pointer<F>::type
                                                    >::arity ,   // number of arguments
                              f,                                    // function itself
                              typename dc_impl::is_rpc_call<F>::type  // whether it is an RPC style call
                              >::dispatch_request_fn();
                              
    portable_dispatch_call_map[c] = (dc_impl::dispatch_type)
              dc_impl::portable_detail::find_dispatcher<F,        // function type
                              FRESULT,                            // result
                              boost::function_traits<               
                                    typename boost::remove_pointer<F>::type
                                                    >::arity ,   // number of arguments
                              f,                                    // function itself
                              typename dc_impl::is_rpc_call<F>::type  // whether it is an RPC style call
                              >::dispatch_call_fn();
  }


  inline size_t register_object(void* v, dc_impl::dc_dist_object_base *rmiinstance) {
    ASSERT_NE(v, (void*)NULL);
    registered_objects.push_back(v);
    registered_rmi_instance.push_back(rmiinstance);
    return registered_objects.size() - 1;
  }


  inline void* get_registered_object(size_t id) {
    while(id >= num_registered_objects()) sched_yield();
    ASSERT_NE(registered_objects[id], (void*)NULL);
    return registered_objects[id];
  }

  inline dc_impl::dc_dist_object_base* get_rmi_instance(size_t id) {
    while(id >= num_registered_objects()) sched_yield();
    ASSERT_NE(registered_rmi_instance[id], (void*)NULL);
    return registered_rmi_instance[id];
  }  
  inline void clear_registered_object(size_t id) {
    registered_objects[id] = (void*)NULL;
    registered_rmi_instance[id] = NULL;
  }
  
  
  /**
  This is depreated. Use the public functions. In particular
  services().full_barrier() may not work as expected
  */
  __attribute__((__deprecated__)) dc_services& services();
  
  /**
   This comm barrier is not a true "barrier" but is
   essentially a sequentialization point. It guarantees that
   all calls from this machine to the target machine performed
   before the comm_barrier() call are completed before any call
   sent after the comm barrier() call.
  */
  void comm_barrier(procid_t targetmachine);
  
  /**
    This is a convenience function which broadcasts a comm_barrier()
    \note having all machines call the comm barrier does not guarantee
    that all calls have been processed. Basically 'p' local barriers
    do not result in a global barrier.
  */
  void comm_barrier();



// Temp hack.
  long long int total_bytes_sent;
  long long int get_total_bytes_sent() {
     return total_bytes_sent;
  }

  template <typename U>
  inline void send_to(procid_t target, U& t, bool control = false);
  
  template <typename U>
  inline void recv_from(procid_t source, U& t, bool control = false);

  template <typename U>
  inline void broadcast(U& data, bool originator, bool control = false);

  template <typename U>
  inline void gather(std::vector<U>& data, procid_t sendto, bool control = false);

  template <typename U>
  inline void all_gather(std::vector<U>& data, bool control = false);

  template <typename U>
  inline void gather_partition(const std::vector<U>& local_contribution,
                        std::vector< std::vector<U> >& ret_partition,
                        bool control = false);
  

  void barrier();
  



 /*****************************************************************************
                      Implementation of Full Barrier
*****************************************************************************/
  void full_barrier();
 private:
  mutex full_barrier_lock;
  conditional full_barrier_cond;
  std::vector<size_t> calls_to_receive;
  // used to inform the counter that the full barrier
  // is in effect and all modifications to the calls_recv
  // counter will need to lock and signal
  bool full_barrier_in_effect;
  
  /** number of 'source' processor counts which have
  not achieved the right recv count */
  atomic<size_t> num_proc_recvs_incomplete; 
                                      
  /// Marked as 1 if the proc is complete
  dense_bitset procs_complete;
  
 /*****************************************************************************
                      Collection of Statistics
*****************************************************************************/
 
 private:
  struct collected_statistics {
    size_t callssent;
    size_t bytessent;
    size_t network_bytessent;
    collected_statistics(): callssent(0), bytessent(0), network_bytessent(0) { }
    void save(oarchive &oarc) const {
      oarc << callssent << bytessent << network_bytessent;
    }
    void load(iarchive &iarc) {
      iarc >> callssent >> bytessent >> network_bytessent;
    }
  };
 public:
  /** Gather RPC statistics. All machines must call 
   this function at the same time. However, only proc 0 will
   return values */
  std::map<std::string, size_t> gather_statistics();

  /** Fills metrics information. All machines must call this
   * function simultaneously. Only proc 0 will have metrics
   */
  void fill_metrics();

  /** returns metrics information collected by fill_metrics
   *  Only proc 0 will have metrics
   */
  inline metrics get_metrics() {
    return rpc_metrics;
  }

  inline void reset_metrics() {
    logstream(LOG_WARNING) << "Metrics cannot be reset on distributed control" << std::endl;
  }
  
  /** Dumps the metric information to a reporter
   * Only proc 0 will have metrics
   */
  inline void report_metrics(imetrics_reporter &reporter) {
    rpc_metrics.report(reporter);
  }

};




} // namespace graphlab

#define REGISTER_RPC(dc, f) dc.register_rpc<typeof(f)*, f>(std::string(BOOST_PP_STRINGIZE(f))) 

#include <graphlab/rpc/function_arg_types_undef.hpp>
#include <graphlab/rpc/function_call_dispatch.hpp>
#include <graphlab/rpc/request_dispatch.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/dc_services.hpp>

namespace graphlab {

template <typename U>
inline void distributed_control::send_to(procid_t target, U& t, bool control) {
  distributed_services->send_to(target, t, control);
}

template <typename U>
inline void distributed_control::recv_from(procid_t source, U& t, bool control) {
  distributed_services->recv_from(source, t, control);
}

template <typename U>
inline void distributed_control::broadcast(U& data, bool originator, bool control) { 
  distributed_services->broadcast(data, originator, control);
}

template <typename U>
inline void distributed_control::gather(std::vector<U>& data, procid_t sendto, bool control) {
  distributed_services->gather(data, sendto, control);
}

template <typename U>
inline void distributed_control::all_gather(std::vector<U>& data, bool control) {
  distributed_services->all_gather(data, control);
}

template <typename U>
inline void distributed_control::gather_partition(const std::vector<U>& local_contribution,
                      std::vector< std::vector<U> >& ret_partition,
                      bool control) {
  distributed_services->gather_partition(local_contribution, ret_partition, control);
}


}
#endif
