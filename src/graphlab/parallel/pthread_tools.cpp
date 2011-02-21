

#include <graphlab/parallel/pthread_tools.hpp>
namespace graphlab { 

  // Some magic to ensure that keys are created at program startup =========>
  void destroy_thread_specific_data(void* ptr);
  struct thread_keys {
    pthread_key_t GRAPHLAB_TSD_ID;
    thread_keys() : GRAPHLAB_TSD_ID(0) { 
      pthread_key_create(&GRAPHLAB_TSD_ID,
                         destroy_thread_specific_data);
    }
  };
  static const thread_keys keys;
  // END MAGIC =============================================================>


  
  

  /**
   * Create thread specific data
   */
  thread_specific_data* create_thread_specific_data(size_t thread_id = 0) {
    // Require that the data not yet exist
    assert(pthread_getspecific(keys.GRAPHLAB_TSD_ID) == NULL);
    // Create the data
    thread_specific_data* data =
      new thread_specific_data(thread_id);
    assert(data != NULL);
    // Set the data
    pthread_setspecific(keys.GRAPHLAB_TSD_ID, data);
    // Return the associated tsd
    return data;
  } // end create the thread specific data

  /**
   * This function tries to get the thread specific data.  If no
   * thread specific data has been associated with the thread than it
   * is created.
   */
  thread_specific_data& thread::get_thread_specific_data() {
    // get the tsd
    thread_specific_data* tsd =
      reinterpret_cast<thread_specific_data*>
      (pthread_getspecific(keys.GRAPHLAB_TSD_ID));
    // If no tsd be has been associated, create one
    if(tsd == NULL) tsd = create_thread_specific_data();
    assert(tsd != NULL);
    return *tsd;
  } // end of get thread specific data

  
  /**
   * Create thread specific data
   */
  void destroy_thread_specific_data(void* ptr) {
    thread_specific_data* tsd =
      reinterpret_cast<thread_specific_data*>(ptr);
    if(tsd != NULL) {
      delete tsd;
    }
  } // end destroy the thread specific data

 


  //! Little helper function used to launch runnable objects      
  void* thread::invoke(void *_args) {
    invoke_args* args = static_cast<invoke_args*>(_args);
    // Create the graphlab thread specific data
    create_thread_specific_data(args->m_thread_id);    
    //! Run the users thread code
    args->m_runnable->run();
    //! Delete the arguments 
    delete args;
    
    //! Properly kill the thread
    thread_destroy_callback();
    pthread_exit(NULL);
  } // end of invoke

  
  /** Get the id of the calling thread.  This will typically be the
      index in the thread group. Between 0 to ncpus. */
  size_t thread::thread_id() {
    return get_thread_specific_data().thread_id();
  }

  

  /**
   * This static method joins the invoking thread with the other
   * thread object.  This thread will not return from the join
   * routine until the other thread complets it run.
   */
  void thread::join(thread& other) {
    void *status;
    // joint the first element
    int error = 0;
    if(other.active())
      error = pthread_join( other.m_p_thread, &status);
    if(error) {
      std::cout << "Major error in join" << std::endl;
      std::cout << "pthread_join() returned error " << error << std::endl;
      exit(EXIT_FAILURE);
    }
  } // end of join


  /**
   * Return the number processing units (individual cores) on this
   * system
   */
  size_t thread::cpu_count() {
#if defined __linux__
    return sysconf(_SC_NPROCESSORS_CONF);
#elif defined(__MACH__) && defined(_SC_NPROCESSORS_ONLN)
    return sysconf (_SC_NPROCESSORS_ONLN);
#elif defined(__MACH__) && defined(HW_NCPU)
    int ncpus = 1;
    size_t len = sizeof(ncpus);
    sysctl((int[2]) {CTL_HW, HW_NCPU}, 2, &ncpus, &len, NULL, 0);
    return ncpus;
#else
    return 0;
#endif
  } // end of cpu count
    
   /**
     * Allow defining a callback when thread is destroyed.
     * This is needed at least from Java JNI, where we have to detach
     * thread from JVM before it dies.
     */
   void (*__thr_callback)()  = NULL;

   void thread::thread_destroy_callback() {
     if (__thr_callback != NULL) __thr_callback();
   }
   
   void thread::set_thread_destroy_callback(void (*callback)()) {
     __thr_callback = callback;
   }

  
  
}




// void create_generators(int seed_value) {
//   //! create a random number generator for the thread
//   if(pthread_getspecific(keys.RAND_SRC_ID) == NULL) {
//     thread::rand_src_type* rand_gen = 
//       new thread::rand_src_type((int)seed_value);
//     // rand_gen->engine().seed((int)seed_value);
//     pthread_setspecific(keys.RAND_SRC_ID, rand_gen);
//   }
//   if(pthread_getspecific(keys.RAND01_ID) == NULL) {
//     thread::rand_src_type* rand_gen = 
//       reinterpret_cast<thread::rand_src_type*>
//       (pthread_getspecific(keys.RAND_SRC_ID));
//     assert(rand_gen != NULL);
//     rand_01_type* rand01 = new rand_01_type(*rand_gen, dist_01_type(0,1));
//     pthread_setspecific(keys.RAND01_ID, rand01);
//   }
// }

// void create_generators() {
//   //! create a random number generator for the thread
//   if(pthread_getspecific(keys.RAND_SRC_ID) == NULL) {
//     thread::rand_src_type* rand_gen = new thread::rand_src_type();
//     pthread_setspecific(keys.RAND_SRC_ID, rand_gen);
//   }
//   if(pthread_getspecific(keys.RAND01_ID) == NULL) {
//     thread::rand_src_type* rand_gen = 
//       reinterpret_cast<thread::rand_src_type*>
//       (pthread_getspecific(keys.RAND_SRC_ID));
//     assert(rand_gen != NULL);
//     rand_01_type* rand01 = new rand_01_type(*rand_gen, dist_01_type(0,1));
//     pthread_setspecific(keys.RAND01_ID, rand01);
//   }
// }

