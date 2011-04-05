

#include <graphlab/parallel/pthread_tools.hpp>
namespace graphlab { 

  // Some magic to ensure that keys are created at program startup =========>
  void destroy_tls_data(void* ptr);
  struct thread_keys {
    pthread_key_t GRAPHLAB_TSD_ID;
    thread_keys() : GRAPHLAB_TSD_ID(0) { 
      pthread_key_create(&GRAPHLAB_TSD_ID,
                         destroy_tls_data);
    }
  };
  static const thread_keys keys;
  // END MAGIC =============================================================>


  
  

  /**
   * Create thread specific data
   */
  thread::tls_data* create_tls_data(size_t thread_id = 0) {
    // Require that the data not yet exist
    assert(pthread_getspecific(keys.GRAPHLAB_TSD_ID) == NULL);
    // Create the data
    thread::tls_data* data =
      new thread::tls_data(thread_id);
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
  thread::tls_data& thread::get_tls_data() {
    // get the tsd
    tls_data* tsd =
      reinterpret_cast<tls_data*>
      (pthread_getspecific(keys.GRAPHLAB_TSD_ID));
    // If no tsd be has been associated, create one
    if(tsd == NULL) tsd = create_tls_data();
    assert(tsd != NULL);
    return *tsd;
  } // end of get thread specific data

  
  /**
   * Create thread specific data
   */
  void destroy_tls_data(void* ptr) {
    thread::tls_data* tsd =
      reinterpret_cast<thread::tls_data*>(ptr);
    if(tsd != NULL) {
      delete tsd;
    }
  } // end destroy the thread specific data

 


  //! Little helper function used to launch runnable objects      
  void* thread::invoke(void *_args) {
    invoke_args* args = static_cast<invoke_args*>(_args);
    // Create the graphlab thread specific data
    create_tls_data(args->m_thread_id);    
    //! Run the users thread code
    args->m_runnable->run();
    //! Delete the arguments 
    delete args;
    
    //! Properly kill the thread
    thread_destroy_callback();
    pthread_exit(NULL);
  } // end of invoke

  

  

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

