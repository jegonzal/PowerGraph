#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/macros_def.hpp>

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

// -----------------------------------------------------------------
//                 Thread Object Static Members 
// -----------------------------------------------------------------
  

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

 


  //! Little helper function used to launch threads
  void* thread::invoke(void *_args) {
    thread::invoke_args* args = static_cast<thread::invoke_args*>(_args);
    // Create the graphlab thread specific data
    create_tls_data(args->m_thread_id);    
    //! Run the users thread code
    args->spawn_routine();
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


// -----------------------------------------------------------------
//                 Thread Object Public Members 
// -----------------------------------------------------------------

  
  void thread::launch(const boost::function<void (void)> &spawn_routine) {
    ASSERT_FALSE(thread_started);
    // fill in the thread attributes
    pthread_attr_t attr;
    int error = 0;
    error = pthread_attr_init(&attr);
    ASSERT_TRUE(!error);
    error = pthread_attr_setstacksize(&attr, m_stack_size);
    ASSERT_TRUE(!error);
    error = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    ASSERT_TRUE(!error);
    


    error = pthread_create(&m_p_thread, 
                            &attr, 
                            invoke,  
                            static_cast<void*>(new invoke_args(m_thread_id, spawn_routine)) );
    thread_started = true;
    if(error) {
      std::cout << "Major error in thread_group.launch (pthread_create). Error: " << error << std::endl;
      exit(EXIT_FAILURE);
    }
    // destroy the attribute object
    error = pthread_attr_destroy(&attr);
    ASSERT_TRUE(!error);
  }
  
  void thread::launch(const boost::function<void (void)> &spawn_routine, size_t cpu_id){
      // if this is not a linux based system simply invoke start and
      // return;
#ifndef __linux__
      launch(spawn_routine);
      return;
#else
      ASSERT_FALSE(thread_started);
      if (cpu_id >= cpu_count() && cpu_count() > 0) {
        // definitely invalid cpu_id
        std::cout << "Invalid cpu id passed on thread_ground.launch()" 
                  << std::endl;
        std::cout << "CPU " << cpu_id + 1 << " requested, but only " 
                  << cpu_count() << " CPUs available" << std::endl;
        exit(EXIT_FAILURE);
      }
      
      // fill in the thread attributes
      pthread_attr_t attr;
      int error = 0;
      error = pthread_attr_init(&attr);
      ASSERT_TRUE(!error);
      error = pthread_attr_setstacksize(&attr, m_stack_size);
      ASSERT_TRUE(!error);
      error = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      ASSERT_TRUE(!error);
      
      // Set Processor Affinity masks (linux only)
      cpu_set_t cpu_set;
      CPU_ZERO(&cpu_set);
      CPU_SET(cpu_id % CPU_SETSIZE, &cpu_set);
      pthread_attr_setaffinity_np(&attr, sizeof(cpu_set), &cpu_set);

      // Launch the thread
      error = pthread_create(&m_p_thread, 
                             &attr, 
                             invoke,
                             static_cast<void*>(new invoke_args(m_thread_id, spawn_routine)) );
      thread_started = true;
      if(error) {
        std::cout << "Major error in thread_group.launch" << std::endl;
        std::cout << "pthread_create() returned error " << error << std::endl;
        exit(EXIT_FAILURE);
      }
      
      
      
      // destroy the attribute object
      error = pthread_attr_destroy(&attr);
      ASSERT_TRUE(!error);
#endif
    }
    
// -----------------------------------------------------------------
//                 Thread Group Object Public Members 
// -----------------------------------------------------------------

  void thread_group::launch(const boost::function<void (void)> &spawn_function) {
    // Create a thread object
    thread local_thread(m_thread_counter++);
    local_thread.launch(spawn_function);
    // keep a local copy of the thread
    m_threads.push_back(local_thread);
  } 


  void thread_group::launch(const boost::function<void (void)> &spawn_function, 
                            size_t cpu_id) {
    // Create a thread object
    thread local_thread(m_thread_counter++);
    local_thread.launch(spawn_function, cpu_id);
    // keep a local copy of the thread
    m_threads.push_back(local_thread);
  }
  
  void thread_group::join() {
    while(!m_threads.empty()) {
      m_threads.front().join(); // Join the first thread
      m_threads.pop_front(); // remove the first element
    }
  }

  void thread_group::signalall(int sig) {
    foreach (thread& t, m_threads) {
      pthread_kill(t.pthreadid(), sig);
    }
  }
}

