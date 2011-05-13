/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <graphlab/parallel/pthread_tools.hpp>
#include <boost/bind.hpp>
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
    void* retval = NULL;
    thread::invoke_args* args = static_cast<thread::invoke_args*>(_args);
    // Create the graphlab thread specific data
    create_tls_data(args->m_thread_id);    
    //! Run the users thread code
    try {
      args->spawn_routine();
    }
    catch (const char* msg) {
      retval = (void*)msg;
    }
    //! Delete the arguments 
    delete args;
    
    //! Properly kill the thread
    thread_destroy_callback();
    return retval;
  } // end of invoke

  

  

  /**
   * This static method joins the invoking thread with the other
   * thread object.  This thread will not return from the join
   * routine until the other thread complets it run.
   */
  void thread::join(thread& other) {
    void *status = NULL;
    // joint the first element
    int error = 0;
    if(other.active()) {
      error = pthread_join( other.m_p_thread, &status);
      if (status != NULL) {
        const char* strstatus = (const char*) status;
        throw strstatus;
      }
    }
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
      if (cpu_id  == size_t(-1)) {
        launch(spawn_routine);
        return;
      }
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

#ifdef HAS_SET_AFFINITY
      // Set Processor Affinity masks (linux only)
      cpu_set_t cpu_set;
      CPU_ZERO(&cpu_set);
      CPU_SET(cpu_id % CPU_SETSIZE, &cpu_set);

      pthread_attr_setaffinity_np(&attr, sizeof(cpu_set), &cpu_set);
#endif
          
      // Launch the thread
      error = pthread_create(&m_p_thread, 
                             &attr, 
                             invoke,
                             static_cast<void*>(new invoke_args(m_thread_id, spawn_routine)));
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
  // thread group exception forwarding is a little more complicated
  // because it has to be able to catch it on a bunch of threads

  void thread_group::invoke(boost::function<void (void)> spawn_function,
                        thread_group *group) {
    const char* retval = NULL;
    try {
      spawn_function();
    }
    catch (const char* c) {
      // signal the thread group to join this thread
      retval = c;
    }
      group->mut.lock();
      group->joinqueue.push(std::make_pair(pthread_self(), retval));
      group->cond.signal();
      group->mut.unlock();

  }
                        

  void thread_group::launch(const boost::function<void (void)> &spawn_function) {
    // Create a thread object and launch it. 
    // We do not need to keep a copy of the thread around
    thread local_thread(m_thread_counter++);
    mut.lock();
    threads_running++;
    mut.unlock();
    local_thread.launch(boost::bind(thread_group::invoke, spawn_function, this));
  } 


  void thread_group::launch(const boost::function<void (void)> &spawn_function, 
                            size_t cpu_id) {
    if (cpu_id == size_t(-1)) {
      launch(spawn_function);
      return;
    }
    // Create a thread object
    thread local_thread(m_thread_counter++);
    mut.lock();
    threads_running++;
    mut.unlock();
    local_thread.launch(boost::bind(thread_group::invoke, spawn_function, this), cpu_id);
  }

  void thread_group::join() {
    mut.lock();
    while(threads_running > 0) {
      // if no threads are joining. wait
      while (joinqueue.empty()) cond.wait(mut);      
      // a thread is joining
      std::pair<pthread_t, const char*> joining_thread = joinqueue.front();
      joinqueue.pop();
      threads_running--;
      
      // unlock here since I might be in join for a little while
      mut.unlock();
      void *unusedstatus = NULL;
      pthread_join(joining_thread.first, &unusedstatus);
      // if there is a return value
      // throw it. It is safe to throw here since I have the mutex unlocked.
      if (joining_thread.second) {
        throw(joining_thread.second);
      }
      mut.lock();
    }
    mut.unlock();
  }


}

