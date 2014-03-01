/**
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

#ifndef GRAPHLAB_MUTEX_HPP
#define GRAPHLAB_MUTEX_HPP


#include <pthread.h>
#include <graphlab/logger/assertions.hpp>


namespace graphlab {

  /**
   * \ingroup util
   *
   * Simple wrapper around pthread's mutex.
   * Before you use, see \ref parallel_object_intricacies.
   */
  class mutex {
  public:
    // mutable not actually needed
    mutable pthread_mutex_t m_mut;
    /// constructs a mutex
    mutex() {
      int error = pthread_mutex_init(&m_mut, NULL);
      ASSERT_TRUE(!error);
    }
    /** Copy constructor which does not copy. Do not use!
        Required for compatibility with some STL implementations (LLVM).
        which use the copy constructor for vector resize,
        rather than the standard constructor.    */
    mutex(const mutex&) {
      int error = pthread_mutex_init(&m_mut, NULL);
      ASSERT_TRUE(!error);
    }

    ~mutex(){
      int error = pthread_mutex_destroy( &m_mut );
      ASSERT_TRUE(!error);
    }

    // not copyable
    void operator=(const mutex& m) { }

    /// Acquires a lock on the mutex
    inline void lock() const {
      int error = pthread_mutex_lock( &m_mut  );
      // if (error) std::cout << "mutex.lock() error: " << error << std::endl;
      ASSERT_TRUE(!error);
    }
    /// Releases a lock on the mutex
    inline void unlock() const {
      int error = pthread_mutex_unlock( &m_mut );
      ASSERT_TRUE(!error);
    }
    /// Non-blocking attempt to acquire a lock on the mutex
    inline bool try_lock() const {
      return pthread_mutex_trylock( &m_mut ) == 0;
    }
    friend class conditional;
  }; // End of Mutex




  /**
   * \ingroup util
   *
   * Simple wrapper around pthread's recursive mutex.
   * Before you use, see \ref parallel_object_intricacies.
   */
  class recursive_mutex {
  public:
    // mutable not actually needed
    mutable pthread_mutex_t m_mut;
    /// constructs a mutex
    recursive_mutex() {
      pthread_mutexattr_t attr;
      int error = pthread_mutexattr_init(&attr);
      ASSERT_TRUE(!error);
      error = pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
      ASSERT_TRUE(!error);
      error = pthread_mutex_init(&m_mut, &attr);
      ASSERT_TRUE(!error);
      pthread_mutexattr_destroy(&attr);
    }
    /** Copy constructor which does not copy. Do not use!
        Required for compatibility with some STL implementations (LLVM).
        which use the copy constructor for vector resize,
        rather than the standard constructor.    */
    recursive_mutex(const recursive_mutex&) {
      pthread_mutexattr_t attr;
      int error = pthread_mutexattr_init(&attr);
      ASSERT_TRUE(!error);
      error = pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
      ASSERT_TRUE(!error);
      error = pthread_mutex_init(&m_mut, &attr);
      ASSERT_TRUE(!error);
      pthread_mutexattr_destroy(&attr);
    }

    ~recursive_mutex(){
      int error = pthread_mutex_destroy( &m_mut );
      ASSERT_TRUE(!error);
    }

    // not copyable
    void operator=(const recursive_mutex& m) { }

    /// Acquires a lock on the mutex
    inline void lock() const {
      int error = pthread_mutex_lock( &m_mut  );
      // if (error) std::cout << "mutex.lock() error: " << error << std::endl;
      ASSERT_TRUE(!error);
    }
    /// Releases a lock on the mutex
    inline void unlock() const {
      int error = pthread_mutex_unlock( &m_mut );
      ASSERT_TRUE(!error);
    }
    /// Non-blocking attempt to acquire a lock on the mutex
    inline bool try_lock() const {
      return pthread_mutex_trylock( &m_mut ) == 0;
    }
    friend class conditional;
  }; // End of Mutex




} // end of graphlab namespace


#endif
