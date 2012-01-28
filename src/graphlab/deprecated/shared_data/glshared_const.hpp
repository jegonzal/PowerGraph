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


#ifndef GRAPHLAB_GLSHARED_CONST_HPP
#define GRAPHLAB_GLSHARED_CONST_HPP


#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/generics/any.hpp>



namespace graphlab {


  /**
   * \brief A constant shared data entry.  
   *
   * The glshared_const<datatype> container is used to store constants
   * that are to be accessed by update functions.  This interface
   * should be used instead of just creating regular global variables
   * as this behavior can be further optimized in distributed and NUMA 
   * settings
   *
   */
  template <typename T>
  class glshared_const  {

  public:

    typedef T value_type;
    typedef T& ref_type;

  private:

    //! The internal contents of the glshared_const object
    T content;

    bool finalized;

  public:
    glshared_const() : finalized(false) { }

    /**
     * Set the value of the constant.  This should only be called once
     * by a single thread.
     */
    void set(const T& c) { 
      ASSERT_FALSE(finalized);
      content = c; 
      finalized = true;
    }

    /**
     * Get the constant.
     */
    const T& get() const { 
      ASSERT_TRUE(finalized);
      return content; 
    }
    
    /**
     * Get the constant. Same as get().
     * Provided to get some semblance of interface similarity with glshared
     */
    const T& get_val() const { 
      ASSERT_TRUE(finalized);
      return content; 
    }
  };



}
#endif

