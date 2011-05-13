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

