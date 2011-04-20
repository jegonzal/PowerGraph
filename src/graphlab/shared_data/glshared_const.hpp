#ifndef GRAPHLAB_GLSHARED_CONST_HPP
#define GRAPHLAB_GLSHARED_CONST_HPP


#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/generics/any.hpp>



namespace graphlab {


  /**
   * A constant shared data entry.  
   * 
   * glshared_const<datatype> variable; 
   *
   * The glshared_const<datatype> container is used to store constants
   * that are to be accessed by update functions.  The glshared_const
   * serves the purpose of automatically managing the physical
   * placement of constants in the NUMA and cluster setting but should
   * be used in all settings for portability.
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
    //! Construct initial shared pointers
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
  };



}
#endif
