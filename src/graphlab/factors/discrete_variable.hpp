#ifndef GRAPHLAB_VARIABLE_HPP
#define GRAPHLAB_VARIABLE_HPP
/**
 * This file contains the definitions of a simple discrete variable.
 *
 *  \author Joseph Gonzalez
 */

#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {

  /** represents a discrete variable */
  class discrete_variable {
  public:
    //! The type used to store the variable id
    typedef uint32_t id_type;
    //! the type used to index the variable assignments
    typedef uint32_t index_type;
    /** construct a discrte variable with a given id and number of
        assignments */
    inline discrete_variable(id_type id = 0, index_type nasgs = 0) : 
      id_(id), nasgs_(nasgs) { }
    //! get the variable id
    inline id_type& id() { return id_; }
    //! get the variable id
    inline const id_type& id() const { return id_; }
    //! get the number of assignments the variable can take
    inline index_type& size() { return nasgs_; }   
    //! get the number of assignments the variable can take
    inline const index_type& size() const { return nasgs_; }
    //! Compare two variables
    inline bool operator<(const discrete_variable& other) const { 
      return id_ < other.id_; 
    }
    //! test equality between two varaibles
    inline bool operator==(const discrete_variable& other) const {
      return id_ == other.id_;     
    }
    //! Test inequality between two variables
    inline bool operator!=(const discrete_variable& other) const { 
      return id_ != other.id_; 
    }
    //! load the variable from an archive
    void load(iarchive& arc) { arc >> id_ >> nasgs_; }
    //! save the variable to an archive
    void save(oarchive& arc) const { arc << id_ << nasgs_; }
    
  private:
    //! The variable id
    id_type id_;
    //! The number of assignments the variable takes
    index_type nasgs_;

  };
}; // end of namespace graphlab


std::ostream& operator<<(std::ostream& out, 
                         const graphlab::discrete_variable& var);


#endif
