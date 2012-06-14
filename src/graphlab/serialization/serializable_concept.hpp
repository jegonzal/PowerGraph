#ifndef GRAPHLAB_SERIALIZABLE
#define GRAPHLAB_SERIALIZABLE
#include <boost/concept/assert.hpp>
#include <boost/concept/requires.hpp>
#include <boost/concept_check.hpp>
#include <sstream>
#include <graphlab/serialization/serialize.hpp>
namespace graphlab {

  /**
   * \brief Concept checks if a type T is \ref serializable.
   */
  template <typename T>
  class Serializable : boost::Assignable<T>, boost::DefaultConstructible<T> {
   public:
    BOOST_CONCEPT_USAGE(Serializable) {
      std::stringstream strm;
      oarchive oarc(strm);
      iarchive iarc(strm);
      T t;
      oarc << t;
      iarc >> t;
    }
  };

} // namespace graphlab
#endif
