#ifndef GRAPHLAB_UNSUPPORTED_SERIALIZE_HPP
#define GRAPHLAB_UNSUPPORTED_SERIALIZE_HPP

#include <serialization/iarchive.hpp>
#include <serialization/oarchive.hpp>
#include <logger/logger.hpp>

namespace graphlab {

  struct unsupported_serialize {
    void save(oarchive& archive) const {      
      ASSERT_MSG(false, "trying to serialize an unserializable object");
    }
    void load(iarchive& archive) {
      ASSERT_MSG(false, "trying to deserialize an unserializable object");
    }
  }; // end of struct
};

#define GRAPHLAB_UNSERIALIZABLE(TNAME)  \
  graphlab::oarchive& operator<<(graphlab::oarchive&,   \
                                const TNAME &mat) {ASSERT_MSG(false,"trying to serialize an unserializable object");} \
  graphlab::iarchive& operator>>(graphlab::iarchive&,  \
                                 TNAME &vec) {ASSERT_MSG(false,"trying to deserialize an unserializable object");}



#endif
