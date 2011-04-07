#ifndef SERIALIZABLE_POD_HPP
#define SERIALIZABLE_POD_HPP

#include <graphlab/serialization/is_pod.hpp>

#define SERIALIZABLE_POD(tname)                   \
namespace graphlab {                              \
    template <>                                   \
    struct gl_is_pod<tname> {                     \
      BOOST_STATIC_CONSTANT(bool, value = true);  \
    };                                            \
}

#endif

