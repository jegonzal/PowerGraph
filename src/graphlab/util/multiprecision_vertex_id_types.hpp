#include <boost/multiprecision/cpp_int.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

#define MULTIPRECISION_VERTEX_ID_SAVE std::stringstream ss;\
  ss << x;\
  std::string s;\
  ss >> s;\
  oarc << s;

#define MULTIPRECISION_VERTEX_ID_LOAD std::string s;\
  iarc >> s;\
  x = 0;\
  x.assign(s);

// Integer of 128 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int128_t , x)
  MULTIPRECISION_VERTEX_ID_SAVE
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int128_t , x)
  MULTIPRECISION_VERTEX_ID_LOAD
END_OUT_OF_PLACE_LOAD()

// Integer of 256 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int256_t , x)
  MULTIPRECISION_VERTEX_ID_SAVE
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int256_t , x)
  MULTIPRECISION_VERTEX_ID_LOAD
END_OUT_OF_PLACE_LOAD()

// Integer of 512 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int512_t , x)
  MULTIPRECISION_VERTEX_ID_SAVE
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int512_t , x)
  MULTIPRECISION_VERTEX_ID_LOAD
END_OUT_OF_PLACE_LOAD()

// Integer of 1024 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int1024_t , x)
  MULTIPRECISION_VERTEX_ID_SAVE
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int1024_t , x)
  MULTIPRECISION_VERTEX_ID_LOAD
END_OUT_OF_PLACE_LOAD()

namespace boost {
  namespace multiprecision {
    template <typename T>
    inline size_t hash_value (const T &x)
    {
      // not a very good hash function, but I just want to get it working first!
      return static_cast<int>(x);
    }
  }
}



