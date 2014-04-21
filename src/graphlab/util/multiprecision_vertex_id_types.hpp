#include <boost/multiprecision/cpp_int.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

// Integer of 128 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int128_t , x)
  std::stringstream ss;
  ss << "0x" << std::hex << x;
  std::string s;
  ss >> s;
  oarc << s;
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int128_t , x)
  std::string s;
  iarc >> s;
  x.assign(s);
END_OUT_OF_PLACE_LOAD()

namespace boost {
  namespace multiprecision {
    inline size_t hash_value (const boost::multiprecision::int128_t &x)
    {
      // not a very good hash function, but I just want to get it working first!
      return static_cast<size_t>(x);
    }
  }
}

// Integer of 256 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int256_t , x)
  std::stringstream ss;
  ss << "0x" << std::hex << x;
  std::string s;
  ss >> s;
  oarc << s;
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int256_t , x)
  std::string s;
  iarc >> s;
  x.assign(s);
END_OUT_OF_PLACE_LOAD()

namespace boost {
  namespace multiprecision {
    inline size_t hash_value (const boost::multiprecision::int256_t &x)
    {
      // not a very good hash function, but I just want to get it working first!
      return static_cast<size_t>(x);
    }
  }
}

// Integer of 512 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int512_t , x)
  std::stringstream ss;
  ss << "0x" << std::hex << x;
  std::string s;
  ss >> s;
  oarc << s;
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int512_t , x)
  std::string s;
  iarc >> s;
  x.assign(s);
END_OUT_OF_PLACE_LOAD()

namespace boost {
  namespace multiprecision {
    inline size_t hash_value (const boost::multiprecision::int512_t &x)
    {
      // not a very good hash function, but I just want to get it working first!
      return static_cast<size_t>(x);
    }
  }
}

// Integer of 1024 bits.

BEGIN_OUT_OF_PLACE_SAVE(oarc, boost::multiprecision::int1024_t , x)
  std::stringstream ss;
  ss << "0x" << std::hex << x;
  std::string s;
  ss >> s;
  oarc << s;
END_OUT_OF_PLACE_SAVE()

BEGIN_OUT_OF_PLACE_LOAD(iarc, boost::multiprecision::int1024_t , x)
  std::string s;
  iarc >> s;
  x.assign(s);
END_OUT_OF_PLACE_LOAD()

namespace boost {
  namespace multiprecision {
    inline size_t hash_value (const boost::multiprecision::int1024_t &x)
    {
      // not a very good hash function, but I just want to get it working first!
      return static_cast<size_t>(x);
    }
  }
}


